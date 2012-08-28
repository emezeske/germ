(ns germs.core
  (:refer-clojure :exclude [force])
  (:require
    [clojure.pprint :refer [pprint]]
    [germs.record-util :refer [defrecord+]]
    [quil
      [applet :as applet]
      [core :as quil]]))

(def ^:private window-size [500 500])
(def ^:private window-center (mapv #(/ % 2) window-size))
(def ^:private gravity [0.0 1500.0])
(def ^:private framerate 75)
(def ^:private air-friction 7.0)

(defrecord+ PhysicsObject [mass position velocity])

(defrecord+ GermSkinPoint
  [^PhysicsObject physics
   rest-angle
   rest-radius
   rest-distance
   tortional-k
   radial-k
   linear-k])

; TODO: Experiment with making the edge low mass and actually
;       simulating a high-mass core. This seems MUCH better.

; TODO: It seems like the angle springs are not necessary.

(defn- make-points [n]
  (for [i (range n)]
    (let [theta (* (/ i n) 2.0 Math/PI)
          radius 100
          position [(Math/cos theta) (Math/sin theta)]
          position (mapv #(* % radius) position)
          position (mapv + position window-center)
          neighbor-angle (/ (* 2.0 Math/PI ) n)
          rest-angle (- Math/PI neighbor-angle)
          rest-distance (* 1.0 radius (Math/sin neighbor-angle))]
      (make-germ-skin-point
        {:physics
          (make-physics-object
            {:mass 0.02
             :position position
             :velocity [0.0 0.0]})
         :rest-angle rest-angle
         :rest-radius radius
         :rest-distance rest-distance
         :tortional-k 60.0
         :radial-k 60.0
         :linear-k 60.0}))))

(def ^:private points (atom (make-points 11)))

(defn- dot-product [a b]
  (apply + (mapv * a b)))

(defn- magnitude [a]
  (->> a
    (mapv #(Math/pow % 2))
    (apply +)
    (Math/sqrt)))

(defn- normalize [a]
  (mapv #(/ % (magnitude a)) a))

(defn- vector-angle [a b]
  (->>
    (mapv normalize [a b])
    (apply dot-product)
    (Math/acos)))

(defn- get-centroid [objects]
  (->> objects
    (reduce (fn [c o]
              (let [p (get-in o [:physics :position])]
                (mapv + c p)))
            [0.0 0.0])
    (mapv #(/ % (count objects)))))

(defn- resolve-gravity [{{:keys [mass]} :physics}]
  (mapv #(* % mass) gravity))

; TODO Is this even working?
(defn- resolve-tortional-spring
  [{{left :position} :physics}
   {{right :position} :physics}
   {{point :position} :physics
    :keys [rest-angle tortional-k]}]
  (let [a (mapv - left point)
        b (mapv - right point)
        angle (vector-angle a b)
        delta (- rest-angle angle)
        bisector (normalize (mapv + a b))
        force (mapv #(* % delta tortional-k) bisector)]
    (mapv #(* % delta tortional-k) bisector)))

(defn- spring-force [a b neutral k]
  (let [displacement (mapv - b a)
        distance (magnitude displacement)
        delta (- distance neutral)
        direction (mapv #(/ % distance) displacement)]
    (mapv #(* % delta k) direction)))

(defn- resolve-linear-spring
  [{{neighbor :position} :physics}
   {{point :position} :physics
    :keys [rest-distance linear-k]}]
  (spring-force point neighbor rest-distance linear-k))

(defn- resolve-centroid-spring
  [centroid
   {{point :position} :physics
    :keys [rest-radius radial-k]}]
  (spring-force point centroid rest-radius radial-k))

(defn- resolve-force [point left right centroid]
  (reduce (partial map +)
    ((juxt
      resolve-gravity
      (partial resolve-tortional-spring left right)
      (partial resolve-linear-spring left)
      (partial resolve-linear-spring right)
      (partial resolve-centroid-spring centroid))
     point)))

(defn- resolve-air-friction
  [{{:keys [velocity]} :physics}
   force]
  (let [speed (magnitude velocity)]
    (if (> speed 0.0)
      (let [direction (mapv #(- (/ % speed)) velocity)
            friction-force (mapv #(* % air-friction) direction)]
        (mapv
          (fn [f ff]
            (if (< (Math/abs f) (Math/abs ff))
              0.0
              ff))
          force
          friction-force))
      [0.0 0.0])))

(defn- resolve-acceleration
  [{{:keys [mass]} :physics}
   force]
  (mapv #(/ % mass) force))

(defn- resolve-velocity
  [{{:keys [velocity]} :physics}
   acceleration
   dt]
  (mapv #(+ %1 (* %2 dt)) velocity acceleration))

(defn- resolve-position
  [{{:keys [position]} :physics}
   velocity
   dt]
  (mapv #(+ %1 (* %2 dt)) position velocity))

(defn- collide [point component condition extreme]
  (let [position-path [:physics :position component]]
    (if (condition (get-in point position-path) extreme)
      (-> point
        (assoc-in position-path extreme)
        (update-in [:physics :velocity component] -)
        ; TODO: Make the collision damping configurable
        (update-in [:physics :velocity] (fn [v] (mapv #(* % 0.3) v))))
      point)))

(defn- resolve-collisions [point]
  (-> point
    (collide 0 < 0)
    (collide 0 > 500)
    (collide 1 < 0)
    (collide 1 > 500)))

; FIXME: For debug:
(pprint @points)

(defn- simulate-points [dt points]
  (let [points-wrapped (concat [(last points)] points [(first points)])
        point-groups (partition 3 1 points-wrapped)
        centroid (get-centroid points)
        ; FIXME Testing
        ;centroid (mapv + centroid [0.0 (* dt 1000.0)])
        ]
    (for [[left point right] point-groups]
      (let [force (resolve-force point left right centroid)
            force (map + force (resolve-air-friction point force))
            acceleration (resolve-acceleration point force)
            velocity (resolve-velocity point acceleration dt)
            position (resolve-position point velocity dt)]
        (-> point
          (assoc-in [:physics :velocity] velocity)
          (assoc-in [:physics :position] position)
          (resolve-collisions))))))

(defn- draw []
  (quil/background 180)
  (quil/stroke-weight 1)
  (quil/fill 100 255)
  ;(doseq [{{[x y] :position} :physics} @points]
  ;  (quil/ellipse x y 10 10))
  (quil/begin-shape)
  (let [head [(last @points)]
        tail [(first @points) (second @points)]
        points-with-ends (concat head @points tail)]
    (doseq [{{p :position} :physics} points-with-ends]
      (apply quil/curve-vertex p)))
  (quil/end-shape)
  (swap! points
    (fn [points]
      (let [steps 10]
        (nth
          (iterate
            (partial simulate-points (/ 1.0 framerate steps))
            points)
          steps)))))

(defn- setup []
  (quil/frame-rate framerate)
  (quil/smooth))

; Don't use quil/defsketch because it e.g. doesn't let you pass in a Var
; for the :size parameter.
(defonce sketch
  (applet/applet
    :title "Germ"
    :setup #'setup
    :draw #'draw
    :size window-size))
