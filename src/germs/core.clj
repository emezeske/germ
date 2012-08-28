(ns germs.core
  (:refer-clojure :exclude [force])
  (:require
    [quil.core :as quil]))

(def ^:private window-size [500 500])
(def ^:private window-center (mapv #(/ % 2) window-size))
(def ^:private gravity [0.0 1500.0])
(def ^:private framerate 75)
(def ^:private air-friction 7.0)

(defrecord Absolutes [mass position])

(defn- make-absolutes [init]
  (merge
    (Absolutes. nil nil)
    init))

(defrecord Derivatives [velocity force])

(defn- make-derivatives [init]
  (merge
    (Derivatives. nil nil)
    init))

(defrecord PhysicsObject [^Absolutes absolutes ^Derivatives derivatives])

(defn- make-physics-object [init]
  (merge
    (PhysicsObject. nil nil)
    init))

(defrecord GermSkinPoint
  [^PhysicsObject physics
   rest-angle
   rest-radius
   rest-distance
   tortional-k
   radial-k
   linear-k])

(defn- make-germ-skin-point [init]
  (merge
    (GermSkinPoint. nil nil nil nil nil nil nil)
    init))

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
            {:absolutes
              (make-absolutes {:mass 0.02 :position position})
             :derivatives
              (make-derivatives {:velocity [0.0 0.0] :force [0.0 0.0]})})
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
              (let [p (get-in o [:physics :absolutes :position])]
                (mapv + c p)))
            [0.0 0.0])
    (mapv #(/ % (count objects)))))

(defn- reset-force [object]
  (assoc-in object [:physics :derivatives :force] [0.0 0.0]))

(defn- add-force [point force]
  (update-in point [:physics :derivatives :force] #(mapv + % force)))

(defn- apply-gravity [point]
  (let [mass (get-in point [:physics :absolutes :mass])]
    (add-force point (mapv #(* % mass) gravity))))

; TODO Is this even working?
(defn- apply-tortional-spring
  [{{{center :position} :absolutes} :physics
    :keys [rest-angle tortional-k]
    :as point}
   {{{neighbor-left :position} :absolutes} :physics}
   {{{neighbor-right :position} :absolutes} :physics}]
  (let [a (mapv - neighbor-left center)
        b (mapv - neighbor-right center)
        angle (vector-angle a b)
        delta (- rest-angle angle)
        bisector (normalize (mapv + a b))
        force (mapv #(* % delta tortional-k) bisector)]
    (add-force point force)))

(defn- spring-force [a b neutral k]
  (let [displacement (mapv - b a)
        distance (magnitude displacement)
        delta (- distance neutral)
        direction (mapv #(/ % distance) displacement)]
    (mapv #(* % delta k) direction)))

(defn- apply-linear-spring
  [{{{pa :position} :absolutes} :physics
    :keys [rest-distance linear-k]
    :as point}
   {{{pb :position} :absolutes} :physics}]
  (add-force point (spring-force pa pb rest-distance linear-k)))

(defn- apply-centroid-spring
  [{{{pa :position} :absolutes} :physics
    :keys [rest-radius radial-k]
    :as point}
   centroid]
  (add-force point (spring-force pa centroid rest-radius radial-k)))

(defn- integrate [point field del dt]
  (let [delta (mapv #(* % dt) del)]
    (update-in point field #(mapv + % delta))))

(defn- apply-acceleration
  [{{{:keys [mass]} :absolutes
     {:keys [force]} :derivatives} :physics
    :as point}
   dt]
  (let [acceleration (mapv #(/ % mass) force)]
    (integrate point [:physics :derivatives :velocity] acceleration dt)))

(defn- apply-velocity
  [{{{:keys [velocity]} :derivatives} :physics
    :as point}
   dt]
  (integrate point [:physics :absolutes :position] velocity dt))

(defn- apply-air-friction
  [{{{:keys [velocity]} :derivatives} :physics
    :as point}]
  (let [speed (magnitude velocity)]
    (if (> speed 0.0)
      (let [direction (mapv #(- (/ % speed)) velocity)
            friction-force (mapv #(* % air-friction) direction)]
        (update-in point [:physics :derivatives :force]
          (fn [force]
            (mapv
              (fn [f ff]
                (if (< (Math/abs f) (Math/abs ff))
                  0
                  (+ f ff)))
              force
              friction-force))))
      point)))

(defn- collide [point component condition extreme]
  (let [position-path [:physics :absolutes :position component]]
    (if (condition (get-in point position-path) extreme)
      (-> point
        (assoc-in position-path extreme)
        (update-in [:physics :derivatives :velocity component] -)
        ; TODO: Make the collision damping configurable
        (update-in [:physics :derivatives :velocity] (fn [v] (mapv #(* % 0.3) v))))
      point)))

(defn- apply-collisions [point]
  (-> point
    (collide 0 < 0)
    (collide 0 > 500)
    (collide 1 < 0)
    (collide 1 > 500)))

(defn- debug [point]
  (when (:debug? point)
   ;(println point)
   ;(flush)
    )
  point)

(defn- simulate-points [dt points]
  (let [points-wrapped (concat [(last points)] points [(first points)])
        point-groups (partition 3 1 points-wrapped)
        centroid (get-centroid points)
        ; FIXME Testing
    ;    centroid (mapv + centroid [0.0 (* dt 1000.0)])
        ]
    (for [[neighbor-left point neighbor-right] point-groups]
      (-> point
        (reset-force)
        (apply-gravity)
        (apply-tortional-spring neighbor-left neighbor-right)
        (apply-linear-spring neighbor-right)
        (apply-linear-spring neighbor-left)
        (apply-centroid-spring centroid)
        (apply-air-friction)
        (apply-acceleration dt)
        (apply-velocity dt)
        (apply-collisions)
        (debug)))))

(defn- draw []
  (quil/background 180)
  (quil/stroke-weight 1)
  (quil/fill 100 255)
  ;(doseq [{{{[x y] :position} :absolutes} :physics} @points]
  ;  (quil/ellipse x y 10 10))
  (quil/begin-shape)
  (let [head [(last @points)]
        tail [(first @points) (second @points)]
        points-with-ends (concat head @points tail)]
    (doseq [{{{p :position} :absolutes} :physics} points-with-ends]
      (apply quil/curve-vertex p)))
  (quil/end-shape)
  (swap! points
    (fn [points]
      (let [steps 5]
        (nth
          (iterate
            (partial simulate-points (/ 1.0 framerate steps))
            points)
          steps)))))

(defn- setup []
  (quil/frame-rate framerate)
  (quil/smooth))

(defonce sketch
  (quil/defsketch germs
    :title "Germ"
    :setup setup
    :draw draw
    ; FIXME: Seriously, size can't be a Var...?
    :size [500 500]
    :keep-on-top true))
