(ns germs.physics
  (:refer-clojure :exclude [force])
  (:require
    [germs.record-util :refer [defrecord+]]))

(defprotocol ISimulated
  (calculate-derivatives [self dt])
  (integrate-euler [self velocities dt])
  (resolve-collisions [self]))

(defrecord+ PhysicsObject [mass position velocity])

(defrecord+ Derivatives [acceleration velocity])

(def ^:private gravity [0.0 1200.0])
(def ^:private air-friction 7.0)

(defn dot-product [a b]
  (apply + (map * a b)))

(defn magnitude [a]
  (->> a
    (map #(Math/pow % 2))
    (apply +)
    (Math/sqrt)))

(defn normalize [a]
  (map #(/ % (magnitude a)) a))

(defn vector-angle [a b]
  (->>
    (map normalize [a b])
    (apply dot-product)
    (Math/acos)))

(defn get-centroid [objects]
  (->> objects
    (reduce (fn [c o]
              (let [p (get-in o [:physics :position])]
                (map + c p)))
            [0.0 0.0])
    (map #(/ % (count objects)))))

(defn resolve-gravity [{{:keys [mass]} :physics}]
  (map #(* % mass) gravity))

; TODO Is this even working?
(defn resolve-tortional-spring
  [{{left :position} :physics}
   {{right :position} :physics}
   {{point :position} :physics
    :keys [rest-angle tortional-k]}]
  (let [a (map - left point)
        b (map - right point)
        angle (vector-angle a b)
        delta (- rest-angle angle)
        bisector (normalize (map + a b))
        force (map #(* % delta tortional-k) bisector)]
    (map #(* % delta tortional-k) bisector)))

(defn- spring-force [a b neutral k]
  (let [displacement (map - b a)
        distance (magnitude displacement)
        delta (- distance neutral)
        direction (map #(/ % distance) displacement)]
    (map #(* % delta k) direction)))

(defn resolve-linear-spring
  [{{neighbor :position} :physics}
   {{point :position} :physics
    :keys [rest-distance linear-k]}]
  (spring-force point neighbor rest-distance linear-k))

(defn resolve-centroid-spring
  [centroid
   {{point :position} :physics
    :keys [rest-radius radial-k]}]
  (spring-force point centroid rest-radius radial-k))

(defn resolve-air-friction
  [{{:keys [velocity]} :physics}
   force]
  (let [speed (magnitude velocity)]
    (if (> speed 0.0)
      (let [direction (map #(- (/ % speed)) velocity)
            friction-force (map #(* % air-friction) direction)]
        (map
          (fn [f ff]
            (if (< (Math/abs f) (Math/abs ff))
              0.0
              ff))
          force
          friction-force))
      [0.0 0.0])))

(defn resolve-acceleration
  [{{:keys [mass]} :physics}
   force]
  (map #(/ % mass) force))

(defn resolve-velocity
  [{{:keys [velocity]} :physics}
   acceleration
   dt]
  (map #(+ %1 (* %2 dt)) velocity acceleration))

(defn resolve-position
  [{{:keys [position]} :physics}
   velocity
   dt]
  (map #(+ %1 (* %2 dt)) position velocity))

(defn collide [point component condition extreme]
  (let [position-path [:physics :position component]]
    (if (condition (get-in point position-path) extreme)
      (-> point
        (assoc-in position-path extreme)
        (update-in [:physics :velocity component] -)
        ; TODO: Make the collision damping configurable
        (update-in [:physics :velocity] (fn [v] (mapv #(* % 0.3) v))))
      point)))

(defn- rk4-weight [derivatives property]
  (let [[d0 d1 d2 d3] (map property derivatives)]
    (map #(* % (/ 1.0 6.0))
      (map +
        d0
        (map #(* % 2.0)
          (map + d1 d2))
        d3))))

(defn- rk4-weights [ds0 ds1 ds2 ds3]
  (map
    (fn [& derivatives]
      (Derivatives.
        (rk4-weight derivatives :acceleration)
        (rk4-weight derivatives :velocity)))
    ds0 ds1 ds2 ds3))

(defn simulate [dt objects]
  ; If objects affect one another's acceleration, they need to be integrated
  ; in lockstep instead of independently.
  (let [[h0 h1 h2 h3] (map #(* % dt) [0.0 0.5 0.5 1.0])]
    (pmap
      (fn [object]
        (let [integrator (fn [ds h] (integrate-euler object (map :velocity ds) h))
              ds0 (calculate-derivatives object h0)
              ds1 (calculate-derivatives (integrator ds0 h1) h1)
              ds2 (calculate-derivatives (integrator ds1 h2) h2)
              ds3 (calculate-derivatives (integrator ds2 h3) h3)
              ds-rk4 (rk4-weights ds0 ds1 ds2 ds3)
              integrated (integrator ds-rk4 dt)]
           (resolve-collisions integrated)))
      objects)))
