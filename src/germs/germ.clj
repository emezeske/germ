(ns germs.germ
  (:refer-clojure :exclude [force])
  (:require
    [germs.graphics :as graphics]
    [germs.physics :as physics]
    [germs.record-util :refer [defrecord+]]
    [quil.core :as quil])
  (:import [germs.physics Derivatives]))

(defrecord+ GermSkinPoint
  [physics
   rest-angle
   rest-radius
   rest-distance
   tortional-k
   radial-k
   linear-k])

; TODO: Experiment with making the edge low mass and actually
;       simulating a high-mass core. This seems MUCH better.

; TODO: It seems like the angle springs are not necessary.

(defn- resolve-force [point left right centroid]
  (reduce (partial map +)
    ((juxt
      physics/resolve-gravity
      (partial physics/resolve-tortional-spring left right)
      (partial physics/resolve-linear-spring left)
      (partial physics/resolve-linear-spring right)
      (partial physics/resolve-centroid-spring centroid))
     point)))

(defn- calculate-derivative [point left right centroid dt]
  (let [force (resolve-force point left right centroid)
        force (map + force (physics/resolve-air-friction point force))
        acceleration (physics/resolve-acceleration point force)
        velocity (physics/resolve-velocity point acceleration dt)]
    (Derivatives. acceleration velocity)))

(defrecord+ Germ [points]
  physics/ISimulated
  (calculate-derivatives [self dt]
    (let [centroid (physics/get-centroid points)
          npoints (count points)]
      (for [i (range npoints)]
        (let [before (mod (dec i) npoints)
              after (mod (inc i) npoints)]
          (calculate-derivative
            (nth points i) (nth points before) (nth points after) centroid dt)))))

  (integrate-euler [self velocities dt]
    (Germ.
      (map
        (fn [point velocity]
          (let [position (physics/resolve-position point velocity dt)]
            (-> point
              (assoc-in [:physics :velocity] (vec velocity))
              (assoc-in [:physics :position] (vec position)))))
         points velocities)))

  (resolve-collisions [self]
    (Germ.
      (for [point points]
        (-> point
          (physics/collide 0 < 0)
          (physics/collide 0 > (first graphics/window-size))
          (physics/collide 1 < 0)
          (physics/collide 1 > (second graphics/window-size))))))

  graphics/IDrawable
  (draw [self]
    ;(doseq [{{[x y] :position} :physics} points]
    ;  (quil/ellipse x y 10 10))
    (quil/begin-shape)
    (let [npoints (count points)
          indices (concat [(dec npoints)] (range 0 npoints) [0 1])]
      (doseq [i indices]
        (apply quil/curve-vertex (:position (:physics (nth points i))))))
    (quil/end-shape)))

(defn make-germ [point-count center]
  (Germ.
    (for [i (range point-count)]
      (let [theta (* (/ i point-count) 2.0 Math/PI)
            radius 100
            position [(Math/cos theta) (Math/sin theta)]
            position (map #(* % radius) position)
            position (map + position center)
            neighbor-angle (/ (* 2.0 Math/PI ) point-count)
            rest-angle (- Math/PI neighbor-angle)
            rest-distance (* 1.0 radius (Math/sin neighbor-angle))]
        (make-germ-skin-point
          {:physics
            (physics/make-physics-object
              {:mass 0.02
               :position (vec position)
               :velocity [0.0 0.0]})
           :rest-angle rest-angle
           :rest-radius radius
           :rest-distance rest-distance
           :tortional-k 10.0
           :radial-k 20.0
           :linear-k 20.0})))))
