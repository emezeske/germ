(ns germs.core
  (:require
    [clojure.pprint :refer [pprint]]
    [germs.germ :as germ]
    [germs.graphics :as graphics]
    [germs.physics :as physics]
    [germs.record-util :refer [defrecord+]]
    [quil.applet :as applet]))

(def ^:private germs
  (atom
    [(germ/make-germ 11 [125 110])
     (germ/make-germ 11 [375 110])]))

; FIXME: For debug:

;(pprint @germs)

;(pprint (collide (first (.-points (first @germs))) 0 < 0))

;(time
;  (dotimes [_ 1000]
;    (doall (physics/simulate (/ 1.0 framerate) @germs))))

(defn- step []
  (swap! germs
    (fn [germs]
      (let [steps 3]
        (nth
          (iterate
            (partial physics/simulate (/ 1.0 graphics/framerate steps))
            germs)
          steps))))
  (graphics/draw @germs))

; Don't use quil/defsketch because it e.g. doesn't let you pass in a Var
; for the :size parameter.
(defonce sketch
  (applet/applet
    :title "Germ"
    :setup #'graphics/setup
    :draw #'step
    :size graphics/window-size))
