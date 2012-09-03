(ns germs.graphics
  (:require
    [quil.core :as quil]))

(defprotocol IDrawable
  (draw [self]))

(def window-size [500 500])
(def framerate 75)

(defn setup []
  (quil/frame-rate framerate)
  (quil/smooth))

(defn draw [germs]
  (quil/background 180)
  (quil/stroke-weight 1)
  (quil/fill 100 255)
  (doseq [germ germs]
    (draw germ)))
