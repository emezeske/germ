(ns germs.record-util)

(defn- camel-to-dash [camel]
  (apply str
    (map
      #(if (Character/isUpperCase %)
        (str "-" (Character/toLowerCase %))
        %)
      camel)))

(defmacro defrecord+
  "Just like defrecord, except it generates a make-* constructor function
   that accepts a map of members."
  [record-name fields & args]
  (let [make (symbol (str "make" (camel-to-dash (name record-name))))
        nils (repeat (count fields) nil)]
    `(do
       (defrecord ~record-name ~fields ~@args)
       (defn- ~make [init#]
         (merge
           (new ~record-name ~@nils)
           init#)))))
