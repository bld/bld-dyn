(asdf:defsystem :bld-dyn
  :author "Ben Diedrich"
  :license "MIT"
  :description "Dynamics library employing geometric algebra"
  :components
  ((:file "package")
   (:file "dyn" :depends-on ("package")))
  :depends-on ("bld-ga" "bld-e2" "bld-e3" "bld-ode" "bld-utils" "bld-gen" "lol"))
