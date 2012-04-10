(asdf:defsystem :bld-orbit
  :author "Ben Diedrich"
  :version "0.0.1"
  :license "MIT"
  :description "Dynamics library employing geometric algebra"
  :components
  ((:file "package-orbit")
   (:file "orbit" :depends-on ("package-orbit")))
  :depends-on ("bld-ga" "bld-e2" "bld-e3" "bld-ode" "bld-utils" "bld-gen"))
