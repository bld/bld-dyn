(defpackage :bld.dyn.system
  (:use :asdf :cl))
(in-package :bld.dyn.system)
(defsystem :bld-dyn
  :name "bld-dyn"
  :author "Ben Diedrich"
  :version "0.0.1"
  :maintainer "Ben Diedrich"
  :license "MIT"
  :description "Dynamics library employing geometric algebra"
  :components
  ((:file "package")
   (:file "dyn" :depends-on ("package")))
  :depends-on ("bld-ga" "bld-e2" "bld-e3" "bld-ode" "bld-utils" "bld-gen"))
