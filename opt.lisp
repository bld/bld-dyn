;; Optimized BLD-DYN functions
(in-package :bld-dyn)

;; Ensure rv2dudt returns an RE3 for 3D problem
(defgamethod e3 rv2dudt bld-e3::*spec* ve3 re3 ve3)

;; Optimize KSHEOM sub-functions
#|(defgamethod e3 u bld-e3::*spec* re3 re3 number number)
(defgamethod e3 duds bld-e3::*spec* re3 re3 number number)
(defgamethod e3 dalphads bld-e3::*spec* ve3 number number)
(defgamethod e3 dbetads bld-e3::*spec* ve3 number number)
(defgamethod e3 deds bld-e3::*spec* ve3 re3 ve3 re3)
(defgamethod e3 dtmds bld-e3::*spec* re3)
|#