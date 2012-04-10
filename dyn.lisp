(in-package :bld-dyn)

;; Define infinity norm method for BLD-ODE Runge-Kutta method
(defmethod norminfx ((x g))
  (norminf x))

;; Solar sail force functions
(defun make-idealsail (&key (area 10000d0) (w 1368d0) (stype 'flat) (c 2.99792458d8) (du 1.49598d11))
  "Make ideal solar sail force function of normal vector and position vector wrt the sun.
AREA: sail area
W: solar intensity (power/area) at distance unit
STYPE: 'flat or 'compound
C: speed of light
DU: distance unit where intensity defined"
  (assert (or (eql stype 'flat) (eql stype 'compound)) (stype) "Sail type must be one of 'FLAT or 'COMPOUND")
  (let ((constant (* 2 area (/ w c) (expt du 2)))
	(cosfn (case stype
		 ('flat (lambda (n rs) (expt (scalar (*i n (unitg rs))) 2d0)))
		 ('compound (lambda (n rs) (scalar (*i n (unitg rs))))))))
    (lambda (n rs)
      (* constant
	 (funcall cosfn n rs)
	 (/ (norme2 rs))
	 n))))

;; 3D inertia linear function from inertia matrix
(defun make-inertia (tensor basis)
  "Make 3D inertia function of bivectors given inertia tensor and basis"
  (lambda (bv)
    (let ((result (bve3)))
      (dotimes (j 3)
	(dotimes (k 3)
	  (setf result (- result
			  (* (aref tensor j k)
			     (*s bv (dual (nth k basis)))
			     (dual (nth j basis)))))))
      result)))

(defun make-inertiav (tensor basis)
  "Make an inertia linear function of a vector given the tensor (3x3 array) and basis (list of 3 E3 vectors)"
  (lambda (v)
    (let ((result (ve3)))
      (dotimes (j 3)
	(dotimes (k 3)
	  (setf result (+ result
			  (* (aref tensor j k) ; Ijk
			     (*s v (nth k basis)) ; vk = v . basisk
			     (nth j basis)))))) ; basisj
      result)))

;; Euler equations of motion
