(in-package :bld-orbit)

;; Define infinity norm method for BLD-ODE Runge-Kutta method
(defmethod norminfx ((x g))
  (norminf x))

;; Sail force functions
(defvar *lightness* 0.1 "Sail lightness number")
(defvar *mu* 1 "Gravitational parameter")
(defmethod sail-normal-force ((rv g) mu lightness)
  (* (unitg rv) (/ (* lightness mu) (norme2 rv))))

;; Cartesian equations of motion
(defvar *cart-forcefun* #'(lambda (tm x) (ve3)))
(defvar *cart-forcefun-sail-normal*
  #'(lambda (tm x)
      (sail-normal-force (gethash :r x) *mu* *lightness*)))
(defmethod dvdt ((r g))
  "Cartesian gravitational acceleration"
  (* r (/ -1 (expt (norme r) 3))))

(defun carteom (tm x)
  "Cartesian orbital equations of motion"
  (with-keys (r v) x
    (make-hash
     :r v
     :v (+ (dvdt r) (funcall *cart-forcefun* tm x)))))

;; KS equations of motion

(defvar *ksh-forcefun* #'(lambda (s x) (ve3)) "KSH force function")
(defvar *ksh-sigma0* (ve3 :c1 1) "KSH reference unit position vector")

(defun w0 (e)
  "Average orbit angular velocity given energy"
  (sqrt (- (/ e 2))))
(defmethod u ((alpha g) (beta g) w0 s)
  "Spinor"
  (+ (* alpha (cos (* w0 s)))
     (* beta (sin (* w0 s)))))
(defmethod alpha ((u g) (duds g) w0 s)
  "Alpha (U0) given U, DUDS, W0, and S"
  (- (* u (cos (* w0 s)))
     (* duds (/ (sin (* w0 s))
		w0))))
(defmethod beta ((u g) (duds g) w0 s)
  "Beta (dU/ds/w0) given U, DUDS, w0, and s"
  (+ (* u (sin (* w0 s)))
     (* duds (/ (cos (* w0 s))
		w0))))
(defmethod duds ((beta g) w0)
  "s-derivative of spinor"
  (* beta w0))
(defmethod dalphads ((ff g) w0 s)
  "s-derivative of alpha"
  (* ff (- (/ (sin (* w0 s)) w0))))
(defmethod dbetads ((ff g) w0 s)
  "s-derivative of beta"
  (* ff (/ (cos (* w0 s)) w0)))
(defmethod deds ((f g) (duds g) (sigma g) (u g))
  "s-derivative of energy"
  (scalar (*i f (*g3 duds sigma (revg u)))))
(defmethod dtmds ((u g))
  "s-derivative of time"
  (norme2 u))
(defun ksheom (s x)
  "KS equations of motion in Hestenes GA form. 
Orbit elements are: alpha (U0), beta (dU0/ds/w0), e (orbit energy), tm (time) 
Expects a variable *data* containing sigma0 (initial orbit frame vector) and forcefun (force function of s, x, and *data*)"
  (with-keys (alpha beta e tm) x
    (let* ((w0 (w0 e))
	   (u (u alpha beta w0 s))
	   (duds (duds beta w0))
	   (sigma (spin *ksh-sigma0* alpha))
	   (r (spin sigma u))
	   (f (funcall *ksh-forcefun* s x))
	   (ff (*g3 f r u)))
      (make-hash
       :alpha (dalphads ff w0 s)
       :beta (dbetads ff w0 s)
       :e (deds f duds sigma u)
       :tm (dtmds u)))))

;; Spinor functions
(defun recoverrotor3d (fs es)
  "Recover a basis given new and original bases"
  (let ((psi (+ (apply #'+ (mapcar #'*g fs (apply #'recipbvs es))) 1)))
    (if (zerogp psi)
	(error "zero psi (180 deg rotation) unsupported")
	(unitg psi))))
(defun recoverspinor3d (r fs es)
  "Recover a spinor given radius, new basis, and original basis"
  (* (recoverrotor3d fs es) (sqrt r)))
(defun rvbasis (rv vv)
  "Return a set of basis vectors derived from position and velocity"
  (let* ((mombv (*o rv vv))
	 (x (unitg rv))
	 (y (unitg (*i rv mombv)))
	 (z (when (= 3 (dimension rv) (dimension vv))
	      (*x2 x y))))
    (if (= 2 (dimension rv) (dimension vv)) ; 2D or 3D?
	(list x y)
	(list x y z))))
(defmethod rv2u ((rv g) (vv g) (basis list))
  (recoverspinor3d (norme rv) (rvbasis rv vv) basis))
(defmethod rv2dudt ((vv g) (u g) (basis1 g))
  (* (*g3 vv u basis1)
     (/ (* 2 (norme2 u)))))
(defun rv2spinors (rv vv basis)
  "Convert position and velocity vectors to spinor and spinor time derivative"
  (let ((u (rv2u rv vv basis)))
    (make-hash
     :u u
     :dudt (rv2dudt vv u (first basis)))))
(defun duds2dt (duds u)
  "Convert spinor time derivative (also given spinor) to s derivative"
  (/ duds (norme2 u)))
(defun dudt2ds (dudt u)
  "Convert spinor s derivative (also given spinor) to t derivative"
  (* dudt (norme2 u)))
(defun spinor2r (u basis1)
  "Given spinor and 1st basis vector, return corresponding position vector"
  (spin basis1 u))
(defun spinors2v (dudt u basis1)
  "Given spinor time derivative, spinor, and 1st basis vector, return corresponding velocity vector"
  (* (*g3 dudt basis1 (revg u)) 2d0))
(defun spinors2energy (u duds mu)
  "Given spinor, spinor s-derivative, and gravitational parameter, return orbit energy"
  (/ (- (* 2 (norme2 duds)) mu) (norme2 u)))

;; KSH conversion functions
(defun rv2ksh (rv vv basis mu)
  "Given position, velocity, list of basis vectors, and gravitational parameter, return KSH state initialized as time=0 and s=0"
  (with-keys (u dudt) (rv2spinors rv vv basis)
    (let* ((duds (dudt2ds dudt u))
	   (e (spinors2energy u duds mu)))
      (make-hash
       :alpha u 
       :beta (/ duds (sqrt (- (/ e 2))))
       :e e
       :tm 0))))
(defun ksh2spinors (x s)
  "Convert KSH state to spinor & spinor s-derivative"
  (with-keys (alpha beta e tm) x
    (let ((w0 (w0 e)))
      (make-hash
       :u (u alpha beta w0 s) ; spinor
       :duds (duds beta w0))))) ; spinor s-derivative
(defun ksh2rv (x s sigma0)
  "Calculate position & velocity vectors from KSH state, s, and initial orbit position unit vector (sigma0)"
  (with-keys (alpha beta e tm) x
    (with-keys (u duds) (ksh2spinors x s)
      (let ((sigma (spin sigma0 alpha)))
	(make-hash
	 :r (spinor2r u sigma) ; position
	 :v (spinors2v (duds2dt duds u) u sigma)))))) ; velocity

;; Classical orbit elements
(defun coe2rv (sma ecc inc raan aop truan mu basis)
  "Convert cassical orbital elements to position & velocity vectors.
SMA semi major axis
ECC eccentricity
INC inclination
RAAN right ascension of ascending node
AOP argument of perigee
TRUAN true anomaly
MU gravitational parameter
BASIS list of 3 orthogonal basis vectors to express position & velocity in"
  (let* ((r-raan (rotor (*o (first basis) (second basis)) raan))
	 (basis-raan (mapcar #'(lambda (x) (rot x r-raan)) basis))
	 (r-inc (rotor (*o (second basis-raan) (third basis-raan)) inc))
	 (basis-inc (mapcar #'(lambda (x) (rot x r-inc)) basis-raan))
	 (r-aop (rotor (*o (first basis-inc) (second basis-inc)) aop))
	 (basis-aop (mapcar #'(lambda (x) (rot x r-aop)) basis-inc))
	 (r-truan (rotor (*o (first basis-aop) (second basis-aop)) truan))
	 (basis-truan (mapcar #'(lambda (x) (rot x r-truan)) basis-aop))
	 (p (* sma (- 1 (* ecc ecc))))
	 (r (/ p (+ 1 (* ecc (cos truan)))))
	 (ruv (first basis-truan))
	 (rv (* ruv r))
	 (angm (sqrt (/ p mu)))
	 (angmbv (* (*o (first basis-truan) (second basis-truan)) angm))
	 (eccuv (first basis-aop))
	 (eccv (* eccuv ecc))
	 (vv (* (*g (invv angmbv) (+ eccv ruv))
		  mu)))
    (values (graden rv 1) (graden vv 1))))

;; position & velocity to classical orbit elements
(defun energy-rv (r v mu)
  "Orbit energy from position, velocity, & gravitational parameter"
  (- (/ (* v v) 2)
     (/ mu r)))
(defun sma-rv (r v mu)
  "semi-major axis from position, velocity, & gravitational parameter"
  (/ mu (* -2 (energy-rv r v mu))))
(defun eccv-rv (rv vv mu)
  "eccentricity vector from position, velocity, & gravitational parameter"
  (/ (- (* rv (- (norme2 vv) (/ mu (norme rv))))
	(* vv (scalar (*i rv vv))))
       mu))
(defun mombv-rv (rv vv)
  "orbital momentum bivector from position & velocity"
  (*o rv vv))
(defun nodev-rv (mombv basis)
  "ascending node vector from momentum and basis"
  (- (*i (third basis) mombv)))
(defun inc-rv (mombv basis)
  "inclination from momentum bivector and basis"
  (acos (/ (scalar (*i (third basis) (dual mombv)))
	   (norme mombv))))
(defun raan-rv (nodev basis)
  "right ascension of ascending node from node vector and basis"
  (let ((tmp (acos (scalar (*i (first basis) (unitg nodev))))))
    (if (< 0 (scalar (*i (second basis) nodev)))
	(- (* 2 pi) tmp)
	tmp)))
(defun aop-rv (nodev eccv basis)
  "argument of perigee from node vector, eccentricity vector, and basis"
  (let ((tmp (acos (scalar (*i (unitg nodev) (unitg eccv))))))
    (if (< 0 (scalar (*i (third basis) eccv)))
	(- (* 2 pi) tmp)
	tmp)))
(defun truan-rv (eccv rv vv)
  "true anomaly from eccentricity, position, and velocity vectors"
  (let ((tmp (acos (scalar (*i (unitg eccv) (unitg rv))))))
    (if (< 0 (scalar (*i rv vv)))
	(- (* 2 pi) tmp)
	tmp)))
(defun rv2coe (rv vv mu basis)
  "Return classical orbital elements given position & velocity vectors, gravitational parameter, and a list of 3 basis vectors"
  (let* ((mombv (mombv-rv rv vv))
	 (nodev (nodev-rv mombv basis))
	 (eccv (eccv-rv rv vv mu)))
    (make-hash
     :sma (sma-rv (norme rv) (norme vv) mu)
     :ecc (norme eccv)
     :inc (inc-rv mombv basis)
     :raan (raan-rv nodev basis)
     :aop (aop-rv nodev eccv basis)
     :truan (truan-rv eccv rv vv))))

;; Test KSH equations of motion
(defparameter *kshtest*
  (make-hash*
   basis (list (ve3 :c1 1) (ve3 :c10 1) (ve3 :c100 1))
   sigma0 (first basis)
   forcefun #'(lambda (s x) (ve3))
   r0 (ve3 :c1 1)
   v0 (ve3 :c10 1.1)
   mu 1
   x0 (rv2ksh r0 v0 basis mu)
   s0 0
   sf (* pi 2))
  "Test data for KSH equations of motion")

(defun testksh (data)
  (with-keys (s0 sf x0) data
    (let ((*kshparam* data))
      (rka #'ksheom s0 sf x0))))

(defparameter *kshsailtest*
  (make-hash*
   basis (list (ve3 :c1 1) (ve3 :c10 1) (ve3 :c100 1))
   sigma0 (first basis)
   alpha 0
   delta 0
   mu 1
   beta 0.1
   forcefun #'(lambda (s x) 
		(with-keys (r v) (ksh2rv x s sigma0)
		  (destructuring-bind (posuv tanuv orbuv) (rvbasis r v)
		    (let ((normuv (+ (* posuv (cos alpha))
				     (* orbuv (* (sin alpha) (cos delta)))
				     (* tanuv (* (sin alpha) (sin delta))))))
		      (* normuv
			 (/ (* beta mu (expt (scalar (*i posuv normuv)) 2))
			    (norme2 r)))))))
   r0 (ve3 :c1 1)
   v0 (ve3 :c10 1)
   x0 (rv2ksh r0 v0 basis mu)
   s0 0
   sf (* pi 2))
  "Test data for KSH equations of motion with solar sail")

(defparameter *ksh-forcefun-sail-normal*
  #'(lambda (s x)
      (with-keys (r v) (ksh2rv x s *ksh-sigma0*)
	(sail-normal-force r *mu* *lightness*))))

(defvar *ephemeris*
  (make-hash
   :alpha (re2 :c0 1)
   :beta (re2 :c11 -1)
   :e -0.5))

(defun planet-ksh-eom (s x)
  (with-keys (alpha beta e) *ephemeris*
    (dtmds (u alpha beta (w0 e) s))))

