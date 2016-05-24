
#	*		*		*	EKF QUATERNION *	*		*

ekf.calibration.indirect  <- function(m, initMean=NA)
{
#[m_, coefs_, tr_] = ekf.calibration.indirect(m)
#Performs the estimation of the calibration coefs by complementary EKF
#   Input:
#      m - Calibration data, recorded while rotating corresponding
#      sensor in 3D
#      initMean - Initial guess for coefs
#   Output:
#      coefs[1x12] - vector of sensor's calibration coeffs
#      m_ - calibrated data
#      tr_ - Covariance matrix trace

# Filter Parameters
covW <- diag(12) * 1e-7
covV <- 1e0
if (is.na(initMean)) initMean = matrix(0,nrow=12,ncol=1)

initCov = diag(12)*1e-2

#EKF parameters
theFltr.Q <- (covW)
theFltr.R <- (covV)
theFltr.x <- (initMean)
theFltr.P <- (initCov)

#Data logs
Nsim <- dim(m)[1]
coefs_ <- matrix(0,nrow=Nsim,ncol=12)
m_ <- matrix(0,nrow=Nsim,ncol=3)
tr_ <- matrix(0,nrow=Nsim,ncol=1)

#Initial values
dB <- matrix(0,nrow=3,ncol=1)
dB0 <- matrix(0,nrow=3,ncol=1)
dG <- matrix(0,nrow=6,ncol=1)

#Filtering
H <- matrix(0,nrow=1,ncol=12)
for (n in 1:Nsim)
    {
    dB[1]  <- dB[1]+theFltr.x[1]
    dB[2]  <- dB[2]+theFltr.x[2]
    dB[3]  <- dB[3]+theFltr.x[3]
    dG[1,1] <- dG[1,1]+theFltr.x[4]
    dG[2,1] <- dG[2,1]+theFltr.x[5]
    dG[3,1] <- dG[3,1]+theFltr.x[6]
    dG[4,1] <- dG[4,1]+theFltr.x[7]
    dG[5,1] <- dG[5,1]+theFltr.x[8]
    dG[6,1] <- dG[6,1]+theFltr.x[9]
    dB0[1,1] <- dB0[1,1]+theFltr.x[10]
    dB0[2,1] <- dB0[2,1]+theFltr.x[11]
    dB0[3,1] <- dB0[3,1]+theFltr.x[12]
        
    B <- matrix(c(dB[1], dG[1], dG[2], dG[3], dB[2], dG[4], dG[5], dG[6], dB[3]),ncol=3,nrow=3,byrow=TRUE)
    
    c <- (diag(3)-B) %*% (matrix(m[n,],nrow=3,ncol=1) - dB0)
    
    I <- diag(12)
    Z <- 1.0-sqrt(c[1]*c[1]+c[2]*c[2]+c[3]*c[3])
    
    c1 <- dB[1]
    c2 <- dB[2]
    c3 <- dB[3]
    c4 <- dG[1]
    c5 <- dG[2]
    c6 <- dG[3]    
    c7 <- dG[4]
    c8 <- dG[5]
    c9 <- dG[6]        
    c10 <- dB0[1,1]
    c11 <- dB0[2,1]
    c12 <- dB0[3,1]
    
    x <- c[1]
    y <- c[2]
    z <- c[3]

    H[1,1]  <- ((c10 - x)*(c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1)))/((c4*(c11 - y) + c5*(c12 - z) +
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2)
    H[1,2]  <- ((c11 - y)*(c6*(c10 - x) + c7*(c12 - z) +
        (c11 - y)*(c2 - 1)))/((c4*(c11 - y) + c5*(c12 - z) +
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) +
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) +
        (c12 - z)*(c3 - 1))^2)^(1/2)
    H[1,3]  <- ((c12 - z)*(c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1)))/((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2)
    H[1,4]  <- ((c11 - y)*(c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1)))/((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2)
    H[1,5]  <- ((c12 - z)*(c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1)))/((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2)
    H[1,6]  <- ((c10 - x)*(c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1)))/((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2)
    H[1,7]  <- ((c12 - z)*(c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1)))/((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2)
    H[1,8]  <- ((c10 - x)*(c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1)))/((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2) 
    H[1,9]  <- ((c11 - y)*(c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1)))/((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2)
    H[1,10] <- (2*c6*(c6*(c10 - x) + c7*(c12 - z) + (c11 - y)*(c2 - 1)) +
        2*c8*(c8*(c10 - x) + c9*(c11 - y) + (c12 - z)*(c3 - 1)) +
        2*(c1 - 1)*(c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1)))/(2*((c4*(c11 - y) + c5*(c12 - z) +
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2))
    H[1,11] <- (2*c4*(c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1)) + 2*c9*(c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1)) + 2*(c2 - 1)*(c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1)))/(2*((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2))
    H[1,12] <- (2*c5*(c4*(c11 - y) + c5*(c12 - z) + (c10 - x)*(c1 - 1)) +
        2*c7*(c6*(c10 - x) + c7*(c12 - z) + (c11 - y)*(c2 - 1)) +
        2*(c3 - 1)*(c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1)))/(2*((c4*(c11 - y) + c5*(c12 - z) + 
        (c10 - x)*(c1 - 1))^2 + (c6*(c10 - x) + c7*(c12 - z) + 
        (c11 - y)*(c2 - 1))^2 + (c8*(c10 - x) + c9*(c11 - y) + 
        (c12 - z)*(c3 - 1))^2)^(1/2))
    
    K <- mrdivide((theFltr.P %*% t(H)) , (H %*% theFltr.P %*% t(H) + theFltr.R))
    theFltr.P <- (I - K %*% H) %*% theFltr.P %*% t(I - K %*% H) + K %*% theFltr.R %*% t(K)
    theFltr.x <- K %*% Z
    
    coefs_[n,] <-  c(dB[1], dB[2], dB[3], dG[1], dG[2], dG[3], dG[4], dG[5], dG[6], dB0[1], dB0[2], dB0[3])

    tr_[n,] <- sum(diag(theFltr.P))
    m_[n,] <- c
}

list(m_=m_, coefs_=coefs_, tr_=tr_)
}

#	*		*		*	svd calibration *	*		*

svd.calibration <- function(X)
{
#X=m

#[m_, coefs_, tr_] = svd.calibration(m)
#Performs the estimation of the calibration coefs by complementary EKF
#   Input:
#      X - Calibration data, recorded while rotating corresponding
#      sensor in 3D
#   Output:
#      coefs[1x12] - vector of sensor's calibration coeffs
#      X_ - calibrated data

# using Merayo technique with a non iterative algoritm
# J.Merayo et al. "Scalar calibration of vector magnemoters"
# Meas. Sci. Technol. 11 (2000) 120-132.
#              
# The calibration tries to find the best 3D ellipsoid that fits the data set
# and returns the parameters of this ellipsoid
#
# Ellipsoid equation : (v-c)'*(U'*U)(v-c) = 1 
# with v a rough triaxes sensor  measurement
#
# calibrated measurement w = U*(v-c)
#
# author : Alain Barraud, Suzanne Lesecq 2008
#
#########################################

coefs = matrix(0,nrow=12,ncol=1)

N = dim(X)[1]
if (N<=10)     stop('not enough data for calibration !!')

#Array for the calibrated result
X_ = matrix(0,nrow=N,ncol=3)

# # write  the ellipsoid equation as D*p=0
# # the best parameter is the solution of min||D*p|| with ||p||=1
# # form D matrix from X measurements

x = X[,1]; y = X[,2]; z = X[,3] 
D = cbind(x^2, y^2, z^2, x * y, x * z, y * z, x, y, z, 1)
#matrix(1,nrow=N,ncol=1)

d2<-qr(D)
#D=upper.tri(d2[["qr"]], diag = TRUE)#avoids to compute the svd of a large matrix
D_short = d2[["qr"]][1:10,]#D[1:10,]
D_short[lower.tri(D_short, diag = FALSE)]=0
D_short = matrix(as.numeric(D_short),nrow=10,ncol=10,byrow=FALSE)
V = svd(D_short,0)$v#because usually N may be very large

p = V[,dim(V)[1]]
if (p[1]<0) p <- -p
p = matrix(p,ncol=1)
# the following matrix A(p) must be positive definite
# The optimization done by svd does not include such a constraint
# With "good" data the constraint is allways satisfied
# With too poor data A may fail to be positive definite
# In this case the calibration fails
#
A = matrix(c(p[1], p[4]/2, p[5]/2, p[4]/2, p[2], p[6]/2, p[5]/2, p[6]/2, p[3]),nrow=3,ncol=3,byrow=TRUE)

tmp = fchol.svdCalib(A)
U<-tmp$A
ok<-tmp$ok
if (!ok) stop('calibration fails too poor data!!')

b = matrix(c(p[7],p[8],p[9]),nrow=3,ncol=1)

v = Utsolve.svdCalib(U,b/2)

d = p[10]
s = 1 / sqrt(v %*% t(v)-d)

c =-t(Usolve.svdCalib(U,v))#ellipsoid center

U =  c(s) * U #shape ellipsoid parameter
coefs[1] = 1-U[1,1]
coefs[2] = 1-U[2,2]
coefs[3] = 1-U[3,3]
coefs[4] = -U[1,2]
coefs[5] = -U[1,3]
coefs[7] = -U[2,3]
coefs[6] = 0
coefs[8] = 0
coefs[9] = 0
coefs[10] = c[1]
coefs[11] = c[2]
coefs[12] = c[3]

#Apply the estimated calibration coefficients
B = diag(3)-matrix(c(coefs[1], coefs[4], coefs[5], coefs[6], coefs[2], coefs[7], coefs[8], coefs[9], coefs[3]),nrow=3,ncol=3,byrow=TRUE)
B0 = matrix(c(coefs[10], coefs[11], coefs[12]),nrow=3,ncol=1,byrow=TRUE)

for (n in 1:N) X_[n,] = t(B %*% (matrix(c(X[n,]),nrow=3,ncol=1,byrow=TRUE)-B0))

list(X_=X_, coefs=coefs)
}

fchol.svdCalib<-function(A)
{
# performs Cholesky factoristation
A[1,1:3] = A[1,1:3] / sqrt(A[1,1])
A[2:3,1] = 0
A[2,2:3] = A[2,2:3] - t(A[1,2]) %*% A[1,2:3]
if (A[2,2]<=0) return (list(A=NA,ok=FALSE))#A is not positive definite
A[2,2:3] = A[2,2:3] / sqrt(A[2,2])
A[3,2] = 0
A[3,3:3] = A[3,3:3] - t(A[1:2,3]) %*% A[1:2,3:3]
if (A[3,3]<=0) return (list(A=NA,ok=FALSE))#A is not positive definite
A[3,3:3] = A[3,3:3] / sqrt(A[3,3])
#A[3+1:3,3] = 0
ok=TRUE
list(A=A,ok=ok)
}

Utsolve.svdCalib<-function(U,b)
{
# solves U'*x=b
x = matrix(0,nrow=1,ncol=3)
x[1] = b[1]/U[1,1]
x[2] = (b[2]-x[1] %*% U[1,2]) / U[2,2]
x[3] = (b[3]-x[1:2] %*% U[1:2,3]) / U[3,3]
x
}

Usolve.svdCalib<-function(U,b)
{
# solves U*x=b
x = matrix(0,nrow=1,ncol=3)
x[3] = b[3]/U[3,3]
x[2] = (b[2]-U[2,3] %*% t(x[3]))/U[2,2]
x[1] = (b[1]-U[1,2:3] %*% matrix(x[2:3],ncol=1))/U[1,1]
x
}

###########################################

#	*		*		*	EKF QUATERNION *	*		*

ahrs.EKF.QUATERNION <- function(Filter, Sensors, Parameters)
{
#Filter <- ahrs_EKF_QUATERNION(Filter, Sensors, Parameters)
# Function implements the EKF-based AHRS algorithm based on measurements
# from three-component accelerometer with orthogonal axes, vector
# magnetometer and three-axis gyroscope.
# Estimates the current quaternion attitude
#
#   Input arguments:
#   Filter - data structure for Extended Kalman Filter
#   Filter.x   State vector [3x1]
#   Filter.P   Covariance matrix [3x3]
#   Filter.Q   System noise matrix [3x3]
#   Filter.R   Measurement noise matrix [6x6]
#   Q - initial attitude quaternion value [1x4]
#   Parameters -  AHRS Parameters
#   Parameters.mn      Magnetic Field Vector In Navigation Frame [3x1], |m|
#   <- 1
#   Parameters.an      Acceleration vector In Navigation Frame [3x1], g
#   Parameters.dt      Sampling period, 1/Hz
#   Sensors - sensors data structure
#   Sensors.w    current calibrated gyroscope measurement [3x1], rad/sec
#   Sensors.a    current calibrated accelerometer measurement [3x1], g
#   Sensors.m    current calibrated magnetometer measurement [3x1], |m| <- 1
#
#   Output arguments:
#   Q - estimated attitude quaternion [1x4]
#   Filter - data structure for an Extended Kalman Filter
#   dw - estimated gyroscopes bias [1x3]

ab <- Sensors$a/norm(Sensors$a,'f')
mb <- Sensors$m/norm(Sensors$m,'f')
wb <- Sensors$w

an <- Parameters$an
mn <- Parameters$mn
dT <- Parameters$dT

A <- DfDx(Filter$x,wb,dT)

Filter$P <- A %*% Filter$P %*% t(A) + Filter$Q

Filter$x <- systemEKF(Filter$x, wb,dT)

dH <- DhDx(Filter$x, an, mn)
Yhat <- measureEKF(Filter$x,an,mn)

Y <- matrix(c(ab, mb),ncol=1)
v <- Y-t(Yhat)
    
nFilter <- KF.cholesky.update(Filter$x, Filter$P,v, Filter$R, dH)
Filter$P <- nFilter$P
Filter$x <- nFilter$x

return( Filter )
}

systemEKF <- function(x,wb,dT)
{
q <- c(x[1], x[2], x[3], x[4])
dw <- c(x[5],x[6],x[7])

#Correct angular rate
What <- wb-dw

What <-matrix(What,ncol=3,byrow=FALSE)

#calculate quaternion
q  <- Qrot(q, What, dT)

q <- Qnormalize(q)

x[1:4] <- q
x
}

measureEKF <- function(x,an,mn)
{
#Measurements
q <- c(x[1], x[2], x[3], x[4])
ab_hat <- vectQrot( q, t(an))
mb_hat <- vectQrot( q, t(mn))
y <- t(c(ab_hat, mb_hat))
y
}

DfDx <- function(x,wb,dT)
{

#x=Filter$x

q1 <- x[1]
q2 <- x[2]
q3 <- x[3]
q4 <- x[4]

dw1 <- x[5]
dw2 <- x[6]
dw3 <- x[7]

w1 <- wb[1]
w2 <- wb[2]
w3 <- wb[3]

#System Jacobian
matrix(c(1,  dT*(dw1/2 - w1/2),  dT*(dw2/2 - w2/2),  dT*(dw3/2 - w3/2),  (dT*q2)/2,  (dT*q3)/2,  (dT*q4)/2, 
    -dT*(dw1/2 - w1/2), 1, -dT*(dw3/2 - w3/2),  dT*(dw2/2 - w2/2), -(dT*q1)/2,  (dT*q4)/2, -(dT*q3)/2, 
    -dT*(dw2/2 - w2/2),  dT*(dw3/2 - w3/2), 1, -dT*(dw1/2 - w1/2), -(dT*q4)/2, -(dT*q1)/2,  (dT*q2)/2, 
    -dT*(dw3/2 - w3/2), -dT*(dw2/2 - w2/2),  dT*(dw1/2 - w1/2), 1,  (dT*q3)/2, -(dT*q2)/2, -(dT*q1)/2, 
    0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 0, 1),nrow=length(x),byrow=TRUE)
}

DhDx <- function(x, an, mn)
{
mx = mn[1]
my = mn[2]
mz = mn[3]

ax = an[1]
ay = an[2]
az = an[3]

q1 = x[1]
q2 = x[2]
q3 = x[3]
q4 = x[4]
#Measurement Jacobian

matrix(c( 2*q1*ax + 2*q4*ay - 2*q3*az, 2*q2*ax + 2*q3*ay + 2*q4*az,
    2*q2*ay - 2*q3*ax - 2*q1*az, 2*q1*ay - 2*q4*ax + 2*q2*az, 0, 0, 0 ,
    2*q1*ay - 2*q4*ax + 2*q2*az, 2*q3*ax - 2*q2*ay + 2*q1*az,
    2*q2*ax + 2*q3*ay + 2*q4*az, 2*q3*az - 2*q4*ay - 2*q1*ax, 0, 0, 0 ,
    2*q3*ax - 2*q2*ay + 2*q1*az, 2*q4*ax - 2*q1*ay - 2*q2*az,
    2*q1*ax + 2*q4*ay - 2*q3*az, 2*q2*ax + 2*q3*ay + 2*q4*az, 0, 0, 0  ,   
    2*q1*mx + 2*q4*my - 2*q3*mz, 2*q2*mx + 2*q3*my + 2*q4*mz,
    2*q2*my - 2*q3*mx - 2*q1*mz, 2*q1*my - 2*q4*mx + 2*q2*mz, 0, 0, 0 ,
    2*q1*my - 2*q4*mx + 2*q2*mz, 2*q3*mx - 2*q2*my + 2*q1*mz,
    2*q2*mx + 2*q3*my + 2*q4*mz, 2*q3*mz - 2*q4*my - 2*q1*mx, 0, 0, 0 ,
    2*q3*mx - 2*q2*my + 2*q1*mz, 2*q4*mx - 2*q1*my - 2*q2*mz,
    2*q1*mx + 2*q4*my - 2*q3*mz, 2*q2*mx + 2*q3*my + 2*q4*mz, 0, 0, 0 ) ,ncol=length(x),byrow=TRUE)
}

KF.cholesky.update <- function(x,P,v,R,H)
# Calculate the KF (or EKF) update given the prior state [x,P]
# the innovation [v,R] and the (linearised) observation model H.
# The result is calculated using Cholesky factorisation, which
# is more numerically stable than a naive implementation.
# Tim Bailey 2003
# Developed by Jose Guivant
{

#x=Filter$x; P=Filter$P; R=Filter$R;H= dH
PHt <- P %*% t(H)
S <- H %*% PHt + R

S <- (S+t(S)) * 0.5 # make symmetric
SChol <- chol(S)

W1 <-  mrdivide(PHt, SChol )   #PHt / SChol    right division
W <- W1 %*% solve(t(SChol))

x <- x + W %*% (v) # update
#write(format(c(Filter$x), scientific = TRUE,digits=8), file='EKF quaternion FilterXsystemEKF_R.txt',ncolumns =7,append=TRUE)
P <- P - W1 %*% t(W1)
list(x=x,P=P)
}


#	*		*		*	LKF EULER *	*		*

ahrs.LKF.EULER <- function(Sensors, State, Parameters)
{

## Get previous state
q <- State$q
if (any(is.na(q))) stop('Error! Q is NA.')
dB <- State$dB
dG <- State$dG
dw <- State$dw
P <- State$P

Wb <- Sensors$w
Ab <- Sensors$a
Mb <- Sensors$m

dT <- Parameters$dT

#Correct Gyroscopes for estimated biases and scale factors
B <- matrix(c(dB[1], dG[1], dG[2], dG[3], dB[2], dG[4], dG[5], dG[6], dB[3]),nrow=3,ncol=3,byrow=TRUE)

#cat('B=',B, ' wb=', c(Wb), ' dw=', dw,'\n')
OmegaBib <- (diag(3)-B) %*% Wb-dw

## Calculate quaternion
q  <- Qrot(q, matrix(OmegaBib,1), dT)

#Q=q;w=matrix(OmegaBib,1)


#Calculate DCM
Cnb <- Q2DCM(Qnormalize(q))
Cbn <- t(Cnb)

#Calculate Gyro Angles
tmp <- DCM2EA(Cnb)
Psi <- tmp[1]
Theta <- tmp[2]
Gamma <- tmp[3]

#Accelerometer angles
ThetaAcc <-  atan2(Ab[1],sqrt(Ab[2]^2+Ab[3]^2))
GammaAcc <- -atan2(Ab[2],sqrt(Ab[1]^2+Ab[3]^2))

#Horizontal projection of magnetic field 
CbnH <- EA2DCM(c(0,ThetaAcc,GammaAcc))
Mh <- t(CbnH) %*% Mb
# Magnetic Heading
PsiMgn <- -atan2(Mh[2],Mh[1])+Parameters$declination


#System matrix
A <- matrix(0,ncol=15,nrow=15)
I <- diag(15)
G <- matrix(c(OmegaBib[2], OmegaBib[3], 0, 0, 0, 0, 0, 0, OmegaBib[1], OmegaBib[3], 0, 0, 0, 0, 0, 0, OmegaBib[1], OmegaBib[2]),nrow=3,ncol=6,byrow=TRUE)
 
A[1:3,4:6]   <- -Cbn
A[1:3,7:9]   <- -Cbn %*% diag(c(OmegaBib))
A[1:3,10:15] <- -Cbn %*% G
A[4:6,4:6] <- -diag(c(1/350, 1/350, 1/350))
A[7:9,7:9] <- -diag(c(1/500, 1/500, 1/500))
A[10:15,10:15] <- -diag(c(1/500, 1/500, 1/500, 1/500, 1/500, 1/500))
F <- I+A * dT

#Measurement matrix
dPsi <- PsiMgn-Psi
if (dPsi > pi) dPsi <- dPsi - (2 * pi)
if (dPsi < -pi) dPsi <- dPsi + (2 * pi)

dTheta <- ThetaAcc-Theta
if (dTheta > pi) dTheta <- dTheta-(2*pi)
if (dTheta < -pi) dTheta <- dTheta+(2*pi)

dGamma <- GammaAcc-Gamma
if (dGamma > pi) dGamma <- dGamma-(2*pi)
if (dGamma < -pi) dGamma <- dGamma+(2*pi)

z <- matrix(c(dGamma, dTheta, dPsi),nrow=3,ncol=1)

H <- matrix(0,nrow=3,ncol=15)
if (abs(sqrt(Ab[1]^2+Ab[2]^2+Ab[3]^2) - 1) < 0.1) H[1:2,1:2] <- diag(2)

if (abs(sqrt(Mb[1]^2+Mb[2]^2+Mb[3]^2) - 1) < 0.5) H[3,3] <- 1

#Kalman Filter
Q <- Parameters$Q
R <- Parameters$R

#P is the problem!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

P <- F %*% P %*% t(F) + Q
# K <- (P %*% t(H))/(H %*% P %*% t(H) +R)
# P <- (I-K %*% H) %*% P %*% t(I-K %*% H) + K %*% R %*% t(K)
# xf <- K %*% z
tmp <- KF.cholesky.update(0,P,z,R,H)
xf <- tmp$x
P <- tmp$P
df_hat <- xf[1:3]
dw_hat <- xf[4:6]
dB_hat <- xf[7:9]
dG_hat <- xf[10:15]

#Update gyroscope biases
dw <- dw+dw_hat

#Update gyroscope scale factors
dB <- dB+dB_hat

#Update gyroscope axes misalignments
dG <- dG+dG_hat

#DCM error estimated
dCbn <- matrix(c ( 0, df_hat[3], -df_hat[2], -df_hat[3], 0, df_hat[1], df_hat[2], -df_hat[1], 0 ), nrow=3,ncol=3,byrow=TRUE)
Cbn <- Cbn %*% (diag(3)-dCbn)
q <- DCM2Q(t(Cbn)) #rotMat2quatern(t(Cbn))

#DCM=t(Cbn);tol = 10 * .Machine$double.eps; ichk=FALSE; ignoreAllChk=FALSE




tmp <- DCM2EA(t(Cbn))
Psi <- tmp[1]
Theta <- tmp[2]
Gamma <- tmp[3]
Attitude <- matrix(c (Psi, Theta, Gamma, PsiMgn, ThetaAcc, GammaAcc), nrow=6,ncol=1)

## Set Current State
State$q <- q
State$dG <- dG
State$dB <- dB
State$dw <- dw
State$P <- P

list(Attitude=Attitude, State=State)
}

#	*		*		*	LKF QUATERNION *	*		*

ahrs.LKF.QUATERNION <- function(Filter, Sensors, q, Parameters, dw)
{#dw<-dw_hat
I <- diag(6)
A <-matrix(0,ncol=6,nrow=6)
z <-matrix(0,ncol=1,nrow=6)
H <-matrix(0,ncol=6,nrow=6)

## Correct angular rate
W_hat <- Sensors$w - t(dw)

## Calculate quaternion
q <- Qrot(q, W_hat, Parameters$dt)

## System matrix
A[1:3,1:3] <- matrix(c(0, W_hat[3], -W_hat[2], -W_hat[3], 0, W_hat[1], W_hat[2], -W_hat[1], 0),nrow=3,ncol=3,byrow=TRUE)
A[1:3,4:6] <- -diag(3)
F <- I + A * Parameters$dt

## Estimated vector values
ab_hat = vectQrot( q, Parameters$an)
mb_hat = vectQrot( q, Parameters$mn)

## Measurement vector
dA = (Sensors$a)/norm(Sensors$a,'f') - ab_hat
dM = (Sensors$m)/norm(Sensors$m,'f') - mb_hat
z[1:3,1] = dA
z[4:6,1] = dM

## Measurement matrix with gating
threshold_a =  abs(sqrt(Sensors$a[1]^2+Sensors$a[2]^2+Sensors$a[3]^2)-1.0)
if (threshold_a < 0.1)
{
H[1,2] = -ab_hat[3]
H[1,3] =  ab_hat[2]
H[2,1] =  ab_hat[3]
H[2,3] = -ab_hat[1]
H[3,1] = -ab_hat[2]
H[3,2] =  ab_hat[1]
}

threshold_m = abs(sqrt(Sensors$m[1]^2+Sensors$m[2]^2+Sensors$m[3]^2)-1.0)
if (threshold_m < 0.5)
{
H[4,2] = -mb_hat[3]
H[4,3] =  mb_hat[2]
H[5,1] =  mb_hat[3]
H[5,3] = -mb_hat[1]
H[6,1] = -mb_hat[2]
H[6,2] =  mb_hat[1]
}

## Kalman Filter
Filter$P = F %*% Filter$P %*% t(F) + Filter$Q
K = mrdivide((Filter$P %*% t(H)), (H %*% Filter$P %*% t(H) + Filter$R))
Filter$P = (I-K %*% H) %*% Filter$P %*% t(I-K %*% H) + K %*% Filter$R %*% t(K)
x = K %*% z

#Update Gyroscopes Bias estimation
dw = dw + x[4:6]

## Correct quaternion
qe  = matrix(c(x[1]/2, x[2]/2, x[3]/2, 0),nrow=1, ncol=4) # ncol=3
qqe = matrix(c(sqrt(1-(qe[1]^2+qe[2]^2+qe[3]^2)),qe[1], qe[2], qe[3]),nrow=1, ncol=4)
q = q %Q*% qqe

list( Filter=Filter, q=q, dw=dw )
}

#	*		*		*	LKF VMATCH *	*		*

ahrs.LKF.VMATCH <- function(Filter, Sensors, q, Parameters)
{
I <- diag(3)
F <-diag(c(Parameters$lambda, Parameters$lambda, Parameters$lambda))

# Estimated vector values
mn_hat = vectQrot(q,(Sensors$m)/norm(Sensors$m,'f'))
an_hat = vectQrot(q,(Sensors$a)/norm(Sensors$a,'f'))

## Measurement vector
z = matrix(0,ncol=1,nrow=6)
dm = Parameters$mn - mn_hat
da = Parameters$an - an_hat
z[1:3,1] = da
z[4:6,1] = dm

## Measurement matrix
H = matrix(0,ncol=3,nrow=6)
H[1,2] = -an_hat[3]
H[1,3] =  an_hat[2]
H[2,1] =  an_hat[3]
H[2,3] = -an_hat[1]
H[3,1] = -an_hat[2]
H[3,2] =  an_hat[1]
H[4,2] = -mn_hat[3]
H[4,3] =  mn_hat[2]
H[5,1] =  mn_hat[3]
H[5,3] = -mn_hat[1]
H[6,1] = -mn_hat[2]
H[6,2] =  mn_hat[1]

## Kalman Filter
Filter$P = F %*% Filter$P %*% t(F) + Filter$Q
K = mrdivide((Filter$P %*% t(H)), (H %*% Filter$P %*% t(H) + Filter$R))
Filter$P = (I-K %*% H) %*% Filter$P %*% t(I-K %*% H) + K %*% Filter$R %*% t(K)
x = K %*% z

## Correct quaternion
qe  = matrix(c(x[1]/2, x[2]/2, x[3]/2, 0), ncol=4,nrow=1) # ncol=3
qqe = matrix(c(sqrt(1-(qe[1]^2+qe[2]^2+qe[3]^2)), qe[1], qe[2], qe[3]), ncol=4,nrow=1) #ncol=1,nrow=4 
q   = q %Q*% qqe

list( Filter=Filter, q=q)
}

#	*		*		*	UKF QUATERNION *	*		*

ahrs.UKF.QUATERNION <- function(Filter, Sensors, Parameters)
{
#Filter = ahrs_UKF_QUATERNION(Filter, Sensors, Parameters)
# Function implements the UKF-based AHRS algorithm based on measurements
# from three-component accelerometer with orthogonal axes, vector
# magnetometer and three-axis gyroscope.
# Estimates the current quaternion attitude
#
#   Input arguments:
#   Filter - data structure for Unscented Kalman Filter
#   Filter.x   State vector [3x1]
#   Filter.P   Covariance matrix [3x3]
#   Filter.Q   System noise matrix [3x3]
#   Filter.R   Measurement noise matrix [6x6]

#   Parameters -  AHRS Parameters
#   Parameters.mn      Magnetic Field Vector In Navigation Frame [3x1], |m|
#   = 1
#   Parameters.an      Acceleration vector In Navigation Frame [3x1], g
#   Parameters.dt      Sampling period, 1/Hz
#   Sensors - sensors data structure
#   Sensors.w    current calibrated gyroscope measurement [3x1], rad/sec
#   Sensors.a    current calibrated accelerometer measurement [3x1], g
#   Sensors.m    current calibrated magnetometer measurement [3x1], |m| = 1
#
#   Output arguments:
#   Filter - data structure for Unscented Kalman Filter

Ab <- Sensors$a/norm(Sensors$a,'f')
Mb <- Sensors$m/norm(Sensors$m,'f')
Wb <- Sensors$w

An <- Parameters$an
Mn <- Parameters$mn
dT <- Parameters$dT

#Measurements vector
z <- matrix(c(Ab, Mb),ncol=1)

#Unscented Kalman Filter
Filter <- updateUKF(Filter,z,Wb,An,Mn,dT)
return( Filter )
}

updateUKF <- function(Filter,z,Wb,An,Mn,dT)
# UKF   Unscented Kalman Filter for nonlinear dynamic systems

# Reference: Julier, SJ. and Uhlmann, J.K., Unscented Filtering and
# Nonlinear Estimation, Proceedings of the IEEE, Vol. 92, No. 3,
# pp.401-422, 2004.
#
# By Yi Cao at Cranfield University, 04/01/2008
{
P <- Filter$P
x <- Filter$x
Q <- Filter$Q
R <- Filter$R

#number of states
L <- length(x) 
#number of measurements
m <- length(z) 
#default, tunable
ki <- 0   
#default, tunable
alpha  <-  1e0  
#default, tunable
beta <- 2                                     
#scaling factor
lambda <- alpha^2*(L+ki)-L    
#scaling factor
c <- L+lambda    
#weights for means 
Wm <- matrix(c(lambda/c, rep(0.5/c, times=2*L) ),nrow=1)

Wc <- Wm
#weights for covariance
Wc[1] <- Wc[1]+(1-alpha^2+beta)               
c <- sqrt(c)
#sigma points around x
X <- sigmas(x,P,c)

#Predict
#unscented transformation of process
u1 <- ut.system(systemUKF,X,Wm,Wc,L,Q,Wb,dT)
x1 <- u1$x
X1 <- u1$Y
P1 <- u1$P
X2 <- u1$Y1
#Update
#unscented transformation of measurments
u2 <- ut.meas(measureUKF, X1, Wm, Wc, m, R, An, Mn)
z1 <- u2$y
P2 <- u2$P
Z2 <- u2$Y1

#transformed cross-covariance
P12 <- X2 %*% diag(c(Wc)) %*% t(Z2)
K <- mrdivide(P12, P2)
#state update
x <- x1+K %*% (z-z1)
#covariance update
P <- P1-K %*% t(P12)
#Assign new filter state
Filter$P <- P
Filter$x <- x
Filter
}

sigmas <- function(x,P,c)
#Sigma points around reference point
#Inputs:
#       x: reference point
#       P: covariance
#       c: coefficient
#Output:
#       X: Sigma points
{
A <- c*t(chol(P))
Y <-matrix( rep(x ,dim(x)[1]),ncol=dim(x)[1])
X <- cbind(x, Y+A, Y-A)
#matrix(c(x, Y+A, Y-A),nrow=1)
X
}

ut.system <- function(f,X,Wm,Wc,n,Q,Wb,dT)
#Unscented Transformation
#Input:
#        f: nonlinear map
#        X: sigma points
#       Wm: weights for mean
#       Wc: weights for covraiance
#        n: number of outputs of f
#        R: additive covariance
#Output:
#        y: transformed mean
#        Y: transformed sampling points
#        P: transformed covariance
#       Y1: transformed deviations
{
L <- dim(X)[2]
x <- matrix(0,nrow=n,ncol=1)
Y <- matrix(0,nrow=n,ncol=L)

for (k in 1:L)
{
    Y[,k] <- systemUKF(X[,k],Wb,dT)#f
    x <- x+Wm[k]*Y[,k]
}
Y1 <- Y-matrix( rep(x ,dim(x)[1]*L),nrow=dim(x)[1],ncol=L)
P <- Y1 %*% diag(c(Wc)) %*% t(Y1)+Q
list(x=x,Y=Y,P=P,Y1=Y1)#[x,Y,P,Y1]
}

ut.meas <- function(f,X,Wm,Wc,n,R,An,Mn)
#Unscented Transformation
#Input:
#        f: nonlinear map
#        X: sigma points
#       Wm: weights for mean
#       Wc: weights for covraiance
#        n: number of outputs of f
#        R: additive covariance
#Output:
#        y: transformed mean
#        Y: transformed sampling points
#        P: transformed covariance
#       Y1: transformed deviations
{
L <- dim(X)[2]
y <- matrix(0,nrow=n,ncol=1)
Y <- matrix(0,nrow=n,ncol=L)
for (k in 1:L)
{
    Y[,k] <- measureUKF(X[,k],An,Mn)#f
    y <- y+Wm[k]*Y[,k]
}
Y1 <- Y- matrix( rep(y ,dim(y)[1]*L),nrow=dim(y)[1],ncol=L)#y(:,ones(1,L))
P <- Y1 %*% diag(c(Wc)) %*% t(Y1)+R
list(y=y,Y=Y,P=P,Y1=Y1)#[y,Y,P,Y1]
}

systemUKF <- function(x,wb,dT)
{
q <- matrix(x[1:4],nrow=1)
dw <- matrix(x[5:7],ncol=1)
#Correct angular rate
W_hat <- wb-dw
#calculate quaternion
W_hat <- matrix(W_hat,ncol=3,byrow=FALSE)
q <- Qrot(q, W_hat, dT)
q <- Qnormalize(q)
x[1:4] <- q
x
}

measureUKF <- function(x,an,mn)
#Measurements
{
q <- matrix(x[1:4],nrow=1)
ab_hat <- vectQrot( q, t(an))
mb_hat <- vectQrot( q, t(mn))
y <- matrix(c(ab_hat, mb_hat),ncol=1)
y
}


MadgwickAHRS<-function(MSamplePeriod, MBeta, q, Gyroscope, Accelerometer, Magnetometer)
{
#error management
if (any(is.na(c(MSamplePeriod, MBeta, q, Gyroscope, Accelerometer, Magnetometer)))) stop('Error! Input must be numeric.')
if (!is.null(dim(MSamplePeriod))) stop('Error! SamplePeriod must be a scalar.')
if (!is.null(dim(MBeta))) stop('Error! Beta must be a scalar.')
# if (!is.null(dim(q))) stop('Error! q must be a vector.')
# if (!is.null(dim(Gyroscope))) stop('Error! Gyroscope must be a vector.')
# if (!is.null(dim(Accelerometer))) stop('Error! Accelerometer must be a vector.')
# if (!is.null(dim(Magnetometer))) stop('Error! Magnetometer must be a vector.')
if (length(MSamplePeriod) !=1) stop('Error! SamplePeriod must be a scalar.')
if (length(MBeta) !=1) stop('Error! Beta must be a scalar.')
if (length(q) !=4) stop('Error! q must be a vector with 4 elements.')
if (length(Gyroscope) !=3) stop('Error! Gyroscope must be a vector with 3 elements.')
if (length(Accelerometer) !=3) stop('Error! Accelerometer must be a vector with 3 elements.')
if (length(Magnetometer) !=3) stop('Error! Magnetometer must be a vector with 3 elements.')
if (!is.numeric(MSamplePeriod))stop('Error! SamplePeriod must be numeric.')
if (!is.numeric(MBeta))stop('Error! Beta must be numeric.')
if (!is.numeric(q))stop('Error! q must be numeric.')
if (!is.numeric(Gyroscope))stop('Error! Gyroscope must be numeric.')
if (!is.numeric(Accelerometer))stop('Error! Accelerometer must be numeric.')
if (!is.numeric(Magnetometer))stop('Error! Magnetometer must be numeric.')

# Normalise accelerometer measurement
if(norm(Accelerometer,'f') == 0) stop('Error! Accelerometer norm is zero.')
Accelerometer <- (Accelerometer/norm(Accelerometer,'f'))# normalise magnitude
# Normalise magnetometer measurement
if(norm(Magnetometer,'f') == 0) stop('Error! Magnetometer norm is zero.')
Magnetometer = (Magnetometer/norm(Magnetometer,'f'))# normalise magnitude
# Reference direction of Earth's magnetic field
h <-  q %Q*% (c(0, Magnetometer) %Q*% Qconj(q))
b <- c(0,norm(cbind(h[2], h[3]),'f'), 0, h[4])
# Gradient decent algorithm corrective step
if (!is.matrix(q)) q <- matrix(q,ncol=4,byrow=FALSE)
F2<-c(2*(q[2]*q[4] - q[1]*q[3 ]) - Accelerometer[1],
2*(q[1]*q[2] + q[3]*q[4 ]) - Accelerometer[2],
2*(0.5 - q[2]^2 - q[3]^2) - Accelerometer[3],
2*b[2]*(0.5 - q[3]^2 - q[4]^2) + 2*b[4]*(q[2]*q[4] - q[1]*q[3 ]) - Magnetometer[1],
2*b[2]*(q[2]*q[3] - q[1]*q[4 ]) + 2*b[4]*(q[1]*q[2] + q[3]*q[4 ]) - Magnetometer[2],
2*b[2]*(q[1]*q[3] + q[2]*q[4 ]) + 2*b[4]*(0.5 - q[2]^2 - q[3]^2) - Magnetometer[3])
J2<-c(-2*q[3], 2*q[4], -2*q[1], 2*q[2],
2*q[2], 2*q[1], 2*q[4], 2*q[3],
0, -4*q[2], -4*q[3], 0,
-2*b[4]*q[3], 2*b[4]*q[4], -4*b[2]*q[3]-2*b[4]*q[1], -4*b[2]*q[4]+2*b[4]*q[2],
-2*b[2]*q[4]+2*b[4]*q[2], 2*b[2]*q[3]+2*b[4]*q[1], 2*b[2]*q[2]+2*b[4]*q[4], -2*b[2]*q[1]+2*b[4]*q[3],
2*b[2]*q[3], 2*b[2]*q[4]-4*b[4]*q[2], 2*b[2]*q[1]-4*b[4]*q[3], 2*b[2]*q[2])
F<-matrix(F2,ncol=1,byrow=TRUE)
J<-matrix(J2,ncol=4,byrow=TRUE)
Mstep <- t(J) %*% F
Mstep <- Mstep/norm(Mstep,'f')
# Compute rate of change of quaternion
qDot <- 0.5 * (q %Q*% cbind(0, Gyroscope[1], Gyroscope[2], Gyroscope[3])) - MBeta * matrix(Mstep,nrow=1)
# Integrate to yield quaternion
q <- q + qDot * MSamplePeriod
MQuaternion <- (q)/norm(q,'f')# normalise quaternion
return(c(MQuaternion)) #turn it into a vector
}

MadgwickIMU<-function(MSamplePeriod, MBeta, q, Gyroscope, Accelerometer)
{
#error management
if (any(is.na(c(MSamplePeriod, MBeta, q, Gyroscope, Accelerometer)))) stop('Error! Input must be numeric.')
if (!is.null(dim(MSamplePeriod))) stop('Error! SamplePeriod must be a scalar.')
if (!is.null(dim(MBeta))) stop('Error! Beta must be a scalar.')
if (!is.null(dim(q))) stop('Error! q must be a vector.')
if (!is.null(dim(Gyroscope))) stop('Error! Gyroscope must be a vector.')
if (!is.null(dim(Accelerometer))) stop('Error! Accelerometer must be a vector.')
if (length(MSamplePeriod) !=1) stop('Error! SamplePeriod must be a scalar.')
if (length(MBeta) !=1) stop('Error! Beta must be a scalar.')
if (length(q) !=4) stop('Error! q must be a vector with 4 elements.')
if (length(Gyroscope) !=3) stop('Error! Gyroscope must be a vector with 3 elements.')
if (length(Accelerometer) !=3) stop('Error! Accelerometer must be a vector with 3 elements.')
if (!is.numeric(MSamplePeriod))stop('Error! SamplePeriod must be numeric.')
if (!is.numeric(MBeta))stop('Error! Beta must be numeric.')
if (!is.numeric(q))stop('Error! q must be numeric.')
if (!is.numeric(Gyroscope))stop('Error! Gyroscope must be numeric.')
if (!is.numeric(Accelerometer))stop('Error! Accelerometer must be numeric.')

# Normalise accelerometer measurement
if(Qnorm(Accelerometer) == 0) stop('Error! Accelerometer norm is zero.')
Accelerometer <- Qnormalize(Accelerometer)# normalise magnitude

# Gradient decent algorithm corrective step
F <- c (2*(q[2]*q[4] - q[1]*q[3]) - Accelerometer[1],
2*(q[1]*q[2] + q[3]*q[4]) - Accelerometer[2],
2*(0.5 - q[2]^2 - q[3]^2) - Accelerometer[3])

J <- c(-2*q[3],	2*q[4],    -2*q[1],	2*q[2],
2*q[2],     2*q[1],     2*q[4],	2*q[3],
0,         -4*q[2],    -4*q[3],	0    )

J<-matrix(J,ncol=4,byrow=TRUE)
Mstep <- t(J) %*% F
Mstep <- Qnormalize(Mstep)
# Compute rate of change of quaternion
qDot <- 0.5 * (q %Q*% cbind(0, Gyroscope[1], Gyroscope[2], Gyroscope[3])) - MBeta * t(Mstep)
# Integrate to yield quaternion
q <- q + qDot * MSamplePeriod
MQuaternion <- Qnormalize(q)# normalise quaternion
return(c(MQuaternion)) #turn it into a vector
}

MahonykAHRS<-function(MSamplePeriod, Kp=2.0, Ki=0.005, q, Gyroscope, Accelerometer, Magnetometer)
{
#error management
if (any(is.na(c(MSamplePeriod, Kp, Ki, q, Gyroscope, Accelerometer, Magnetometer)))) stop('Error! Input must be numeric.')
if (!is.null(dim(MSamplePeriod))) stop('Error! SamplePeriod must be a scalar.')
if (!is.null(dim(Kp))) stop('Error! Kp must be a scalar.')
if (!is.null(dim(Ki))) stop('Error! Ki must be a scalar.')
if (!is.null(dim(q))) stop('Error! q must be a vector.')
if (!is.null(dim(Gyroscope))) stop('Error! Gyroscope must be a vector.')
if (!is.null(dim(Accelerometer))) stop('Error! Accelerometer must be a vector.')
if (!is.null(dim(Magnetometer))) stop('Error! Magnetometer must be a vector.')
if (length(MSamplePeriod) !=1) stop('Error! SamplePeriod must be a scalar.')
if (length(Kp) !=1) stop('Error! Kp must be a scalar.')
if (length(Ki) !=1) stop('Error! Ki must be a scalar.')
if (length(q) !=4) stop('Error! q must be a vector with 4 elements.')
if (length(Gyroscope) !=3) stop('Error! Gyroscope must be a vector with 3 elements.')
if (length(Accelerometer) !=3) stop('Error! Accelerometer must be a vector with 3 elements.')
if (length(Magnetometer) !=3) stop('Error! Magnetometer must be a vector with 3 elements.')
if (!is.numeric(MSamplePeriod))stop('Error! SamplePeriod must be numeric.')
if (!is.numeric(Kp))stop('Error! Kp must be numeric.')
if (!is.numeric(Ki))stop('Error! Ki must be numeric.')
if (!is.numeric(q))stop('Error! q must be numeric.')
if (!is.numeric(Gyroscope))stop('Error! Gyroscope must be numeric.')
if (!is.numeric(Accelerometer))stop('Error! Accelerometer must be numeric.')
if (!is.numeric(Magnetometer))stop('Error! Magnetometer must be numeric.')
eInt <- c(0, 0, 0) # integral error
# Normalise accelerometer measurement
if(Qnorm(Accelerometer) == 0) stop('Error! Accelerometer norm is zero.')
Accelerometer <- Qnormalize(Accelerometer)# normalise magnitude
# Normalise magnetometer measurement
if(Qnorm(Magnetometer) == 0) stop('Error! Magnetometer norm is zero.')
Magnetometer = Qnormalize(Magnetometer)# normalise magnitude

# Reference direction of Earth's magnetic field
h <-  q %Q*% (c(0, Magnetometer) %Q*% Qconj(q))
b <- cbind(0,Qnorm(c(h[2], h[3])), 0, h[4])
# Estimated direction of gravity and magnetic field
v <- c(2*(q[2]*q[4] - q[1]*q[3]), 
 2*(q[1]*q[2] + q[3]*q[4]),
 q[1]^2 - q[2]^2 - q[3]^2 + q[4]^2)
w <- c(2*b[2]*(0.5 - q[3]^2 - q[4]^2) + 2*b[4]*(q[2]*q[4] - q[1]*q[3]),
 2*b[2]*(q[2]*q[3] - q[1]*q[4]) + 2*b[4]*(q[1]*q[2] + q[3]*q[4]),
 2*b[2]*(q[1]*q[3] + q[2]*q[4]) + 2*b[4]*(0.5 - q[2]^2 - q[3]^2)) 

# Error is sum of cross product between estimated direction and measured direction of fields
Merror <- crossprod(Accelerometer, v) + crossprod(Magnetometer, w)
if(Ki > 0) eInt <- eInt + Merror * MSamplePeriod
else eInt <- c(0, 0, 0)

# Apply feedback terms
Gyroscope <- Gyroscope + Kp * Merror + Ki * eInt

# Compute rate of change of quaternion
qDot <- 0.5 * (q %Q*% cbind(0, Gyroscope[1], Gyroscope[2], Gyroscope[3]))
 
# Integrate to yield quaternion
q <- q + qDot * MSamplePeriod

MQuaternion <- Qnormalize(q)# normalise quaternion
return(c(MQuaternion)) #turn it into a vector
}

MahonykIMU<-function(MSamplePeriod, Kp=2.0, Ki=0.005, q, Gyroscope, Accelerometer)
{
#error management
if (any(is.na(c(MSamplePeriod, Kp, Ki, q, Gyroscope, Accelerometer)))) stop('Error! Input must be numeric.')
if (!is.null(dim(MSamplePeriod))) stop('Error! SamplePeriod must be a scalar.')
if (!is.null(dim(Kp))) stop('Error! Kp must be a scalar.')
if (!is.null(dim(Ki))) stop('Error! Ki must be a scalar.')
if (!is.null(dim(q))) stop('Error! q must be a vector.')
if (!is.null(dim(Gyroscope))) stop('Error! Gyroscope must be a vector.')
if (!is.null(dim(Accelerometer))) stop('Error! Accelerometer must be a vector.')
if (length(MSamplePeriod) !=1) stop('Error! SamplePeriod must be a scalar.')
if (length(Kp) !=1) stop('Error! Kp must be a scalar.')
if (length(Ki) !=1) stop('Error! Ki must be a scalar.')
if (length(q) !=4) stop('Error! q must be a vector with 4 elements.')
if (length(Gyroscope) !=3) stop('Error! Gyroscope must be a vector with 3 elements.')
if (length(Accelerometer) !=3) stop('Error! Accelerometer must be a vector with 3 elements.')
if (!is.numeric(MSamplePeriod))stop('Error! SamplePeriod must be numeric.')
if (!is.numeric(Kp))stop('Error! Kp must be numeric.')
if (!is.numeric(Ki))stop('Error! Ki must be numeric.')
if (!is.numeric(q))stop('Error! q must be numeric.')
if (!is.numeric(Gyroscope))stop('Error! Gyroscope must be numeric.')
if (!is.numeric(Accelerometer))stop('Error! Accelerometer must be numeric.')
eInt <- c(0, 0, 0) # integral error
# Normalise accelerometer measurement
if(Qnorm(Accelerometer) == 0) stop('Error! Accelerometer norm is zero.')
Accelerometer <- Qnormalize(Accelerometer)# normalise magnitude

# Estimated direction of gravity and magnetic field
v <- c(2*(q[2]*q[4] - q[1]*q[3]), 
 2*(q[1]*q[2] + q[3]*q[4]),
 q[1]^2 - q[2]^2 - q[3]^2 + q[4]^2)

# Error is sum of cross product between estimated direction and measured direction of fields
Merror <- crossprod(Accelerometer, v)
if(Ki > 0) eInt <- eInt + Merror * MSamplePeriod
else eInt <- c(0, 0, 0)

# Apply feedback terms
Gyroscope <- Gyroscope + Kp * Merror + Ki * eInt

# Compute rate of change of quaternion
qDot <- 0.5 * q %Q*% cbind(0, Gyroscope[1], Gyroscope[2], Gyroscope[3])
 
# Integrate to yield quaternion
q <- q + qDot * MSamplePeriod

MQuaternion <- Qnormalize(q)# normalise quaternion
return(c(MQuaternion)) #turn it into a vector
}

#	*	*	*	Function "MahonyAHRSupdateIMU" Package "RAHRS"		*	*	*
MahonyAHRSupdateIMU <- function(gxi,gyi,gzi,axi,ayi,azi,sampleFreqi,twoKpi,twoKii,integralFBxi,integralFByi,integralFBzi,q0i,q1i,q2i,halfex)
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(twoKpi)) stop('Argument <<twoKpi>> should be a number')
if (!is.numeric(twoKii)) stop('Argument <<twoKii>> should be a number')
if (!is.numeric(integralFBxi)) stop('Argument <<integralFBxi>> should be a number')
if (!is.numeric(integralFByi)) stop('Argument <<integralFByi>> should be a number')
if (!is.numeric(integralFBzi)) stop('Argument <<integralFBzi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(halfex)) stop('Argument <<halfex>> should be a number')
ret <- .C("MahonyAHRSupdateIMU",as.single(gxi),as.single(gyi),as.single(gzi),as.single(axi),as.single(ayi),as.single(azi),as.single(sampleFreqi),as.single(twoKpi),as.single(twoKii),as.single(integralFBxi),as.single(integralFByi),as.single(integralFBzi),as.single(q0i),as.single(q1i),as.single(q2i),as.single(halfex),DUP = TRUE, PACKAGE="RAHRS")
}

#	*	*	*	Function "MadgwickAHRSupdateIMU" Package "RAHRS"		*	*	*
MadgwickAHRSupdateIMU <- function(gxi,gyi,gzi,axi,ayi,azi,sampleFreqi,betai,q0i,q1i,q2i,sampleFreq)#twoKpDefi,twoKiDefi,
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(betai)) stop('Argument <<betai>> should be a number')
#if (!is.numeric(twoKpDefi)) stop('Argument <<twoKpDefi>> should be a number')
#if (!is.numeric(twoKiDefi)) stop('Argument <<twoKiDefi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(sampleFreq)) stop('Argument <<sampleFreq>> should be a number')
ret <- .C("MadgwickAHRSupdateIMU",as.single(gxi),as.single(gyi),as.single(gzi),as.single(axi),as.single(ayi),as.single(azi),as.single(sampleFreqi),as.single(betai),as.single(q0i),as.single(q1i),as.single(q2i),as.single(sampleFreq),DUP = TRUE, PACKAGE="RAHRS")#as.single(twoKpDefi),as.single(twoKiDefi),
}

#	*	*	*	Function "MadgwickAHRSupdate" Package "RAHRS"		*	*	*
MadgwickAHRSupdate <- function(gxi,gyi,gzi,axi,ayi,azi,mxi,myi,mzi,sampleFreqi,betai,q0i,q1i,q2i,gz)#twoKpDefi,twoKiDefi,
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(mxi)) stop('Argument <<mxi>> should be a number')
if (!is.numeric(myi)) stop('Argument <<myi>> should be a number')
if (!is.numeric(mzi)) stop('Argument <<mzi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(betai)) stop('Argument <<betai>> should be a number')
#if (!is.numeric(twoKpDefi)) stop('Argument <<twoKpDefi>> should be a number')
#if (!is.numeric(twoKiDefi)) stop('Argument <<twoKiDefi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(gz)) stop('Argument <<gz>> should be a number')
ret <- .C("MadgwickAHRSupdate",as.single(gxi),as.single(gyi),as.single(gzi),as.single(axi),as.single(ayi),as.single(azi),as.single(mxi),as.single(myi),as.single(mzi),as.single(sampleFreqi),as.single(betai),as.single(q0i),as.single(q1i),as.single(q2i),as.single(gz),DUP = TRUE, PACKAGE="RAHRS")#as.single(twoKpDefi),as.single(twoKiDefi),
}

#	*	*	*	Function "MahonyAHRSupdate" Package "RAHRS"		*	*	*
MahonyAHRSupdate <- function(gxi,gyi,gzi,axi,ayi,azi,mxi,myi,mzi,sampleFreqi,twoKpi,twoKii,integralFBxi,integralFByi,integralFBzi,q0i,q1i,q2i,halfex)
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(mxi)) stop('Argument <<mxi>> should be a number')
if (!is.numeric(myi)) stop('Argument <<myi>> should be a number')
if (!is.numeric(mzi)) stop('Argument <<mzi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(twoKpi)) stop('Argument <<twoKpi>> should be a number')
if (!is.numeric(twoKii)) stop('Argument <<twoKii>> should be a number')
if (!is.numeric(integralFBxi)) stop('Argument <<integralFBxi>> should be a number')
if (!is.numeric(integralFByi)) stop('Argument <<integralFByi>> should be a number')
if (!is.numeric(integralFBzi)) stop('Argument <<integralFBzi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(halfex)) stop('Argument <<halfex>> should be a number')
ret <- .C("MahonyAHRSupdate",as.single(gxi),as.single(gyi),as.single(gzi),as.single(axi),as.single(ayi),as.single(azi),as.single(mxi),as.single(myi),as.single(mzi),as.single(sampleFreqi),as.single(twoKpi),as.single(twoKii),as.single(integralFBxi),as.single(integralFByi),as.single(integralFBzi),as.single(q0i),as.single(q1i),as.single(q2i),as.single(halfex),DUP = TRUE, PACKAGE="RAHRS")
}

#	*	*	*	Function "MahonyAHRSupdateIMU2" Package "RAHRS"		*	*	*
MahonyAHRSupdateIMU2 <- function(gxi,gyi,gzi,axi,ayi,azi,sampleFreqi,twoKpi,twoKii,integralFBxi,integralFByi,integralFBzi,q0i,q1i,q2i,halfex)
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(twoKpi)) stop('Argument <<twoKpi>> should be a number')
if (!is.numeric(twoKii)) stop('Argument <<twoKii>> should be a number')
if (!is.numeric(integralFBxi)) stop('Argument <<integralFBxi>> should be a number')
if (!is.numeric(integralFByi)) stop('Argument <<integralFByi>> should be a number')
if (!is.numeric(integralFBzi)) stop('Argument <<integralFBzi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(halfex)) stop('Argument <<halfex>> should be a number')
ret <- .C("MahonyAHRSupdateIMU2",as.single(gxi),as.single(gyi),as.single(gzi),as.single(axi),as.single(ayi),as.single(azi),as.single(sampleFreqi),as.single(twoKpi),as.single(twoKii),as.single(integralFBxi),as.single(integralFByi),as.single(integralFBzi),as.single(q0i),as.single(q1i),as.single(q2i),as.single(halfex),DUP = TRUE, PACKAGE="RAHRS")
}

#	*	*	*	Function "MadgwickAHRSupdateIMU2" Package "RAHRS"		*	*	*
MadgwickAHRSupdateIMU2 <- function(gxi,gyi,gzi,axi,ayi,azi,sampleFreqi,betai,q0i,q1i,q2i,sampleFreq)#twoKpDefi,twoKiDefi,
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(betai)) stop('Argument <<betai>> should be a number')
#if (!is.numeric(twoKpDefi)) stop('Argument <<twoKpDefi>> should be a number')
#if (!is.numeric(twoKiDefi)) stop('Argument <<twoKiDefi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(sampleFreq)) stop('Argument <<sampleFreq>> should be a number')
ret <- .C("MadgwickAHRSupdateIMU2",as.single(gxi),as.single(gyi),as.single(gzi),as.single(axi),as.single(ayi),as.single(azi),as.single(sampleFreqi),as.single(betai),as.single(q0i),as.single(q1i),as.single(q2i),as.single(sampleFreq),DUP = TRUE, PACKAGE="RAHRS")#as.single(twoKpDefi),as.single(twoKiDefi),
}

#	*	*	*	Function "MadgwickAHRSupdate2" Package "RAHRS"		*	*	*
MadgwickAHRSupdate2 <- function(gxi,gyi,gzi,axi,ayi,azi,mxi,myi,mzi,sampleFreqi,betai,q0i,q1i,q2i,gz)#twoKpDefi,twoKiDefi,
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(mxi)) stop('Argument <<mxi>> should be a number')
if (!is.numeric(myi)) stop('Argument <<myi>> should be a number')
if (!is.numeric(mzi)) stop('Argument <<mzi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(betai)) stop('Argument <<betai>> should be a number')
#if (!is.numeric(twoKpDefi)) stop('Argument <<twoKpDefi>> should be a number')
#if (!is.numeric(twoKiDefi)) stop('Argument <<twoKiDefi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(gz)) stop('Argument <<gz>> should be a number')
ret <- .C("MadgwickAHRSupdate2",as.single(gxi),as.single(gyi),as.single(gzi),as.single(axi),as.single(ayi),as.single(azi),as.single(mxi),as.single(myi),as.single(mzi),as.single(sampleFreqi),as.single(betai),as.single(q0i),as.single(q1i),as.single(q2i),as.single(gz),DUP = TRUE, PACKAGE="RAHRS")#as.single(twoKpDefi),as.single(twoKiDefi),
}

#	*	*	*	Function "MahonyAHRSupdate2" Package "RAHRS"		*	*	*
MahonyAHRSupdate2 <- function(gxi,gyi,gzi,axi,ayi,azi,mxi,myi,mzi,sampleFreqi,twoKpi,twoKii,integralFBxi,integralFByi,integralFBzi,q0i,q1i,q2i,halfex)
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(mxi)) stop('Argument <<mxi>> should be a number')
if (!is.numeric(myi)) stop('Argument <<myi>> should be a number')
if (!is.numeric(mzi)) stop('Argument <<mzi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(twoKpi)) stop('Argument <<twoKpi>> should be a number')
if (!is.numeric(twoKii)) stop('Argument <<twoKii>> should be a number')
if (!is.numeric(integralFBxi)) stop('Argument <<integralFBxi>> should be a number')
if (!is.numeric(integralFByi)) stop('Argument <<integralFByi>> should be a number')
if (!is.numeric(integralFBzi)) stop('Argument <<integralFBzi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(halfex)) stop('Argument <<halfex>> should be a number')
ret <- .C("MahonyAHRSupdate2",as.single(gxi),as.single(gyi),as.single(gzi),as.single(axi),as.single(ayi),as.single(azi),as.single(mxi),as.single(myi),as.single(mzi),as.single(sampleFreqi),as.single(twoKpi),as.single(twoKii),as.single(integralFBxi),as.single(integralFByi),as.single(integralFBzi),as.single(q0i),as.single(q1i),as.single(q2i),as.single(halfex),DUP = TRUE, PACKAGE="RAHRS")
}

#	*	*	*	Function "MahonyAHRSupdateIMUDbl" Package "RAHRS"		*	*	*
MahonyAHRSupdateIMUDbl <- function(gxi,gyi,gzi,axi,ayi,azi,sampleFreqi,twoKpi,twoKii,integralFBxi,integralFByi,integralFBzi,q0i,q1i,q2i,halfex)
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(twoKpi)) stop('Argument <<twoKpi>> should be a number')
if (!is.numeric(twoKii)) stop('Argument <<twoKii>> should be a number')
if (!is.numeric(integralFBxi)) stop('Argument <<integralFBxi>> should be a number')
if (!is.numeric(integralFByi)) stop('Argument <<integralFByi>> should be a number')
if (!is.numeric(integralFBzi)) stop('Argument <<integralFBzi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(halfex)) stop('Argument <<halfex>> should be a number')
ret <- .C("MahonyAHRSupdateIMUDbl",as.numeric(gxi),as.numeric(gyi),as.numeric(gzi),as.numeric(axi),as.numeric(ayi),as.numeric(azi),as.numeric(sampleFreqi),as.numeric(twoKpi),as.numeric(twoKii),as.numeric(integralFBxi),as.numeric(integralFByi),as.numeric(integralFBzi),as.numeric(q0i),as.numeric(q1i),as.numeric(q2i),as.numeric(halfex),DUP = TRUE, PACKAGE="RAHRS")
}

#	*	*	*	Function "MadgwickAHRSupdateIMUDbl" Package "RAHRS"		*	*	*
MadgwickAHRSupdateIMUDbl <- function(gxi,gyi,gzi,axi,ayi,azi,sampleFreqi,betai,q0i,q1i,q2i,sampleFreq)#twoKpDefi,twoKiDefi,
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(betai)) stop('Argument <<betai>> should be a number')
#if (!is.numeric(twoKpDefi)) stop('Argument <<twoKpDefi>> should be a number')
#if (!is.numeric(twoKiDefi)) stop('Argument <<twoKiDefi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(sampleFreq)) stop('Argument <<sampleFreq>> should be a number')
ret <- .C("MadgwickAHRSupdateIMUDbl",as.numeric(gxi),as.numeric(gyi),as.numeric(gzi),as.numeric(axi),as.numeric(ayi),as.numeric(azi),as.numeric(sampleFreqi),as.numeric(betai),as.numeric(q0i),as.numeric(q1i),as.numeric(q2i),as.numeric(sampleFreq),DUP = TRUE, PACKAGE="RAHRS")#as.numeric(twoKpDefi),as.numeric(twoKiDefi),
}

#	*	*	*	Function "MadgwickAHRSupdateDbl" Package "RAHRS"		*	*	*
MadgwickAHRSupdateDbl <- function(gxi,gyi,gzi,axi,ayi,azi,mxi,myi,mzi,sampleFreqi,betai,q0i,q1i,q2i,gz)#twoKpDefi,twoKiDefi,
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(mxi)) stop('Argument <<mxi>> should be a number')
if (!is.numeric(myi)) stop('Argument <<myi>> should be a number')
if (!is.numeric(mzi)) stop('Argument <<mzi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(betai)) stop('Argument <<betai>> should be a number')
#if (!is.numeric(twoKpDefi)) stop('Argument <<twoKpDefi>> should be a number')
#if (!is.numeric(twoKiDefi)) stop('Argument <<twoKiDefi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(gz)) stop('Argument <<gz>> should be a number')
ret <- .C("MadgwickAHRSupdateDbl",as.numeric(gxi),as.numeric(gyi),as.numeric(gzi),as.numeric(axi),as.numeric(ayi),as.numeric(azi),as.numeric(mxi),as.numeric(myi),as.numeric(mzi),as.numeric(sampleFreqi),as.numeric(betai),as.numeric(q0i),as.numeric(q1i),as.numeric(q2i),as.numeric(gz),DUP = TRUE, PACKAGE="RAHRS")#as.numeric(twoKpDefi),as.numeric(twoKiDefi),
}

#	*	*	*	Function "MahonyAHRSupdateDbl" Package "RAHRS"		*	*	*
MahonyAHRSupdateDbl <- function(gxi,gyi,gzi,axi,ayi,azi,mxi,myi,mzi,sampleFreqi,twoKpi,twoKii,integralFBxi,integralFByi,integralFBzi,q0i,q1i,q2i,halfex)
{
if (!is.numeric(gxi)) stop('Argument <<gxi>> should be a number')
if (!is.numeric(gyi)) stop('Argument <<gyi>> should be a number')
if (!is.numeric(gzi)) stop('Argument <<gzi>> should be a number')
if (!is.numeric(axi)) stop('Argument <<axi>> should be a number')
if (!is.numeric(ayi)) stop('Argument <<ayi>> should be a number')
if (!is.numeric(azi)) stop('Argument <<azi>> should be a number')
if (!is.numeric(mxi)) stop('Argument <<mxi>> should be a number')
if (!is.numeric(myi)) stop('Argument <<myi>> should be a number')
if (!is.numeric(mzi)) stop('Argument <<mzi>> should be a number')
if (!is.numeric(sampleFreqi)) stop('Argument <<sampleFreqi>> should be a number')
if (!is.numeric(twoKpi)) stop('Argument <<twoKpi>> should be a number')
if (!is.numeric(twoKii)) stop('Argument <<twoKii>> should be a number')
if (!is.numeric(integralFBxi)) stop('Argument <<integralFBxi>> should be a number')
if (!is.numeric(integralFByi)) stop('Argument <<integralFByi>> should be a number')
if (!is.numeric(integralFBzi)) stop('Argument <<integralFBzi>> should be a number')
if (!is.numeric(q0i)) stop('Argument <<q0i>> should be a number')
if (!is.numeric(q1i)) stop('Argument <<q1i>> should be a number')
if (!is.numeric(q2i)) stop('Argument <<q2i>> should be a number')
if (!is.numeric(halfex)) stop('Argument <<halfex>> should be a number')
ret <- .C("MahonyAHRSupdateDbl",as.numeric(gxi),as.numeric(gyi),as.numeric(gzi),as.numeric(axi),as.numeric(ayi),as.numeric(azi),as.numeric(mxi),as.numeric(myi),as.numeric(mzi),as.numeric(sampleFreqi),as.numeric(twoKpi),as.numeric(twoKii),as.numeric(integralFBxi),as.numeric(integralFByi),as.numeric(integralFBzi),as.numeric(q0i),as.numeric(q1i),as.numeric(q2i),as.numeric(halfex),DUP = TRUE, PACKAGE="RAHRS")
}


