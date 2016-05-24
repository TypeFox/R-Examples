library(RAHRS)

# Data
data(anglesGyroLib)

accl.coefs <- matrix(c(-0.0771, -0.0817, -0.0769,  0.0066, -0.0032,  0.0001,  0.0090, -0.0023, -0.0035, -0.0281,  0.0430,  0.0306),ncol=1)
magn.coefs <- matrix(c(-0.3801, -0.3108, -0.3351,  0.0342,  0.0109, -0.0066,  0.0070, -0.0696, -0.0492,  0.0552,  0.0565,  0.0190),ncol=1)
gyro.coefs <- vector(mode='numeric',12)

# Number of simulation steps
Nsim <- dim(anglesGyroLib)[1]
m <- matrix(unlist(anglesGyroLib[,c('Mx','My','Mz')]),ncol=3,nrow=Nsim,byrow=FALSE)/3000
a <- matrix(unlist(anglesGyroLib[,c('Ax','Ay','Az')]),ncol=3,nrow=Nsim,byrow=FALSE)/3000
w <- matrix(unlist(anglesGyroLib[,c('Wx','Wy','Wz')]),ncol=3,nrow=Nsim,byrow=FALSE)/3000

#Number of simulation steps
Nsim <- dim(a)[1]

#Define Output data Arrays
R.out <- matrix(0,ncol=3,nrow=Nsim)
dw.out <- matrix(0,ncol=3,nrow=Nsim)
w.out <- matrix(0,ncol=3,nrow=Nsim)
a.out <- matrix(0,ncol=3,nrow=Nsim)
m.out <- matrix(0,ncol=3,nrow=Nsim)
TRIAD.out <- matrix(0,ncol=3,nrow=Nsim)

## Calibrate magnetometers
B <- matrix(as.numeric(c(magn.coefs[1], magn.coefs[4], magn.coefs[5],
                         magn.coefs[6], magn.coefs[2], magn.coefs[7],
                         magn.coefs[8], magn.coefs[9], magn.coefs[3])),nrow=3,ncol=3,byrow=TRUE)
B0 <- matrix(as.numeric(c(magn.coefs[10],magn.coefs[11],magn.coefs[12])),nrow=3,ncol=1)

for (n in 1:Nsim) m.out[n,] <- t((diag(3) - B) %*% (matrix(m[n,],3,1) - B0))

## Calibrate accelerometers
B <- matrix(as.numeric(c(accl.coefs[1], accl.coefs[4], accl.coefs[5],
                         accl.coefs[6], accl.coefs[2], accl.coefs[7],
                         accl.coefs[8], accl.coefs[9], accl.coefs[3])),nrow=3,ncol=3,byrow=TRUE)
B0 <- matrix(as.numeric(c(accl.coefs[10],accl.coefs[11],accl.coefs[12])),nrow=3,ncol=1)

for (n in 1:Nsim) a.out[n,] <- t((diag(3)-B) %*% (matrix(a[n,],3,1) - B0))

## Calibrate gyroscopes
B <- matrix(as.numeric(c(gyro.coefs[1], gyro.coefs[4], gyro.coefs[5],
                         gyro.coefs[6], gyro.coefs[2], gyro.coefs[7],
                         gyro.coefs[8], gyro.coefs[9], gyro.coefs[3])),nrow=3,ncol=3,byrow=TRUE)
B0 <- matrix(as.numeric(c(gyro.coefs[10],gyro.coefs[11],gyro.coefs[12])),nrow=3,ncol=1)

for (n in 1:Nsim) w.out[n,] <- t((diag(3)-B) %*% (matrix(w[n,],3,1) - B0))

## AHRS Parameters
#Magnetic Field Vector In Navigation Frame
Parameters<-list(mn=0,an=0,dt=0)
Parameters$mn = matrix(c(0.315777529635464, 0.057133095826051, -0.947111588535720),nrow=1,ncol=3)
#Acceleration vector In Navigation Frame
Parameters$an = matrix(c(0, 0, -1),nrow=1,ncol=3)
#Sampling Rate 1/Hz
Parameters$dt  = c(1/100)
#Initial attitude quaternion value
q = c( 1.0,0.0,0.0,0.0 ) 
#initial value of estimated bias
dw.outhat = matrix(0,ncol=1,nrow=3)

#Filter parameters and states
Filter <- list(P=0,Q=0,R=0)
Filter$P = diag(c(1e-3, 1e-3, 1e-3, 1e-3, 1e-3, 1e-3))
Filter$Q = diag(c(1e-4, 1e-4, 1e-4, 1e-10, 1e-10, 1e-10))
Filter$R = diag(c(1e-1, 1e-1, 1e-1, 1e-1, 1e-1, 1e-1))
Sensors<-list(w=0,a=0,m=0)

## Main loop
for (n in 1:Nsim)#Nsim
    {#n<-100
    Sensors$w = (matrix(w.out[n,],nrow=1,ncol=3))
    Sensors$a = (matrix(a.out[n,],nrow=1,ncol=3))
    Sensors$m = (matrix(m.out[n,],nrow=1,ncol=3))
    
    ## LKF Part
tmp <- ahrs.LKF.QUATERNION(Filter, Sensors, q, Parameters, dw.outhat)
Filter <- tmp$Filter
q <- tmp$q
dw.outhat <- tmp$dw
R.out[n,] <-Q2EA(q)
dw.out[n,] <- dw.outhat

    ## TRIAD Part
    W1 = matrix(a.out[n,],nrow=1) / norm(matrix(a.out[n,],nrow=1),'f') 
    W2 = matrix(m.out[n,],nrow=1)/norm(matrix(m.out[n,],nrow=1),'f')
    
    V1 = Parameters$an
    V2 = Parameters$mn
    
    Ou1 = W1
    Ou2 = cross(W1,W2)/norm(cross(W1,W2),'f')
    Ou3 = cross((W1),cross(W1,W2))/norm(cross(W1,W2),'f')

    R1 = V1
    R2 = (cross(V1,V2)/norm(cross(V1,V2),'f'))
    R3 = (cross((V1),cross(V1,V2))/norm(cross(V1,V2),'f'))
    
    Mou = cbind(t(Ou1), t(Ou2), t(Ou3))
    Mr = cbind(t(R1), t(R2), t(R3))
    
    A = Mou %*% t(Mr)
    
    # Calculate angles
    TRIAD.out[n,] <- DCM2EA(A)#rotMat2euler(A)
   list(w.out=w.out,a.out=a.out,m.out=m.out,R.out=R.out,dw.out=dw.out, TRIAD.out=TRIAD.out)
}


#postscript('LKF.QUATERNION.pdf')
## Plot Results
#psi - Angle around Z axis
#theta - Angle around Y axis
#gamma - Angle around X axis

plot(TRIAD.out[,1],col='red',ylim=c(-4,4), type='l',main='LKF quaternion TRIAD / AHRS',xlab='Time 10 msec (100 Hz)', ylab='Euler angles')
par(new=T)
plot(TRIAD.out[,2],col='blue',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
par(new=T)
plot(TRIAD.out[,3],col='green',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
par(new=T)
plot(R.out[,1],col='orange',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
par(new=T)
plot(R.out[,2],col='cyan',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
par(new=T)
plot(R.out[,3],col='magenta',ylim=c(-4,4), type='l',main='',xlab='', ylab='')
abline(h=(seq(-4,4,1)), col="lightgray", lty="dotted")
abline(v=(seq(0,16000,1000)), col="lightgray", lty="dotted")
legend("topleft", c( expression(paste(psi,plain(TRIAD))) ,expression(paste(theta,plain(TRIAD))),
expression(paste(gamma,plain(TRIAD))),expression(paste(psi,plain(AHRS))) ,expression(paste(theta,plain(AHRS))),
expression(paste(gamma,plain(AHRS)))),col=c('red','blue','green','orange','cyan','magenta'), lty = c(1, 1, 1, 1, 1, 1),bg='white')

#dev.off()

