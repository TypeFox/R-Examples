

Residuals.parabola <- function(XY,ParG)
{
#   Projecting a given set of points onto a parabola and computing the distances from the points to the parabola
Vertex <- ParG[1:2]
p= ParG[3]  
Angle <- ParG[4]
n <- dim(XY)[1]
XYproj <- matrix(0,n,2)
Iter <- matrix(0,n,1)
tolerance <- 1e-9
tol.p <- tolerance * p/2
pp <- p*p
#  Matrix Q for rotating the points and the parabola
s <- sin(Angle) 
cA <- cos(Angle)
Q <- matrix(c(cA, -s,s, cA),2,2,byrow=TRUE)
#  data points in canonical coordinates
XY0  <- cbind(XY[,1]-Vertex[1], XY[,2]-Vertex[2]) %*% Q   
XYA <- cbind(XY0[,1], abs(XY0[,2]))
#  main loop over the data points
for (i in 1:n){
    u <- XYA[i,1]  
v <- XYA[i,2]
    uu <- u*u      
vv <- v*v   
pu=p*u
    if (v == 0) z2 <- 1 else z2 <- sign(XY0[i,2])
    #       does the point lie on the x-axis?
    if (v<tol.p){
        if (u>p) XYproj[i,] <- cbind(u-p, z2*sqrt(max(2*p*(u-p),0))) else XYproj[i,] <- c(0, 0)
        next
    }
    #  pick the initial point T starting from -1/2 until F(t)>0
    for (j in 1:20){
        Tp <- 1/2^j-1  
        Fp <- vv/(Tp+1)^2 - 2*pu-2*pp*Tp
        if (Fp>0) break
    }
    #      generic case: start the iterative procedure
    for (iter in 1:100){
        F0= vv/(Tp+1)^2
        Fp  <- F0-2*(pu+pp*Tp)
        if (Fp<0) break
        Fder <- 2*(F0/(Tp+1) + pp)
        Ratio <- Fp/Fder
        if (Ratio<tol.p){
Iter[i]=iter
break
}
        Tp <- Tp + Ratio
    }
    #      compute the projection of the point onto the parabola
    XYproj[i,] <- cbind(u+Tp*p, XY0[i,2]/(Tp+1))
} # end the main loop
XYproj <- XYproj %*% t(Q)
XYproj <- cbind(XYproj[,1]+Vertex[1], XYproj[,2]+Vertex[2])
RSS <- norm(XY-XYproj,'F')^2
list(RSS=RSS, XYproj=XYproj)
}   # Residuals.parabola

Residuals.hyperbola <- function(XY,ParG)
{#   Projecting a given set of points onto a hyperbola and computing the distances from the points to the hyperbola
Center <- ParG[1:2]   
Axes <- ParG[3:4]  
Angle <- ParG[5]
n <- dim(XY)[1]
XYproj <- matrix(0,n,2)
tolerance <- 1e-9
a <- Axes[1]
b=Axes[2]
aa <- a^2  
bb <- b^2 
at=sqrt(a)
bt=sqrt(b)
tol.a <- tolerance*a
tol.b <- tolerance*b
tol.aa <- tolerance*aa
#  Matrix Q for rotating the points and the hyperbola
s <- sin(Angle)
cA <- cos(Angle)
Q <- matrix(c(cA, -s,s, cA),2,2,byrow=TRUE)
#  data points in canonical coordinates
XY0  <- (cbind(XY[,1]-Center[1], XY[,2]-Center[2])) %*% Q
XYA <- abs(XY0)
# find the inflection points TInf for all given pairs
XYS <- sqrt(XYA)
XYS1 <- matrix(bb*at*XYS[,1]-aa*bt*XYS[,2],ncol=1)
XYS2 <- matrix(at*XYS[,1]+bt*XYS[,2],ncol=1)
TInf <- matrix(XYS1/XYS2,ncol=1) #inflection points
#  main loop over the data points
for (i in 1:n){
    u <- XYA[i,1]  
v <- XYA[i,2]
    ua <- u*a      
vb <- v*b
    if (u == 0) z1 <- 1 else z1 <- sign(XY0[i,1])
    if (v == 0) z2 <- 1 else z2 <- sign(XY0[i,2])
    #       does the point lie on the major axis?
    if (v<tol.b){
        if (u>a+bb/a){
            xproj <- aa*u/(aa+bb)
            XYproj[i,] <- cbind(z1*xproj, z2*b*sqrt(max((xproj/a)^2-1,0)))
       } else XYproj[i,] <- cbind(z1*a, 0)
        next
    } # end if
    #       does the point lie on the minor axis?
    if (u<tol.a){
        yproj <- bb*v/(aa+bb)
        XYproj[i,] <- cbind(z1*a*sqrt(1+(yproj/b)^2), z2*yproj )
        next
    } # } if
    #     generic case: start the iterative procedure
    T0 <- TInf[i]  # inflection point t
    F0  <- (ua/(T0+aa))^2 - (vb/(-T0+bb))^2 - 1  # the corresponding F(t)
    # if F0>0 then pick the initial point to the right of the root
    # if F0<0 then pick the initial point to the left of the root
    if (F0 > 0){  # case1: pick the initial point T until F(T)<0
        for (j in 1:20){
            Tp<-bb-(-T0+bb)/2^j
            Fp<-(ua/(Tp+aa))^2 - (vb/(-Tp+bb))^2 - 1
            if (Fp<0) break
        } # end for j
        for (iter in 1:100){
            Taa <- Tp + aa
            Tbb <- -Tp + bb
            PP1 <- (ua/Taa)^2
            PP2 <- (vb/Tbb)^2
            Fp  <- PP1 - PP2 - 1
            if (Fp>0) break
            Fder <- 2*(PP1/Taa + PP2/Tbb)
            Ratio <- Fp/Fder
            if (abs(Ratio)<tol.aa) break
            Tp <- Tp + Ratio
        } } else { for (j in 1:20) { # case2: pick the initial point T until F(T)>0
            Tp=-aa+(T0+aa)/2^j
            Fp=(ua/(Tp+aa))^2 - (vb/(-Tp+bb))^2 - 1
            if (Fp>0) break
        } # end for j
        for (iter in 1:100){
            Taa <- Tp + aa
            Tbb <- -Tp + bb
            PP1 <- (ua/Taa)^2
            PP2 <- (vb/Tbb)^2
            Fp  <- PP1 - PP2 - 1
            if (Fp<0)  break
            Fder <- 2*(PP1/Taa + PP2/Tbb)
            Ratio <- Fp/Fder
            if (abs(Ratio)<tol.aa) break
            Tp <- Tp + Ratio
        }
    } # end if
    #   compute the projection of the point onto the hyperbola
    if (Taa < 1e-6){
        yproj <- XY0[i,2]*bb/Tbb
        XYproj[i,] <- cbind(sign(XY0[i,1])*a*sqrt(1+(yproj/b)^2), yproj)
    } else { if (Tbb < 1e-6){
        xproj=XY0[i,1]*aa/Taa
        yproj <- sign(XY0[i,2])*b*sqrt(max((xproj/a)^2-1,0))
        XYproj[i,] <- cbind(xproj, yproj)
   } else XYproj[i,] <- cbind(XY0[i,1]*aa/Taa, XY0[i,2]*bb/Tbb)
    }  # end if
} # } the main for-loop
XYproj <- XYproj %*% t(Q)
XYproj <- cbind(XYproj[,1]+Center[1], XYproj[,2]+Center[2])
RSS <- norm(XY-XYproj,'F')^2
list(RSS=RSS, XYproj=XYproj)
}   # Residuals.hyperbola


fit.conicLMA <- function(XY,ParAini,LambdaIni=1,epsilonP=0.0000000001,epsilonF=0.0000000000001,IterMAX=2000000)
{
# Fitting a conic to a given set of points (Implicit method) using algebraic parameter
#epsilonP <- tolerance (small threshold)
#epsilonF <- tolerance (small threshold)
#IterMAX <- maximal number of (main) iterations usually 10-20 suffice
lambda.sqrt <- sqrt(LambdaIni)   #  sqrt(Lambda) is actually used by the code
tmp <- AtoG(ParAini)
ParGini <- tmp$ParG
exitCode <- tmp$exitCode
if (exitCode == 1){
    tmp <- Residuals.ellipse(XY,ParGini)
  Fp <- tmp$RSS
XYproj <- tmp$XYproj  
} else { if (exitCode == 2){
tmp <- Residuals.hyperbola(XY,ParGini)
Fp <- tmp$RSS
XYproj <- tmp$J
} else stop('invalid initial parameter')
}
tmp <- JmatrixLMA(XY,ParAini,XYproj)  # calculate the Jacobian matrix
  Res <- tmp$Res
  J <- tmp$J
# tmp <- HessianLMA(XY,ParAini,XYproj) # calculate the Hessian matrix
#  H <- tmp$H
#  ev <- tmp$ev
ParA <- ParAini
ParG <- ParGini
codeTemp <- exitCode
ParGTemp <- ParG
for (iter in 1:IterMAX)         #  main loop, each run is one (main) iteration
{
    while (TRUE)         #  secondary loop - adjusting Lambda (no limit on cycles)
{
      DelPar <- mldivide(rbind(J, lambda.sqrt*diag(6)), rbind(-Res, matrix(0,6,1))) # step candidate
        ParTemp <- ParA + DelPar
        ParTemp <- ParTemp/norm(ParTemp,'F')
        tmp <- AtoG(ParTemp)
        ParGTemp <- tmp$ParG
        codeTemp <- tmp$exitCode
        if (codeTemp != exitCode) progress <- 1 else progress <- norm(DelPar,'F')
        if (progress < epsilonP)  break               # stopping rule
        if (codeTemp == 1){
            tmp <- Residuals.ellipse(XY,ParGTemp)
      FTemp <- tmp$RSS
      XYprojTemp <- tmp$XYproj
        } else {if (codeTemp == 2){
            tmp <- Residuals.hyperbola(XY,ParGTemp)
      FTemp <- tmp$RSS
      XYprojTemp <- tmp$XYproj
        }else{
            lambda.sqrt <- lambda.sqrt*2 #if it's degenerate, increase lambda, recompute the step
            next
            }
        }
        if (FTemp < Fp*(1-epsilonF/lambda.sqrt)){        #   yes, improvement
            lambda.sqrt <- lambda.sqrt/2   # reduce lambda, move to next iteration
            break
        }else {                           #   no improvement
            lambda.sqrt <- lambda.sqrt*2 # increase lambda, recompute the step
            next
        }
    }   #   while (TRUE), the end of the secondary loop
    if (progress < epsilonP)  break # stopping rule
    tmp <- JmatrixLMA(XY,ParTemp,XYprojTemp)
  Res <- tmp$Res
  J <- tmp$J
    ParA <- ParTemp
 Fp <- FTemp  # update the iteration
    ParG <- ParGTemp
exitCode <- codeTemp
}      #    main loop
RSS <- Fp
iters <- iter
list(ParA=ParA,RSS=RSS,iters=iters,exitCode=exitCode)
}

GtoA<-function(ParG)
{# Conversion of geometric parameters of an ellipse
# ParG <- [Center(1:2), Axes(1:2), Angle]'
# to algebraic parameters
k <- cos(ParG[5])
s <- sin(ParG[5])
a <- ParG[3]
b <- ParG[4]
Xc <- ParG[1]
Yc <- ParG[2]
P <- (k/a)^2 + (s/b)^2
Q <- (s/a)^2 + (k/b)^2
R <- 2*k*s*(1/a^2 - 1/b^2)
cbind(A=P, B=R, C=Q, D=-2*P*Xc-R*Yc, E=-2*Q*Yc-R*Xc, F=P*Xc^2+Q*Yc^2+R*Xc*Yc-1)
}

AtoG <- function(ParA)
{
#  Conversion of algebraic parameters to geometric parameters
tolerance1 <- 1.e-10
tolerance2 <- 1.e-20
if (abs(ParA[1]-ParA[3]) > tolerance1) Angle <- atan(ParA[2]/(ParA[1]-ParA[3]))/2 else Angle <- pi/4
cA <- cos(Angle)
s <- sin(Angle)
Q <- matrix(c(cA, s,-s, cA),2,2,byrow=TRUE)
M <- matrix(c(ParA[1], ParA[2]/2, ParA[2]/2, ParA[3]),2,2,byrow=TRUE)
D <- Q %*% M %*% t(Q)
N <- Q %*% matrix(c(ParA[4], ParA[5]),2,1,byrow=TRUE)
O <- ParA[6]
if ((D[1,1] < 0) && (D[2,2] < 0)){
    D <- -D
    N <- -N
    O <- -O
}
UVcenter <- matrix(c(-N[1,1]/2/D[1,1], -N[2,1]/2/D[2,2]),2,1,byrow=TRUE)
free <- O - UVcenter[1,1]*UVcenter[1,1]*D[1,1] - UVcenter[2,1]*UVcenter[2,1]*D[2,2]
# if the determinant of [A B/2 D/2 B/2 C E/2 D/2 E/2 F]is zero 
# And if K>0,then it's a empty set
# otherwise the conic is degenerate
Deg <- matrix(c(ParA[1],ParA[2]/2,ParA[4]/2, ParA[2]/2,ParA[3],ParA[5]/2, ParA[4]/2,ParA[5]/2,ParA[6]),3,3,byrow=TRUE)
K1 <- matrix(c(ParA[1],ParA[4]/2, ParA[4]/2, ParA[6]),2,2,byrow=TRUE)
K2 <- matrix(c(ParA[3],ParA[5]/2, ParA[5]/2, ParA[6]),2,2,byrow=TRUE)
K <- det(K1)+det(K2)
if ((abs(det(Deg)) < tolerance2)) if ((abs(det(M))<tolerance2) && (K > tolerance2)) { exitCode <- 4 
# empty set(imaginary parellel lines)
   } else {
        exitCode <- -1
# degenerate cases
    } else {
    if (D[1,1]*D[2,2] > tolerance1) if (free < 0) { exitCode <- 1
 # ellipse
       } else {
            exitCode <- 0
 # empty set(imaginary ellipse)
        } else { if (D[1,1]*D[2,2] < - tolerance1){
        exitCode <- 2
  # hyperbola
   } else {
        exitCode <- 3
  # parabola
    }
} }
XYcenter <- t(Q) %*% UVcenter
Axes <- matrix(c(sqrt(abs(free/D[1,1])), sqrt(abs(free/D[2,2]))),2,1)
if (exitCode == 1 && Axes[1]<Axes[2]) {
AA <- Axes[1]
Axes[1] <- Axes[2]
Axes[2] <- AA
Angle <- Angle + pi/2
}
if (exitCode == 2 && free*D[1,1]>0) {
AA <- Axes[1]
Axes[1] <- Axes[2]
Axes[2] <- AA
Angle <- Angle + pi/2
}
while (Angle > pi) Angle <- Angle - pi
while (Angle < 0) Angle <- Angle + pi
ParG <- rbind(XYcenter, Axes, Angle)
dimnames(ParG) <- NULL
list(ParG=ParG, exitCode=exitCode)
}

Residuals.ellipse <- function(XY,ParG)
{
Center <- ParG[1:2]
Axes <- ParG[3:4]
Angle <- ParG[5]
n <- dim(XY)[1]
XYproj <- matrix(0,n,2)
tolerance <- 1e-9
#  First handling the circle case
if (abs((Axes[1]-Axes[2])/Axes[1])<tolerance){
    phiall <- angle(XY[,1]-Center[1] + sqrt(-1) * (XY[,2]-Center[2]))
    XYproj <- matrix(c(Axes[1] * cos(phiall)+Center[1], Axes[2] * sin(phiall)+Center[2]),1,2)
    RSS <- norm(XY-XYproj,'F')^2
    list(RSS=RSS, XYproj=XYproj)
}
# Now dealing with proper ellipses
a <- Axes[1]
b <- Axes[2]
aa <- a^2
bb <- b^2
tol.a <- tolerance * a
tol.b <- tolerance * b
tol.aa <- tolerance * aa
#  Matrix Q for rotating the points and the ellipse to the canonical system
s <- sin(Angle)
cA <- cos(Angle)
Q <- matrix(c(cA, -s,s, cA),2,2,byrow=TRUE)
#  data points in canonical coordinates
XY0  <- cbind(XY[,1]-Center[1], XY[,2]-Center[2]) %*% Q
XYA <- abs(XY0)
Tini <- apply(cbind(a * (XYA[,1]-a),b * (XYA[,2]-b)),1,max)
#  main loop over the data points
for (i in 1:n){
    u <- XYA[i,1]
 v <- XYA[i,2]
    ua <- u * a
     vb <- v * b
    if (u == 0) z1 <- 1 else z1 <- sign(XY0[i,1])
    if (v == 0) z2 <- 1 else z2 <- sign(XY0[i,2])
    #       does the point lie on the minor axis?
    if (u<tol.a){
        if (XY0[i,2]<0) XYproj[i,] <- cbind(0, -b) else XYproj[i,] <- cbind(0, b)
        next
    }
    #       does the point lie on the major axis?
    if (v<tol.b){
        if (u < (a-bb/a)){
            xproj <- aa * u/(aa-bb)
            XYproj[i,] <- matrix(c(z1 * xproj, z2 * b * sqrt(max(1-(xproj/a)^2,0))),1,2)
        } else {
            XYproj[i,] <- cbind(z1 * a, 0)
        }
        next
    }
    #      generic case: start the iterative procedure
    T <- Tini[i]
    for (iter in 1:100){
        Taa <- T + aa
        Tbb <- T + bb
        PP1 <- (ua/Taa)^2
        PP2 <- (vb/Tbb)^2
        Fp  <- PP1 + PP2 - 1
        if (Fp<0) break
        Fder <- 2 * (PP1/Taa + PP2/Tbb)
        Ratio <- Fp/Fder
        if (Ratio<tol.aa) break
        T <- T + Ratio
    }
    #       compute the projection of the point onto the ellipse
    xproj <- XY0[i,1] * aa/Taa
    yproj <- sign(XY0[i,2]) * b * sqrt(max(1-(xproj/a)^2,0))
    XYproj[i,] <- cbind(xproj, yproj)
} # end the main loop
#    rotate back to the original system
XYproj <- XYproj %*% t(Q)
XYproj <- cbind(XYproj[,1]+Center[1], XYproj[,2]+Center[2])
RSS <- norm(XY-XYproj,'F')^2
list(RSS=RSS,XYproj=XYproj)
}   # Residuals.ellipse

JmatrixLMA<-function(XY,ParA,XYproj){
#Compute the Jacobian matrix(Implicit method)using algebraic parameter
n <- dim(XY)[1]
Res <- matrix(0,n,1)
X <- matrix(XY[,1],ncol=1)
Y <- matrix(XY[,2],ncol=1)
Z <- cbind(X^2, X*Y, Y^2, X, Y, matrix(1,n,ncol=1) ) %*% ParA
DD <- XY-XYproj
for (i in 1:n) Res[i]  <- sign(Z[i])*norm(matrix(DD[i,]),'F')
D2 <- matrix(c(XYproj,matrix(1,n,1)),ncol=3,byrow=FALSE)
x <- matrix(XYproj[,1],ncol=1)
y <- matrix(XYproj[,2],ncol=1)
xx <- matrix(x^2,ncol=1)
yy <- matrix(y^2,ncol=1)
xy <- matrix(x*y,ncol=1)
# dPar <- matrix(c(xx,xy,yy,x,y,ones(n,1)),1)       #partial derivative of P wrt ParA
du <- D2 %*% matrix(c(2*ParA[1],ParA[2],ParA[4]),3)   #partial derivative of P wrt x-coordinate
dv <- D2 %*% matrix(c(ParA[2],2*ParA[3],ParA[5]),3)   #partial derivative of P wrt y-coordinate
eA  <- sqrt(du^2+dv^2)
J <- cbind(xx/eA, xy/eA, yy/eA, x/eA, y/eA,matrix(1,n,1)/eA)
list(Res = Res,J = J)
}

ResidualsG<-function(XY,ParG)
{
#   Projecting a given set of points onto an ellipse and computing the distances from the points to the ellipse
Center  <-  ParG[1:2]
Axes  <-  ParG[3:4]
Angle  <-  ParG[5]
n  <-  dim(XY)[1]
XYproj  <-  matrix(0,n,2)
tolerance  <-  1e-9
#  First handling the circle case
if (abs((Axes[1]-Axes[2])/Axes[1])<tolerance){
   phiall  <-  angle(XY[,1]-Center[1] + sqrt(-1)*(XY[,2]-Center[2]))
   XYproj  <-  cbind(Axes[1]*cos(phiall)+Center[1], Axes[2]*sin(phiall)+Center[2])
   RSS  <-  norm(XY-XYproj,'F')^2
   list(RSS=RSS, XYproj=XYproj)
}
#  Now dealing with proper ellipses
a  <-  Axes[1]
b  <-  Axes[2]
aa  <-  a^2
bb  <-  b^2
tol_a  <-  tolerance*a
tol_b  <-  tolerance*b
tol_aa  <-  tolerance*aa
#  Matrix Q for rotating the points and the ellipse
s <- sin(Angle)
cA <- cos(Angle)
Q <- matrix(c(cA, -s,s, cA),2,2,byrow=TRUE)
XY0   <-  cbind(XY[,1]-Center[1], XY[,2]-Center[2]) %*% Q
   #  data points in canonical coordinates
XYA  <-  abs(XY0)
Tini  <- matrix( apply(cbind(a*(XYA[,1]-a),b*(XYA[,2]-b)),1,max),ncol=1)
#  main loop over the data points
for (i in 1:n){
    u  <-  XYA[i,1]
 v  <-  XYA[i,2]
    ua  <-  u * a
     vb  <-  v * b
#       does the point lie on the minor axis?
    if (u<tol_a){
       if (XY0[i,2]<0) XYproj[i,]  <-  cbind(0, -b) else XYproj[i,]  <-  cbind(0, b)
       next
    }
#       does the point lie on the major axis?
    if (v<tol_b){
       if (u<a-bb/a){
          xproj  <-  aa*u/(aa-bb)
          XYproj[i,]  <-  cbind(sign(XY0[i,1])*xproj, sign(XY0[i,2])*b*sqrt(1-(xproj/a)^2))
      } else {
          XYproj[i,]  <-  cbind(sign(XY0[i,1])*a, 0)
       }
       next
    }
#      generic case: start the iterative procedure
    T  <-  Tini[i]
    for (iter in 1:100){
        Taa  <-  T + aa
        Tbb  <-  T + bb
        PP1  <-  (ua/Taa)^2
        PP2  <-  (vb/Tbb)^2
        Fp  <-  PP1 + PP2 - 1
        if (Fp<0) break
        Fder  <-  2*(PP1/Taa + PP2/Tbb)
        Ratio  <-  Fp/Fder
        if (Ratio<tol_aa) break
        T  <-  T + Ratio
    }
#      compute the projection of the point onto the ellipse
    XYproj[i,]  <-  cbind(XY0[i,1]*aa/Taa, XY0[i,2]*bb/Tbb)
}
XYproj  <-  XYproj %*% t(Q)
XYproj  <-  cbind(XYproj[,1]+Center[1], XYproj[,2]+Center[2])
RSS  <-  norm(XY-XYproj,'F')^2
list(RSS=RSS, XYproj=XYproj)
}   # ResidualsG

fit.ellipseLMG <- function(XY,ParGini,LambdaIni=1,epsilon=0.000001,IterMAX = 200,L = 200)
{# fitting ellipse using Implicit method
  #epsilon <- tolerance (small threshold)
  #IterMAX <- maximal number of (main) iterations usually 10-20 suffice
  #L <- boundary for major/minor axis
  lambda.sqrt <- sqrt(LambdaIni)#  sqrt(Lambda) is actually used by the code
  TF <- FALSE
  tmp <- ResidualsG(XY,ParGini)
  Fp <- tmp$RSS
  XYproj <- tmp$XYproj
  tmp <- JmatrixLMG(XY,ParGini,XYproj)
  Res <- tmp$Res
  J <- tmp$J
  ParG <- ParGini
  for (iter in 1:IterMAX)         #  main loop, each run is one (main) iteration
  {
    while (TRUE)         #  secondary loop - adjusting Lambda (no limit on cycles)
    {
      DelPar <- mldivide(rbind(J, lambda.sqrt*diag(5)), rbind(-Res, matrix(0,5,1))) # step candidate
      progress <- norm(DelPar,'F')/(norm(ParG,'F')+epsilon)
      if (progress < epsilon)  break            # stopping rule
      ParTemp <- ParG + DelPar
      if (ParTemp[3]< ParTemp[4]){        #   out of range
        Temp=ParTemp[3]
        ParTemp[3]=ParTemp[4]
        ParTemp[4]=Temp
        ParTemp[5]=ParTemp[5]-sign(ParTemp[5])*pi/2
      }
      tmp  <- ResidualsG(XY,ParTemp)
      FTemp <- tmp$RSS
      XYprojTemp <- tmp$XYproj
      if (FTemp < Fp && ParTemp[3]>0 && ParTemp[4]>0) {       #   yes, improvement
        lambda.sqrt <- lambda.sqrt/2   # reduce lambda, move to next iteration
        break
      } else {                           #   no improvement
        lambda.sqrt <- lambda.sqrt*2 # increase lambda, recompute the step
        next
      }
    }   #   while TRUE, the end of the secondary loop
    if (ParTemp[3] > L){       # diverge
      TF <- TRUE
      break
    }
    if (progress < epsilon) break  
    tmp <- JmatrixLMG(XY,ParTemp,XYprojTemp)
    Res <- tmp$Res
    J <- tmp$J
    ParG <- ParTemp
    Fp <- FTemp  # update the iteration
  }
  RSS <- Fp
  iters <- iter
  # make the angle parameter between 0 and pi
  while(ParG[5] >= pi) ParG[5] <- ParG[5] - pi
  while(ParG[5] < 0) ParG[5] <- ParG[5] + pi
  if (iters >= IterMAX) TF <- TRUE
  list(ParG=ParG,RSS=RSS,iters=iters,TF=TF)
}   #fit.ellipseLMG

fit.ellipseLMG.H <- function(XY,ParGini,LambdaIni=1,epsilon=0.000001,IterMAX = 200,L = 200)
{# fitting ellipse using Implicit method
  #epsilon <- tolerance (small threshold)
  #IterMAX <- maximal number of (main) iterations usually 10-20 suffice
  #L <- boundary for major/minor axis
  lambda.sqrt <- sqrt(LambdaIni)# sqrt(Lambda) is actually used by the code
  TF <- FALSE
  tmp <- ResidualsG(XY,ParGini)
  Fp <- tmp$RSS
  XYproj <- tmp$XYproj
  tmp <- JmatrixLMG(XY,ParGini,XYproj)
  Res <- tmp$Res
  J <- tmp$J
  ParG <- ParGini
  for (iter in 1:IterMAX)         #  main loop, each run is one (main) iteration
  {
    while (TRUE)         #  secondary loop - adjusting Lambda (no limit on cycles)
    {
      DelPar <- mldivide(rbind(J, lambda.sqrt*diag(5)), rbind(-Res, matrix(0,5,1))) # step candidate
      progress <- norm(DelPar,'F')/(norm(ParG,'F')+epsilon)
      if (progress < epsilon)  break            # stopping rule
      ParTemp <- ParG + DelPar
      if (ParTemp[3]< ParTemp[4]){        #   out of range
        Temp=ParTemp[3]
        ParTemp[3]=ParTemp[4]
        ParTemp[4]=Temp
        ParTemp[5]=ParTemp[5]-sign(ParTemp[5])*pi/2
      }
      tmp  <- ResidualsG(XY,ParTemp)
      FTemp <- tmp$RSS
      XYprojTemp <- tmp$XYproj
      if (FTemp < Fp && ParTemp[3]>0 && ParTemp[4]>0) {       #   yes, improvement
        lambda.sqrt <- lambda.sqrt/2   # reduce lambda, move to next iteration
        break
      } else {                           #   no improvement
        lambda.sqrt <- lambda.sqrt*2 # increase lambda, recompute the step
        next
      }
    }   #   while TRUE, the end of the secondary loop
    if (ParTemp[3] > L){       # diverge
      TF <- TRUE
      break
    }
    if (progress < epsilon) break  
    tmp <- JmatrixLMG(XY,ParTemp,XYprojTemp)
    Res <- tmp$Res
    J <- tmp$J
    ParG <- ParTemp
    Fp <- FTemp  # update the iteration
  }
  RSS <- Fp
  iters <- iter
  Jg <- t(J) %*% Res
  H <- t(J) %*% J
  # make the angle parameter between 0 and pi
  while(ParG[5] >= pi) ParG[5] <- ParG[5] - pi
  while(ParG[5] < 0) ParG[5] <- ParG[5] + pi
  if (iters >= IterMAX) TF <- TRUE
  list(ParG=ParG,RSS=RSS,iters=iters,TF=TF,Jg=Jg,H=H)
}   #fit.ellipseLMG.H

JmatrixLMG <- function(XY,A,XYproj){
#Implicit method
n <- dim(XY)[1]
Res  <-  matrix(0,n,1)
s <- sin(A[5])
C <- cos(A[5])
ss <- s*s
cc <- C*C
cs <- C*s
a <- A[3]
b <- A[4]
aa <- a^2
bb <- b^2
ba <- bb-aa
DD <- XY-XYproj
D1 <- cbind(XYproj[,1]-A[1], XYproj[,2]-A[2])
D2 <- cbind(D1[,2], D1[,1])
for (i in 1:n) Res[i]  <- sign(DD[i,] %*% matrix(D1[i,],2)) * norm(matrix(DD[i,]),'F')
du <- D1 %*% rbind(bb*cc+aa*ss, cs*ba)
dv <- D2 %*% rbind(bb*ss+aa*cc, cs*ba)
D3 <- cbind(D1[,1]^2, D1[,2]^2, D1[,1]*D1[,2],matrix(1,n,1))
d3 <-  D3 %*% rbind(ss,cc,-2*cs,-bb) * a
d4 <-  D3 %*% rbind(cc,ss,2*cs,-aa)*b
d5 <-  D3 %*% rbind(-ba*cs,ba*cs,ba*(cc-ss),0)
e  <- sqrt(du^2+dv^2)
J <- cbind(-du/e, -dv/e, d3/e, d4/e, d5/e)
list(Res=Res,J=J)
}

CircleFitByLandau<-function(XY,ParIni=NA,epsilon=0.000001,IterMAX = 50)
{
if (length(ParIni) != 3) ParIni <- NA
if (!is.numeric(ParIni)) ParIni <- NA
if (any(is.na(ParIni))) ParIni <- estimateInitialGuessCircle(XY)
centroid <- apply(XY,2,mean)
X <- XY[,1] - centroid[1]
Y <- XY[,2] - centroid[2]
#  centering the initial guess
ParNew <- matrix(c(ParIni - c(centroid,0)),1)
for (iter in 1:IterMAX)         #  main loop, each run is one iteration
{
ParOld <- ParNew
Dx <- X - ParOld[1]
Dy <- Y - ParOld[2]
D <- sqrt(Dx * Dx + Dy * Dy)
ParNew <- cbind(-mean(Dx/D) , -mean(Dy/D) , 1) * mean(D)
progress <- (norm(ParNew-ParOld,'F'))/(norm(ParOld,'F')+epsilon)
if (progress < epsilon) break #  stopping rule
}   #  the end of the main loop (over iterations)
# decentering the parameter vector for output
Par <- ParOld + c(centroid,0)
Par
}

CircleFitBySpath<-function(XY,ParIni=NA,epsilon=0.000001,IterMAX = 50)
{
if (length(ParIni) != 3) ParIni <- NA
if (!is.numeric(ParIni)) ParIni <- NA
if (any(is.na(ParIni))) ParIni <- estimateInitialGuessCircle(XY)
centroid <- apply(XY,2,mean)
X <- XY[,1] - centroid[1]
Y <- XY[,2] - centroid[2]
#  centering the initial guess
ParNew <- matrix(c(ParIni - c(centroid,0)),1)
for (iter in 1:IterMAX)         #  main loop, each run is one iteration
{
ParOld <- ParNew
Dx <- X - ParOld[1]
Dy <- Y - ParOld[2]
D <- sqrt(Dx * Dx + Dy * Dy)
Mu <- mean(Dx/D)
Mv <- mean(Dy/D)
Mr <- mean(D)
Radius <- (Mu * ParOld[1] + Mv * ParOld[2] + Mr)/(1.0 - Mu * Mu - Mv * Mv)
ParNew <- cbind(-Mu , -Mv , 1) * Radius
progress <- (norm(ParNew-ParOld,'F'))/(norm(ParOld,'F')+epsilon)
if (progress < epsilon) break #  stopping rule
}   #  the end of the main loop (over iterations)
# decentering the parameter vector for output
Par <- ParOld + c(centroid,0)
Par
}

estimateInitialGuessCircle <- function(XY)
{# estimate initial guess for circle LM
x0 <- mean(XY[,1])
y0 <- mean(XY[,2])
c(x0,y0,mean(sqrt((XY[,1]^2+x0^2)+(XY[,2]^2+y0^2))))
}

CurrentIteration<-function(Par,XY){
# computes the objective function F and its derivatives at the current point Par
Dx = matrix(XY[,1] - Par[1],ncol=1)
Dy = matrix(XY[,2] - Par[2],ncol=1)
D = sqrt(Dx^2 + Dy^2)
J = cbind(-Dx / D, -Dy / D,  -1 )
g = matrix(D - Par[3],ncol=1)
F = sqrt(sum(g^2))^2
Radius = mean(D)
list(J=J,g=g,F=F,Radius = Radius)
}

LMcircleFit<-function(XY,ParIni=NA,LambdaIni=1,epsilon=0.000001,IterMAX = 50){
#epsilon = tolerance (small threshold)
#IterMAX = maximal number of (main) iterations; usually 10-20 suffice
if (length(ParIni) != 3) ParIni <- NA
if (!is.numeric(ParIni)) ParIni <- NA
if (is.na(ParIni)) ParIni <- estimateInitialGuessCircle(XY)
lambda_sqrt = sqrt(LambdaIni)#  sqrt(Lambda) is actually used by the code
Par = ParIni# starting with the given initial guess
Jtmp = CurrentIteration(Par,XY)# compute objective function and its derivatives
J <- Jtmp$J
g <- Jtmp$g
F <- Jtmp$F
for (iter in 1:IterMAX){ #  main loop, each run is one (main) iteration
    while (TRUE){ #  secondary loop - adjusting Lambda (no limit on cycles)
        DelPar = mldivide(rbind(J, lambda_sqrt * diag(3)), rbind(g, 0,0,0))   # step candidate
        progress = sqrt(sum(DelPar^2)) / (sqrt(sum(Par^2))+epsilon)
        if (progress < epsilon)  break # stopping rule
        ParTemp = Par - t(DelPar)
        Jtmp = CurrentIteration(ParTemp,XY);  # objective function + derivatives
	JTemp<- Jtmp$J
	gTemp<- Jtmp$g
	FTemp<- Jtmp$F
        if (FTemp < F && ParTemp[3]>0){        #   yes, improvement
           lambda_sqrt = lambda_sqrt/2   # reduce lambda, move to next iteration
           break
        } else {                           #   no improvement
           lambda_sqrt = lambda_sqrt*2 # increase lambda, recompute the step
           next
        }
    }   #   while (1), the end of the secondary loop
    if (progress < epsilon)  break # stopping rule
    Par = ParTemp;  J = JTemp;  g = gTemp;  F = FTemp # update the iteration
}
Par
}

CurrentIterationReduced<-function(Par,XY){
# computes the objective function F and its derivatives at the current point Par
Dx = matrix(XY[,1] - Par[1],ncol=1)
Dy = matrix(XY[,2] - Par[2],ncol=1)
D = sqrt(Dx^2 + Dy^2)
J = cbind(-Dx + mean(Dx), -Dy + mean(Dy) )
Radius = mean(D)
g = matrix(D - Radius,ncol=1)
F = sqrt(sum(g^2))^2
list(J=J,g=g,F=F,Radius = Radius)
}

LMreducedCircleFit<-function(XY,ParIni=NA,LambdaIni=1,epsilon=0.000001,IterMAX = 50){
#epsilon = tolerance (small threshold)
#IterMAX = maximal number of (main) iterations; usually 10-20 suffice
if (length(ParIni) != 3) ParIni <- NA
if (!is.numeric(ParIni)) ParIni <- NA
if (any(is.na(ParIni))) ParIni <- estimateInitialGuessCircle(XY)
lambda_sqrt = sqrt(LambdaIni)#  sqrt(Lambda) is actually used by the code
Par = ParIni[1:2]# starting with the given initial guess
Jtmp = CurrentIterationReduced(Par,XY)# compute objective function and its derivatives
J <- Jtmp$J[,1:2]
g <- Jtmp$g
F <- Jtmp$F
Radius <- Jtmp$Radius
for (iter in 1:IterMAX){ #  main loop, each run is one (main) iteration

    while (TRUE){ #  secondary loop - adjusting Lambda (no limit on cycles)
        DelPar = mldivide(rbind(J, lambda_sqrt * diag(2)), rbind(g, 0,0))   # step candidate
        progress = norm(matrix(DelPar,1),'F')/(Radius+norm(matrix(Par,1),'F')+epsilon)
        if (progress < epsilon)  break # stopping rule
        ParTemp = Par - t(DelPar)
        Jtmp = CurrentIterationReduced(ParTemp,XY)  # objective function + derivatives
	JTemp<- Jtmp$J
	gTemp<- Jtmp$g
	FTemp<- Jtmp$F
	RadiusTemp<- Jtmp$Radius
        if (FTemp < F){        #   yes, improvement
           lambda_sqrt = lambda_sqrt/2   # reduce lambda, move to next iteration
           break
        } else {                           #   no improvement
           lambda_sqrt = lambda_sqrt*2 # increase lambda, recompute the step
           next
        }
    }   #   while (1), the end of the secondary loop
    if (progress < epsilon)  break # stopping rule
    Par = ParTemp;  J = JTemp;  g = gTemp;  F = FTemp;  Radius = RadiusTemp # update the iteration
}
c(Par, Radius) # assembling the full parameter vector "Par" for output
}

calculateCircle<-function(x, y, r, steps=50,sector=c(0,360),randomDist=FALSE, randomFun=runif,noiseFun=NA,...)
{
  points = matrix(0,steps,2)
  if (randomDist) n<-sector[1]+randomFun(steps,...)*(sector[2]-sector[1]) else n<-seq(sector[1],sector[2],length.out=steps)
  if (randomDist) repeat {
  n[which(!(n>=sector[1] & n<=sector[2]))]<-sector[1]+randomFun(sum(!(n>=sector[1] & n<=sector[2])),...)*(sector[2]-sector[1])
  if (all(n>=sector[1] & n<=sector[2])) break
  }
    alpha = n * (pi / 180)
    sinalpha = sin(alpha)
    cosalpha = cos(alpha)
    points[,1]<- x + (r * cosalpha)
    points[,2]<- y + (r * sinalpha)
    if (is.function(noiseFun)) points<-apply(points,1:2,noiseFun)
  points
}

calculateEllipse<-function(x, y, a, b, angle=0, steps=50,sector=c(0,360),randomDist=FALSE, randomFun=runif,noiseFun=NA,...)
{# http://en.wikipedia.org/w/index.php?title=Ellipse&oldid=456212176#Ellipses_in_computer_graphics
  points = matrix(0,steps,2)
  beta = angle * (pi / 180)
  sinbeta = sin(beta)
  cosbeta = cos(beta)
  if (randomDist) n<-sector[1]+randomFun(steps,...)*(sector[2]-sector[1]) else n<-seq(sector[1],sector[2],length.out=steps)
    if (randomDist) repeat {
  n[which(!(n>=sector[1] & n<=sector[2]))]<-sector[1]+randomFun(sum(!(n>=sector[1] & n<=sector[2])),...)*(sector[2]-sector[1])
  if (all(n>=sector[1] & n<=sector[2])) break
  }
    alpha = n * (pi / 180)
    sinalpha = sin(alpha)
    cosalpha = cos(alpha)
    points[,1]<- x + (a * cosalpha * cosbeta - b * sinalpha * sinbeta)
    points[,2]<- y + (a * cosalpha * sinbeta + b * sinalpha * cosbeta)
    if (is.function(noiseFun)) points<-apply(points,1:2,noiseFun)
  points
}

EllipseDirectFit<-function(XY){
# translated to R by Jose Gama 2014
# Original code by: Nikolai Chernov http://www.mathworks.com/matlabcentral/fileexchange/22684-ellipse-fit-direct-method
# A. W. Fitzgibbon, M. Pilu, R. B. Fisher, "Direct Least Squares Fitting of Ellipses", IEEE Trans. PAMI, Vol. 21, pages 476-480 (1999)
# Halir R, Flusser J (1998) Proceedings of the 6th International Conference in Central Europe on Computer Graphics and Visualization, 
# Numerically stable direct least squares fitting of ellipses (WSCG, Plzen, Czech Republic), pp 125–132.
centroid <- apply(XY,2,mean)
D1 <- cbind((XY[,1]-centroid[1])^2, (XY[,1]-centroid[1])*(XY[,2]-centroid[2]), (XY[,2]-centroid[2])^2)
D2 <- cbind(XY[,1]-centroid[1], XY[,2]-centroid[2], matrix(1,dim(XY)[1]))
S1 <- t(D1) %*% D1
S2 <- t(D1) %*% D2
S3 <- t(D2) %*% D2
T <- -solve(S3) %*% t(S2)
M <- S1 + S2 %*% T
M <- rbind(M[3,]/2, -M[2,], M[1,]/2)
evec<-eigen(M)$vectors
cond <- 4*evec[1,]*evec[3,]-evec[2,]^2
A1 <- matrix(evec[,which(cond>0)[1]],3)
A <- rbind(A1, T %*% A1)
A4 <- A[4]-2*A[1]*centroid[1]-A[2]*centroid[2]
A5 <- A[5]-2*A[3]*centroid[2]-A[2]*centroid[1]
A6 <- A[6]+A[1]*centroid[1]^2+A[3]*centroid[2]^2+ A[2]*centroid[1]*centroid[2]-A[4]*centroid[1]-A[5]*centroid[2]
A[4] <- A4;  A[5] <- A5;  A[6] <- A6
A <- A / norm(A,'F')
# # general-form conic equation  ax^2 + bxy + cy^2 +dx + ey + f = 0
# a <- A[1];b <- A[2]/2;C <- A[3];d <- A[4]/2;E <- A[5]/2;f <- A[6]
# x0 <- (C*d-b*E)/(b*b-a*C)
# y0 <- (a*E-b*d)/(b*b-a*C)
# semiaxis.a <- sqrt(2*(a*E^2+C*d^2+f*b^2-2*b*d*E-a*C*f)/((b^2-a*C)*(sqrt((a-C)^2+4*b^2)-(a+C))))
# semiaxis.b <- sqrt(2*(a*E^2+C*d^2+f*b^2-2*b*d*E-a*C*f)/((b^2-a*C)*(-sqrt((a-C)^2+4*b^2)-(a+C))))
# #, x0=x0,y0=y0,semiaxis.a=semiaxis.a,semiaxis.b=semiaxis.b
A
}

EllipseFitByTaubin<-function(XY){
# translated to R by Jose Gama 2014
# Original code by: Nikolai Chernov 
centroid <- apply(XY,2,mean)
Z <- cbind((XY[,1]-centroid[1])^2, (XY[,1]-centroid[1])*(XY[,2]-centroid[2]), (XY[,2]-centroid[2])^2, XY[,1]-centroid[1], XY[,2]-centroid[2], matrix(1,dim(XY)[1]))
M <- t(Z) %*% Z/dim(XY)[1]
P <- rbind(cbind(M[1,1]-M[1,6]^2, M[1,2]-M[1,6]*M[2,6], M[1,3]-M[1,6]*M[3,6], M[1,4], M[1,5]),
cbind(M[1,2]-M[1,6]*M[2,6], M[2,2]-M[2,6]^2, M[2,3]-M[2,6]*M[3,6], M[2,4], M[2,5]),
cbind(M[1,3]-M[1,6]*M[3,6], M[2,3]-M[2,6]*M[3,6], M[3,3]-M[3,6]^2, M[3,4], M[3,5]),
cbind(M[1,4], M[2,4], M[3,4], M[4,4], M[4,5]),
cbind(M[1,5], M[2,5], M[3,5], M[4,5], M[5,5]))
Q <- rbind(cbind(4*M[1,6], 2*M[2,6], 0, 0, 0),
cbind(2*M[2,6], M[1,6]+M[3,6], 2*M[2,6], 0, 0),
cbind(0, 2*M[2,6], 4*M[3,6], 0, 0),
cbind(0, 0, 0, 1, 0),
cbind(0, 0, 0, 0, 1))
V2 <- geigen(P,Q)
V <- V2$vectors
diagD <- V2$values
Dsort <- matrix(sort(diagD),ncol=1)
ID <- matrix(order(diagD),ncol=1)
A <- matrix(V[,ID[1]],ncol=1)
A <- rbind(A, -t(A[1:3]) %*% M[1:3,6])
A4 <- A[4]-2*A[1]*centroid[1]-A[2]*centroid[2]
A5 <- A[5]-2*A[3]*centroid[2]-A[2]*centroid[1]
A6 <- A[6]+A[1]*centroid[1]^2+A[3]*centroid[2]^2+ A[2]*centroid[1]*centroid[2]-A[4]*centroid[1]-A[5]*centroid[2]
A[4] <- A4;  A[5] <- A5;  A[6] <- A6
A <- A / norm(A,'F')
A
}

CircleFitByTaubin<-function(XY){
# translated to R by Jose Gama 2014
# Original code by: Nikolai Chernov 
n <- dim(XY)[1] # number of data points
centroid <- apply(XY,2,mean)
Mxx <- 0; Myy <- 0; Mxy <- 0; Mxz <- 0; Myz <- 0; Mzz <- 0
for (x in 1:n) { 
Xi <- XY[x,1] - centroid[1]
Yi <- XY[x,2] - centroid[2]
Zi <- Xi*Xi + Yi*Yi
Mxy <- Mxy + Xi*Yi
Mxx <- Mxx + Xi*Xi
Myy <- Myy + Yi*Yi
Mxz <- Mxz + Xi*Zi
Myz <- Myz + Yi*Zi
Mzz <- Mzz + Zi*Zi }
Mxx <- Mxx/n
Myy <- Myy/n
Mxy <- Mxy/n
Mxz <- Mxz/n
Myz <- Myz/n
Mzz <- Mzz/n
Mz <- Mxx + Myy
Cov_xy <- Mxx * Myy - Mxy * Mxy
A3 <- 4 * Mz
A2 <- -3 * Mz * Mz - Mzz
A1 <- Mzz * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz - Mz * Mz * Mz
A0 <- Mxz * Mxz * Myy + Myz * Myz * Mxx - Mzz * Cov_xy - 2 * Mxz * Myz * Mxy + Mz * Mz * Cov_xy
A22 <- A2 + A2
A33 <- A3 + A3 + A3
xnew <- 0
ynew <- 1e+20
epsilon <- 1e-12
IterMax <- 20
for (iter in 1:IterMax){
    yold <- ynew
    ynew <- A0 + xnew*(A1 + xnew*(A2 + xnew*A3))
    if (abs(ynew) > abs(yold)) {
       cat('Newton-Taubin goes wrong direction: |ynew| > |yold|\n')
       xnew <- 0
       break
    }
    Dy <- A1 + xnew*(A22 + xnew*A33)
    xold <- xnew
    xnew <- xold - ynew/Dy
    if (abs((xnew-xold)/xnew) < epsilon) break
    if (iter >= IterMax){
        cat('Newton-Taubin will not converge\n')
        xnew <- 0
    }
    if (xnew<0){
        cat(1,'Newton-Taubin negative root:  x=%f\n',xnew,'\n')
        xnew <- 0
    }
}
DET <- xnew*xnew - xnew*Mz + Cov_xy
Center <- cbind(Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy)/DET/2
P <- cbind(Center+centroid , sqrt(Center %*% t(Center)+Mz))
P
}

CircleFitByPratt<-function(XY){
# translated to R by Jose Gama 2014
# Original code by: Nikolai Chernov 
n <- dim(XY)[1] # number of data points
centroid <- apply(XY,2,mean)
Mxx <- 0; Myy <- 0; Mxy <- 0; Mxz <- 0; Myz <- 0; Mzz <- 0
for (x in 1:n) { 
Xi <- XY[x,1] - centroid[1]
Yi <- XY[x,2] - centroid[2]
Zi <- Xi*Xi + Yi*Yi
Mxy <- Mxy + Xi*Yi
Mxx <- Mxx + Xi*Xi
Myy <- Myy + Yi*Yi
Mxz <- Mxz + Xi*Zi
Myz <- Myz + Yi*Zi
Mzz <- Mzz + Zi*Zi }
Mxx <- Mxx/n
Myy <- Myy/n
Mxy <- Mxy/n
Mxz <- Mxz/n
Myz <- Myz/n
Mzz <- Mzz/n
Mz <- Mxx + Myy
Cov_xy <- Mxx * Myy - Mxy * Mxy
Mxz2 <- Mxz * Mxz
Myz2 <- Myz * Myz
A2 <- 4 * Cov_xy - 3 * Mz * Mz - Mzz
A1 <- Mzz * Mz + 4 * Cov_xy * Mz - Mxz2 - Myz2 - Mz * Mz * Mz
A0 <- Mxz2 * Myy + Myz2 * Mxx - Mzz * Cov_xy - 2 * Mxz * Myz * Mxy + Mz * Mz * Cov_xy
A22 <- A2 + A2
epsilon <- 1e-12 
ynew <- 1e+20
IterMax <- 20
xnew <- 0
for (iter in 1:IterMax){
    yold <- ynew
    ynew <- A0 + xnew*(A1 + xnew*(A2 + xnew*xnew*4))
    if (abs(ynew) > abs(yold)) {
       cat('Newton-Pratt goes wrong direction: |ynew| > |yold|\n')
       xnew <- 0
       break
    }
    Dy <- A1 + xnew*(A22 + 16*xnew*xnew)
    xold <- xnew
    xnew <- xold - ynew/Dy
    if (abs((xnew-xold)/xnew) < epsilon) break
    if (iter >= IterMax){
        cat('Newton-Pratt will not converge\n')
        xnew <- 0
    }
    if (xnew<0){
        cat(1,'Newton-Pratt negative root:  x=%f\n',xnew,'\n')
        xnew <- 0
    }
}
DET <- xnew*xnew - xnew*Mz + Cov_xy
Center <- cbind(Mxz*(Myy-xnew)-Myz*Mxy , Myz*(Mxx-xnew)-Mxz*Mxy)/DET/2
P <- cbind(Center+centroid , sqrt(Center %*% t(Center)+Mz+2*xnew))
P
}

CircleFitByKasa<-function(XY){
# translated to R by Jose Gama 2014
# Original code by: Nikolai Chernov 
P <- mldivide(cbind(XY,1), matrix(XY[,1]^2 + XY[,2]^2,ncol=1))
Pout = cbind(P[1]/2 , P[2]/2 , sqrt((P[1]^2+P[2]^2)/4+P[3]))
Pout
}

# http://en.wikipedia.org/wiki/Ellipse#Mathematical_definitions_and_properties
ellipticity <- function(minorAxis,majorAxis) 1 - (minorAxis/majorAxis) # ellipticity = flattening factor
ellipseEccentricity <- function(minorAxis,majorAxis) (sqrt (1 - (minorAxis/majorAxis)^2)) # eccentricity of the ellipse
ellipseFocus <- function(minorAxis,majorAxis) sqrt(minorAxis-majorAxis)^2 # focus of the ellipse
ellipseRa <- function(minorAxis,majorAxis) (1+ellipseEccentricity(minorAxis,majorAxis))*majorAxis # radius at apoapsis (the farthest distance)
ellipseRp <- function(minorAxis,majorAxis) (1-ellipseEccentricity(minorAxis,majorAxis))*majorAxis # radius at periapsis (the closest distance)
ellipse.l <- function(minorAxis,majorAxis)
{# semi-latus rectum l
ra <- ellipseRa(minorAxis,majorAxis)
rp <- ellipseRp(minorAxis,majorAxis)
2*ra*rp/(ra+rp)
}


conic2parametric<-function(A, bv, cv){
# Diagonalise A - find Q, D such at A = Q' * D * Q
# Copyright Richard Brown, this code can be freely used and modified so
# long as this line is retained
# FITELLIPSE : Least squares ellipse fitting demonstration
# Richard Brown, May 28, 2007
# http://www.mathworks.com/matlabcentral/fileexchange/15125-fitellipse-m/content/demo/html/ellipsedemo.html
eTMP<-eigen(A)
D<-diag(eTMP$values)
Q<-eTMP$vectors
Q<-t(Q)
# If the determinant < 0, it's not an ellipse
if (prod(diag(D)) <= 0 )  stop('Linear fit did not produce an ellipse')
# We have b_h' = 2 * t' * A + b'
tV = -0.5 * mldivide(A , bv)
c_h = matrix(tV,1) %*% A %*% tV + t(bv) %*% tV + cv
list(z = tV,a = sqrt(-c_h / D[1,1]),b = sqrt(-c_h / D[2,2]),alpha = atan2(Q[1,2], Q[1,1]))
}

fitbookstein<-function(x){
#FITBOOKSTEIN   Linear ellipse fit using bookstein constraint
#   lambda_1^2 + lambda_2^2 = 1, where lambda_i are the eigenvalues of A
# Copyright Richard Brown, this code can be freely used and modified so
# long as this line is retained
# FITELLIPSE : Least squares ellipse fitting demonstration
# Richard Brown, May 28, 2007
# http://www.mathworks.com/matlabcentral/fileexchange/15125-fitellipse-m/content/demo/html/ellipsedemo.html
# W. Gander, G. H. Golub, R. Strebel, 1994
# Least-Squares Fitting of Circles and Ellipses
# BIT Numerical Mathematics, Springer 
# Convenience variables
m  = dim(x)[1]
x1 = x[, 1]
x2 = x[, 2]
# Define the coefficient matrix B, such that we solve the system
# B *[v; w] = 0, with the constraint norm(w) == 1
B = cbind(x1, x2, rep(1,m), x1^2, sqrt(2) * x1 * x2, x2^2)
# To enforce the constraint, we need to take the QR decomposition
qTMP<-qr(B)
R<-qr.R(qTMP)
Q<-qr.Q(qTMP)
# Decompose R into blocks
R11 = R[1:3, 1:3]
R12 = R[1:3, 4:6]
R22 = R[4:6, 4:6]
# Solve R22 * w = 0 subject to norm(w) == 1
svdTMP<-svd(R22)
U<-svdTMP[["u"]]
S<-diag(svdTMP[["d"]])
V<-svdTMP[["v"]]
w = matrix(V[, 3],3,1)
# Solve for the remaining variables
v = mldivide(-R11 , R12) %*% w
# Fill in the quadratic form
A        = matrix(0,2,2)
A[1]     = w[1]
A[2:3] = 1 / sqrt(2) * w[2]
A[4]     = w[3]
bv       = v[1:2]
c1        = v[3]
# Find the parameters
cTMP<-conic2parametric(A, bv, c1)
list(z=cTMP$z, a=cTMP$a, b=cTMP$b, alpha=cTMP$alpha)
}

fitggk<-function(x){
# Linear least squares with the Euclidean-invariant constraint Trace(A) = 1
# Copyright Richard Brown, this code can be freely used and modified so
# long as this line is retained
# FITELLIPSE : Least squares ellipse fitting demonstration
# Richard Brown, May 28, 2007
# http://www.mathworks.com/matlabcentral/fileexchange/15125-fitellipse-m/content/demo/html/ellipsedemo.html
# W. Gander, G. H. Golub, R. Strebel, 1994
# Least-Squares Fitting of Circles and Ellipses
# BIT Numerical Mathematics, Springer 
# Convenience variables
m  = dim(x)[1]
x1 = x[, 1]
x2 = x[, 2]
# Coefficient matrix
B = cbind(2 * x1 * x2, x2^2 - x1^2, x1, x2, rep(1,m))
v = mldivide(B , -x1^2)
# For clarity, fill in the quadratic form variables
A        = matrix(0,2,2)
A[1]     = 1-v[2]
A[2:3] = v[1]
A[2,2]     = v[2]
bv       = v[3:4]
c1        = v[5]
# Find the parameters
cTMP<-conic2parametric(A, bv, c1)
list(z=cTMP$z, a=cTMP$a, b=cTMP$b, alpha=cTMP$alpha)
}






