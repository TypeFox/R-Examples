require(geophys)

opar=par(no.readonly=TRUE)

s1 = matrix(c(5,2,2,3) , byrow=TRUE, ncol=2, nrow=2)
Stensor = s1

theta = 0*pi/180

rot1  = cbind( c( cos(theta), -sin(theta)), c(sin(theta), cos(theta)))

DoMohrFig1(s1, rot1) 


############  routine visualization



s1 = matrix(c(15,-2,-2,3) , byrow=TRUE, ncol=2, nrow=2)
Stensor = s1

DoMohr(Stensor) 




############ 3D  visualization


Stensor = matrix(c(
15, 0, 0,
0, 10, 0,
0,  0, 5), ncol=3)

M = DoMohr(Stensor)


par(ask=FALSE)

Stensor = matrix(c(
15, 0, 0,
0, 10, 0,
0,  0, 5), ncol=3)

stress(Stensor=Stensor)

############  start with a predefined plane
P1 = c(0.2, 1, 1, 0)
P2 = c(1, 0.1, 1, 0)
P3 = c(1, 1, 0.4, 0)

  S = stressSETUP(P1, P2, P3, xscale=30   )

stress(PPs = S$PPs, Rview =S$Rview, xscale = S$xscale, Stensor=Stensor )


  S = stressSETUP( )

Nvec = NORMvec(S$PPs, S$xscale, S$Rview, S$aglyph , add = FALSE)

Stensor = matrix(c(
15, 0, 0,
0, 8, 0,
0,  0, 5), ncol=3)
Mstress  = Maxstress(Nvec, Stensor)

DoMohr(Stensor)
 axis(1)
axis(2)
points(Mstress$sigNORMmax , Mstress$tauSHEARmax, pch=21, col='blue'  , bg='gold' )

u=par('usr')

segments(0, Mstress$tauSHEARmax, Mstress$sigNORMmax ,
Mstress$tauSHEARmax, lty=2, col='green' , lwd=3 )

text(mean(c(0, Mstress$tauSHEARmax)),  Mstress$tauSHEARmax,
"MaxShear in Plane", pos=3)


segments(Mstress$sigNORMmax , u[3] , Mstress$sigNORMmax ,
Mstress$tauSHEARmax, lty=2, col='purple' , lwd=3  )

text(Mstress$sigNORMmax , u[3], "MaxNormal stress", adj=c(0,-1) )


GG = randpoles(30, 40, 10, opt="norm", N=20)

graphics.off()

par(ask=TRUE)

dev.new(width=10, height=6 )

par(mfrow=c(1,2))
    
RFOC::net()
RFOC::qpoint(30, 40, col = "red", UP=FALSE )
RFOC::qpoint(GG$az, GG$dip, col = "blue", UP=FALSE )
rplane = RFOC::lowplane(30-90, 40, UP=TRUE, col='red')


for(i in 1:length(GG$az))
  {
    rplane = RFOC::lowplane(GG$az[i]-90, GG$dip[i], UP=TRUE, col='blue')
  }

g = RSEIS::TOCART(GG$az, GG$dip)

AA = DoMohr(Stensor)

for(i in 1:length(g$az))
  {
    KVEC = c(g$x[i], g$y[i],g$z[i])
    Mstress  = Maxstress(KVEC, Stensor)
    points(Mstress$sigNORMmax , Mstress$tauSHEARmax, pch=21, col='blue'  , bg='gold' )

  }


par(opar)

