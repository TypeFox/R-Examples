DoMohrFig1 <-
function(Stensor=matrix(c(5,1, 1, 3), ncol=2), rot1=NULL)
{

  if(missing(Stensor)) {Stensor=matrix(c(5,1, 1, 3), ncol=2)  }
  if(missing(rot1))
    {
      s1 = Stensor 
    }
  else
    {
      
      s1 = Stensor %*% rot1
      
    }
  
  ES = eigen(s1)
  
  Rmohr = sqrt(  ((s1[1,1]-s1[2,2])^2)/4 +  s1[1,2]^2 )

Save = (ES$values[1]+ES$values[2])/2

ps1 = ES$values[1]
ps2 = ES$values[2]

ex = Save
why = 0

cmohr = GEOmap::darc( rad=Rmohr, ang1=0, ang2=360, x1=ex, y1=why, n=1)

RNGM = range( cmohr$x)
Prange = c(0 , RNGM[2]+0.05*diff(RNGM))

plot(range( Prange ) , range(cmohr$y), type='n', asp=1, axes=FALSE, ann=FALSE)
abline(v=0, h=0, lty=2)

lines(cmohr, lwd=2)
points(ex, why)

u = par("usr")

segments(ex, 0, ex,   0.9*u[3])

segments(ps1, 0.95*u[4] , ps1,    0)
arrows( 0, 0.90*u[4] ,ps1, 0.90*u[4],    length=0.1, code=2) 
text(mean(c( 0, ps1 ) ), 0.90*u[4], labels=expression(sigma[1]), pos=3, xpd=TRUE)

segments(ps2, 0, ps2,    0.6*u[3])
arrows( 0, 0.58*u[3] ,ps2, 0.58*u[3],    length=0.1, code=2) 
text(mean(c( 0, ps2 ) ),  0.58*u[3],  labels=expression(sigma[2]), pos=3, xpd=TRUE)

segments(ex, 0, ex,   0.9*u[3])
arrows( 0, 0.88*u[3] ,ex, 0.88*u[3],    length=0.1, code=2) 
text(mean(c( 0, ex ) ),  0.88*u[3],  labels=expression(sigma[ave]), pos=3, xpd=TRUE)

points(s1[1,1], -s1[2,1])
points(s1[2,2], s1[2,1])

segments(s1[1,1], -s1[2,1], s1[2,2], s1[2,1])

text(s1[1,1], -s1[2,1], labels="X", adj =c(-1,1), font=2, xpd=TRUE)
text(s1[2,2], s1[2,1], labels="Y", adj =c(1.2,-.7), font=2, xpd=TRUE)

points(ps1, 0)
points(ps2, 0)
text(ps1, 0, labels="A",  xpd=TRUE, adj =c(-.5,0) , font=2)

text(ps2, 0, labels="B",  xpd=TRUE, adj =c(1.5, 0) , font=2)


segments(s1[1,1], -s1[2,1], s1[1,1], u[3])
arrows( 0, 0.98*u[3] , s1[1,1], 0.98*u[3],    length=0.1, code=3) 
text(mean(c( 0, s1[1,1] ) ), 0.98*u[3] ,  labels=expression(sigma[x]), pos=3, xpd=TRUE)


segments(s1[2,2], s1[2,1], s1[2,2], 0.8*u[4])
arrows( 0, 0.78*u[4] , s1[2,2], 0.78*u[4],    length=0.1, code=3) 

sigylabx = mean(c( 0, s1[2,2] ) )
text(sigylabx, 0.78*u[4] ,  labels=expression(sigma[y]), pos=1, xpd=TRUE)

segments(s1[1,1], -s1[2,1], u[2] , -s1[2,1])
arrows( (u[2]+ps1)/2 , 0,  (u[2]+ps1)/2,-s1[2,1], code=3, length=0.1)
text( (u[2]+ps1)/2, mean(c(0, -s1[2,1])),  xpd=TRUE, labels=expression(tau[xy]), pos=4, xpd=TRUE)
 
g1x = mean(c( 0, sigylabx))

segments(s1[2,2], s1[2,1], g1x,  s1[2,1])

arrows(  mean(c( g1x, sigylabx)), 0,  mean(c( g1x, sigylabx)) , s1[2,1], code=3, length=0.1)
text(  mean(c( g1x, sigylabx)), mean(c(0, s1[2,1])),  labels=expression(tau[xy]), pos=2, xpd=TRUE)

}

