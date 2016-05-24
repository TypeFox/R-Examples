`imageSH` <-
function( phiS, del, lam, SCALE=FALSE, UP=FALSE, col=NULL )
{
if(missing(SCALE)) { SCALE=FALSE }
if(missing(UP)) { UP=FALSE }
if(missing(col)) { col=heat.colors(20) }
 DEG2RAD = pi/180
    RAD2DEG = 180/pi

x = seq(-1, 1, 0.01)
y = x

X = matrix(rep(x, length(y)), nrow= length(x))
Y = t(X)

p = RAD2DEG*(pi/2 -atan2(Y, X))
p[p<0] = p[p<0] + 360

R = sqrt(X^2+Y^2)
R[R>1] = NaN
dip =RAD2DEG*2*asin(R/sqrt(2))

if(UP==TRUE) { dip = 180-dip }
G = radSH( del, phiS, lam, dip, p)
image(x,y,G, col = col,asp=1,  xlab='', ylab='', axes=FALSE )
if(SCALE==TRUE) { imageSCALE( G, col = col, x=1.1, labels="breaks", nlab=10) }

}

