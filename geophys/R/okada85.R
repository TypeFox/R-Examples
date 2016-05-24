okada85 <-
function(e=0,n=0,depth=0,strike=20, dip=20, L=5, W=3, rake=20, slip=1, U3=1, nu= 0.25){


  ############   documentation from the matlab function
#          [E,N] = meshgrid(linspace(-10,10,50));
#          [uE,uN,uZ] = okada85(E,N,2,30,70,5,3,-45,1,1,'plot');
#          figure, surf(E,N,uN)
#


#       [uE,uN,uZ,uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...
#          E,N,DEPTH,STRIKE,DIP,LENGTH,WIDTH,RAKE,SLIP,OPEN)
#       computes displacements, tilts and strains at the surface of an elastic
#       half-space, due to a dislocation defined by RAKE, SLIP, and OPEN on a
#       rectangular fault defined by orientation STRIKE and DIP, and size LENGTH and
#       WIDTH. The fault centroid is located (0,0,-DEPTH).
#
#          E,N    : coordinates of observation points in a geographic referential
#                   (East,North,Up) relative to fault centroid (units are described below)
#          DEPTH  : depth of the fault centroid (DEPTH > 0)
#          STRIKE : fault trace direction (0 to 360 relative to North), defined so
#                   that the fault dips to the right side of the trace
#          DIP    : angle between the fault and a horizontal plane (0 to 90)
#          LENGTH : fault length in the STRIKE direction (LENGTH > 0)
#          WIDTH  : fault width in the DIP direction (WIDTH > 0)
#          RAKE   : direction the hanging wall moves during rupture, measured relative
#                   to the fault STRIKE (-180 to 180).
#          SLIP   : dislocation in RAKE direction (length unit)
#          OPEN   : dislocation in tensile component (same unit as SLIP)
#
#       returns the following variables (same matrix size as E and N):
#          uN,uE,uZ        : displacements (unit of SLIP and OPEN)
#          uZE,uZN         : tilts (in rad * FACTOR)
#          uNN,uNE,uEN,uEE : horizontal strains POSITIVE = COMPRESSION (unit of FACTOR)
#
#       Length unit consistency: E, N, DEPTH, LENGTH, and WIDTH must have the same
#       unit (e.g. km) which can be different from that of SLIP and OPEN (e.g. m) but
#       with a possible FACTOR on tilt and strain results (in this case, an
#       amplification of km/m = 1000). To have FACTOR = 1 (tilt in radians and
#       correct strain unit), use the same length unit for all aforesaid variables.
#
#       [...] = OKADA85(...,NU) specifies Poisson's ratio NU (default is 0.25 for
#       an isotropic medium).
#
#       Formulas and notations from Okada [1985] solution excepted for strain
#       convention (here positive strain means compression), and for the fault
#       parameters after Aki & Richards [1980], e.g.:
#             DIP=90, RAKE=0   : left lateral (senestral) strike slip
#             DIP=90, RAKE=180 : right lateral (dextral) strike slip
#             DIP=70, RAKE=90  : reverse fault
#             DIP=70, RAKE=-90 : normal fault
#
#       It is also possible to produce partial outputs, with following syntax:
#          [uE,uN,uZ] = OKADA85(...) for displacements only;
#          [uE,uN,uZ,uZE,uZN] = OKADA85(...) for displacements and tilts;
#          [uE,uN,uZ,uNN,uNE,uEN,uEE] = OKADA85(...) for displacements and strains;
#          [uZE,uZN] = OKADA85(...) for tilts only;
#          [uZE,uZN,uNN,uNE,uEN,uEE] = OKADA85(...) for tilts and strains;
#          [uNN,uNE,uEN,uEE] = OKADA85(...) for strains only.
#
#       Note that vertical strain components can be obtained with following equations:
#          uNZ = -uZN;
#          uEZ = -uZE;
#          uZZ = -(uEE + uNN)*NU/(1-NU);
#
#       [...] = OKADA85(...,'plot') or OKADA85(...) without output argument produces
#       a 3-D figure with fault geometry and dislocation at scale (if all of the fault
#       parameters are scalar).
#
#       Equations are all vectorized excepted for argument DIP which must be
#       a scalar; all other arguments can be scalar or matrix of the same size.
#
#       Example:
#
#          [E,N] = meshgrid(linspace(-10,10,50));
#          [uE,uN,uZ] = okada85(E,N,2,30,70,5,3,-45,1,1,'plot');
#          figure, surf(E,N,uN)
#
#       considers a 5x3 fault at depth 2, N30-strike, 70-dip, and unit dislocation
#       in all directions (reverse, senestral and open). Displacements are computed
#       on a regular grid from -10 to 10, and North displacements are plotted as a
#       surface.
#
#
#       Author: Francois Beauducel <beauducel@ipgp.fr>
#          Institut de Physique du Globe de Paris
#       Created: 1997
#       Updated: 2011-03-08
#
#       References:
#          Aki K., and P. G. Richards, Quantitative seismology, Freemann & Co,
#             New York, 1980.
#          Okada Y., Surface deformation due to shear and tensile faults in a
#             half-space, Bull. Seismol. Soc. Am., 75:4, 1135-1154, 1985.
#
#       Acknowledgments: Dmitry Nicolsky, University of Alaska

#       Development history:
#          [2011-03-08]: help review.
#          [2011-03-06]: new optional argument to plot fault geometry with
#             output arguments, and bug correction for the fault centroid position
#             (in calculation and plot).
#          [2010-11-29]: change coordinates and depth to fault centroid
#             (instead of middle top edge).
#          [2010-09-24]: bugs correction in the syntax of I1, K2 and uyy_tf
#             functions, affecting some components. Detected by Dmitry Nicolsky.
#

  ########   translate to R by Lan
  ########  modified by J. M. Lees  Thu Jun 16 10:03:00 EDT 2011

if(missing(nu))  { nu= 0.25 }
########

uxxds = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return(xi*q/R^3+ J3(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta))
}
########

uxxss = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return(xi^2*q*A(eta,R)- J1(xi,eta,q,delta,nu,R)*sin(delta))
}

########

uxxtf = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return(xi*q^2*A(eta,R)+ J3(xi,eta,q,delta,nu,R)*sin(delta)^2)
}


########

uxyds = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
yb = eta*cos(delta) + q*sin(delta)
return(yb*q/R^3 - sin(delta)/R+ J1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta))
}

########

uxyss = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
return( xi^3*db/(R^3*(eta^2 + q^2))- (xi^3*A(eta,R) + J2(xi,eta,q,delta,nu,R))*sin(delta))
}

########

uxytf = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
return(-db*q/R^3- xi^2*q*A(eta,R)*sin(delta)+ J1(xi,eta,q,delta,nu,R)*sin(delta)^2)
}
########

uyds= function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return((eta*cos(delta) + q*sin(delta))*q/(R*(R + xi))	+ cos(delta)*atan(xi*eta/(q*R))- Ione(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta))
}

########

uyss = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return((eta*cos(delta) + q*sin(delta))*q/(R*(R + eta))+ q*cos(delta)/(R + eta)+ Itwo(eta,q,delta,nu,R)*sin(delta))
}

uytf = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
u = -(eta*sin(delta) - q*cos(delta))*q/(R*(R + xi))- sin(delta)*(xi*q/(R*(R + eta))- atan(xi*eta/(q*R)))- Ione(xi,eta,q,delta,nu,R)*sin(delta)^2
}
########

uyxds= function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
yb = eta*cos(delta) + q*sin(delta)
return(yb*q/R^3+ q*cos(delta)/(R*(R + eta))+ J1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta))
}

########

uyxss = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return(xi*q/R^3*cos(delta)+ (xi*q^2*A(eta,R) - J2(xi,eta,q,delta,nu,R))*sin(delta))
}


uyxtf=function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return(q^2/R^3*cos(delta)+ q^3*A(eta,R)*sin(delta)+ J1(xi,eta,q,delta,nu,R)*sin(delta)^2)
}

########

uyyds=function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
yb = eta*cos(delta) + q*sin(delta)
return(yb^2*q*A(xi,R)- (2*yb/(R*(R + xi)) + xi*cos(delta)/(R*(R + eta)))*sin(delta)	+ J2(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta))
}

########


uyyss = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
yb = eta*cos(delta) + q*sin(delta)
return(yb*q/R^3*cos(delta)+ (q^3*A(eta,R)*sin(delta) - 2*q*sin(delta)/(R*(R + eta))	- (xi^2 + eta^2)/R^3*cos(delta) - J4(xi,eta,q,delta,nu,R))*sin(delta))
}


uyytf=function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
yb = eta*cos(delta) + q*sin(delta)
return((yb*cos(delta) - db*sin(delta))*q^2*A(xi,R)- q*sin(2*delta)/(R*(R + xi))- (xi*A(eta,R) - J2(xi,eta,q,delta,nu,R))*sin(delta)^2)
}

########

uzds = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
return( db*q/(R*(R + xi))+ sin(delta)*atan(xi*eta/(q*R))- Ifive(xi,eta,q,delta,nu,R,db)*sin(delta)*cos(delta))
}


########

uzss = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
return((eta*sin(delta) - q*cos(delta))*q/(R*(R + eta))+ q*sin(delta)/(R + eta) + Ifour(db,eta,q,delta,nu,R)*sin(delta))
}


uztf = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
return((eta*cos(delta) + q*sin(delta))*q/(R*(R + xi)) + cos(delta)*(xi*q/(R*(R + eta))- atan(xi*eta/(q*R)))- Ifive(xi,eta,q,delta,nu,R,db)*sin(delta)^2)
}
########

uzxds = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
return(db*q/R^3+ q*sin(delta)/(R*(R + eta))+ K3(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta))
}


########

uzxss = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return(-xi*q^2*A(eta,R)*cos(delta)+ ((xi*q)/R^3 - K1(xi,eta,q,delta,nu,R))*sin(delta))
}

########

uzxtf = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
return(q^2/R^3*sin(delta)- q^3*A(eta,R)*cos(delta)+ K3(xi,eta,q,delta,nu,R)*sin(delta)^2)
}

########


uzyds = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
yb = eta*cos(delta) + q*sin(delta)
return(yb*db*q*A(xi,R)- (2*db/(R*(R + xi)) + xi*sin(delta)/(R*(R + eta)))*sin(delta)+ K1(xi,eta,q,delta,nu,R)*sin(delta)*cos(delta))
}

########

uzyss = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
yb = eta*cos(delta) + q*sin(delta)
return(db*q/R^3*cos(delta)+ (xi^2*q*A(eta,R)*cos(delta) - sin(delta)/R + yb*q/R^3 - K2(xi,eta,q,delta,nu,R))*sin(delta))
}

########

uzytf = function(xi,eta,q,delta,nu) {
R = sqrt(xi^2 + eta^2 + q^2)
db = eta*sin(delta) - q*cos(delta)
yb = eta*cos(delta) + q*sin(delta)
return((yb*sin(delta) + db*cos(delta))*q^2*A(xi,R)+ xi*q^2*A(eta,R)*sin(delta)*cos(delta)- (2*q/(R*(R + xi)) - K1(xi,eta,q,delta,nu,R))*sin(delta)^2)
}

########


########
uxtf = function(xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
zed = q^2 /(R*(R + eta)) - Ithree(eta,q,delta,nu,R)*sin(delta)^2 
return(zed )
}
########

uxss = function (xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
##  a1=xi*q/(R*(R + eta))
## a2 = atan(xi*eta/(q*R))
## a3 = Ione(xi,eta,q,delta,nu,R)*sin(delta)

##  cat(paste(format(a1), format(a2), format(a3), sep="  "), sep="\n" )


zed=xi*q/(R*(R + eta)) + atan(xi*eta/(q*R)) + Ione(xi,eta,q,delta,nu,R)*sin(delta)


## bed = a1 +a2 + a3

## print(paste("bed zed ", bed, zed))

return(zed)
}

########
uxds= function (xi,eta,q,delta,nu){
R = sqrt(xi^2 + eta^2 + q^2)
zed = q/R - Ithree(eta,q,delta,nu,R) * sin(delta) *cos(delta)
return(zed)
}

########
K3 = function(xi,eta,q,delta,nu,R){
db = eta*sin(delta) - q*cos(delta)
yb = eta*cos(delta) + q*sin(delta)
if (cos(delta) > 2.2204e-16)
	return((1 - 2*nu) * 1/cos(delta) * (q/(R*(R + eta)) - yb/(R*(R + db)))) else
	return((1 - 2*nu) * sin(delta)/(R + db) * (xi^2/(R*(R + db)) - 1))
}
########
K2 = function (xi,eta,q,delta,nu,R){
return((1 - 2*nu) * (-sin(delta)/R + q*cos(delta)/R*(R + eta)) - K3(xi,eta,q,delta,nu,R))
}
########
K1 = function(xi,eta,q,delta,nu,R){
db = eta*sin(delta) - q*cos(delta)
if (cos(delta) > 2.2204e-16)
	return((1 - 2*nu) * xi/cos(delta) * (1/(R*(R + db)) - sin(delta)/(R*(R + eta)))) else
	return((1 - 2*nu) * xi*q/(R + db)^2)
}
########
J4=function(xi,eta,q,delta,nu,R){
return((1 - 2*nu) * (-cos(delta)/R - q*sin(delta)/(R*(R + eta)))- J1(xi,eta,q,delta,nu,R))
}
########
J3 = function(xi,eta,q,delta,nu,R){
return((1 - 2*nu) * -xi/(R*(R + eta))- J2(xi,eta,q,delta,nu,R))
}

########
J2=function(xi,eta,q,delta,nu,R){
db = eta*sin(delta) - q*cos(delta)
yb = eta*cos(delta) + q*sin(delta)
if (cos(delta) > 2.2204e-16)
	return((1 - 2*nu) * 1/cos(delta) * xi*yb/(R*(R + db)^2)- sin(delta)/cos(delta)*K1(xi,eta,q,delta,nu,R)) else
	return((1 - 2*nu)/2 * xi*sin(delta)/(R + db)^2 * (2*q^2/(R*(R + db)) - 1))
}
########
J1 = function(xi,eta,q,delta,nu,R){
db = eta*sin(delta) - q*cos(delta)
if (cos(delta) > 2.2204e-16)
	return((1 - 2*nu) * 1/cos(delta) * (xi^2/(R*(R + db)^2) - 1/(R + db))- sin(delta)/cos(delta)*K3(xi,eta,q,delta,nu,R)) else
	return((1 - 2*nu)/2 * q/(R + db)^2 * (2*xi^2/(R*(R + db)) - 1))
}
########


Itwo = function(eta,q,delta,nu,R){
return((1 - 2*nu) * (-log(R + eta)) - Ithree(eta,q,delta,nu,R))
}

########
Ione <- function(xi,eta,q,delta,nu,R){
db = eta*sin(delta) - q*cos(delta)

if (cos(delta) > 2.2204e-16)
  {
zed  = (1 - 2*nu) * (-xi/(cos(delta)*(R+db))) - sin(delta)/cos(delta) * Ifive(xi,eta,q,delta,nu,R,db)
}
else
  {
zed = ( -(1 - 2*nu)/2 * xi*q/(R + db)^2)
}


return(zed)

}
########
Ithree = function(eta,q,delta,nu,R){
yb = eta*cos(delta) + q*sin(delta)
db = eta*sin(delta) - q*cos(delta)
if (cos(delta) > 2.2204e-16)
  {
    zed = ((1 - 2*nu) * (yb/(cos(delta)*(R + db)) - log(R + eta))
           + sin(delta)/cos(delta) * Ifour(db,eta,q,delta,nu,R))
      }
else
  {
	zed = ((1 - 2*nu)/2 * (eta/(R + db) + yb*q/(R + db)^2 - log(R + eta)))
      }

return(zed)
}
########


Ifour = function(db,eta,q,delta,nu,R){
  if (cos(delta) > 2.2204e-16)
    {
      zed = (1 - 2*nu) * 1/cos(delta) * (log(R + db) - sin(delta)*log(R + eta))
    }
  else
    {
      zed =  -(1 - 2*nu) * q/(R + db)    
    }

  return(zed)
}
########
Ifive = function(xi, eta, q, delta, nu, R, db){
X=sqrt(xi^2+q^2)
if (cos(delta) > 2.2204e-16)
  {
zed = ((1-2*nu)*2/cos(delta)
        *atan((eta*(X+q*cos(delta))+X*(R+X)*sin(delta))/(xi*(R+X)*cos(delta))))
}
else
  {
zed = (-(1-2*nu)*xi*sin(delta)/(R+db))

}

return(zed)

}
########


A = function(x,R){
return((2*R + x)/(R^3*(R + x)^2))
}

chinnery<- function(f,x,p,L,W,q,delta,nu)
{
  fun = match.fun(f)
y = (fun(x,p,q,delta,nu)
    -fun(x,p-W,q,delta,nu)
    -fun(x-L, p, q, delta, nu)
    +fun(x-L, p-W,q,delta,nu));

  

return(y)
}













strike = strike*pi/180
delta = dip*pi/180         ###   dip in matlab code, delta in okada
rake = rake*pi/180

##    Defines dislocation in the fault plane system
U1 = cos(rake)*slip
U2 = sin(rake)*slip


# Converts fault coordinates (E,N,DEPTH) relative to centroid
#  into Okada's reference system (X,Y,D)
# depth = depth+sin(delta)*W
dtop = depth+sin(delta)*W/2

ec = e + cos(strike)*cos(delta)*W/2;
nc = n - sin(strike)*cos(delta)*W/2;

x = cos(strike)*nc +sin(strike)*ec + L/2
y = sin(strike)*nc -cos(strike)*ec +cos(delta)*W


##  Variable substitution (independent from xi and eta)
p = y*cos(delta)+dtop*sin(delta)
q = y*sin(delta) -dtop*cos(delta)


#################   Displacements
###print(c(U1, U2, U3))

###print(c( x,p,L,W,q,delta,nu))

c1 = chinnery(uxss,x,p,L,W,q,delta,nu)

c2 =  chinnery(uxds,x,p,L,W,q,delta,nu)

c3 =  chinnery(uxtf,x,p,L,W,q,delta,nu)

###print("c 1,2,3:")
###print(c(c1, c2, c3))

ux = -U1/(2*pi) *c1 -U2/(2*pi) *c2+U3/(2*pi) * c3



c1 = chinnery(uyss,x,p,L,W,q,delta,nu)

c2 =  chinnery(uyds,x,p,L,W,q,delta,nu)

c3 =  chinnery(uytf,x,p,L,W,q,delta,nu)

uy = -U1/(2*pi) * c1    -U2/(2*pi) * c2  +U3/(2*pi) * c3



c1 =chinnery(uzss,x,p,L,W,q,delta,nu)
c2 =chinnery(uzds,x,p,L,W,q,delta,nu)
c3 =chinnery(uztf,x,p,L,W,q,delta,nu)

uz = -U1/(2*pi) * c1-U2/(2*pi) * c2 +U3/(2*pi) * c3
	
ue = sin(strike)*ux - cos(strike)*uy
un = cos(strike)*ux + sin(strike)*uy

###print("ux, uy, uz, ue, un")
###print(c(ux, uy, uz, ue, un))

c1=chinnery(uzxss,x,p,L,W,q,delta,nu)
c2=chinnery(uzxds,x,p,L,W,q,delta,nu)
c3=chinnery(uzxtf,x,p,L,W,q,delta,nu)


uzx = -U1/(2*pi) * c1  -U2/(2*pi) * c2 +U3/(2*pi)  * c3


c1=chinnery(uzyss,x,p,L,W,q,delta,nu)
c2=chinnery(uzyds,x,p,L,W,q,delta,nu)
  c3=chinnery(uzytf,x,p,L,W,q,delta,nu)

uzy = -U1/(2*pi) * c1 -U2/(2*pi) * c3 + U3/(2*pi)    * c3
	
uze = sin(strike)*uzx - cos(strike)*uzy
uzn = cos(strike)*uzx + sin(strike)*uzy

c1 =chinnery(uxxss,x,p,L,W,q,delta,nu)
c2 =chinnery(uxxds,x,p,L,W,q,delta,nu)
c3 =chinnery(uxxtf,x,p,L,W,q,delta,nu)
    

uxx = -U1/(2*pi) * c1  -U2/(2*pi) * c2 +U3/(2*pi) * c3

c1 =chinnery(uxyss,x,p,L,W,q,delta,nu)
c2 =chinnery(uxyds,x,p,L,W,q,delta,nu)
c3 =chinnery(uxytf,x,p,L,W,q,delta,nu)
    

uxy = -U1/(2*pi) * c1 -U2/(2*pi) * c2 +U3/(2*pi) * c3

c1 =chinnery(uyxss,x,p,L,W,q,delta,nu)
c2 =chinnery(uyxds,x,p,L,W,q,delta,nu)
c3 =chinnery(uyxtf,x,p,L,W,q,delta,nu)
    
uyx = -U1/(2*pi) * c1  -U2/(2*pi) * c2 +U3/(2*pi) * c3


c1 =chinnery(uyyss,x,p,L,W,q,delta,nu)
c2 =chinnery(uyyds,x,p,L,W,q,delta,nu)
c3 =chinnery(uyytf,x,p,L,W,q,delta,nu)
    
uyy = -U1/(2*pi) * c1 -U2/(2*pi) * c2 +U3/(2*pi) * c3
	
unn = cos(strike)^2*uxx + sin(2*strike)*(uxy + uyx)/2 + sin(strike)^2*uyy

une = sin(2*strike)*(uxx - uyy)/2 + sin(strike)^2*uyx - cos(strike)^2*uxy

uen = sin(2*strike)*(uxx - uyy)/2 - cos(strike)^2*uyx + sin(strike)^2*uxy

uee = sin(strike)^2*uxx - sin(2*strike)*(uyx + uxy)/2 + cos(strike)^2*uyy


return(list(uE=ue,uN=un,uZ=uz,uZE=uze,uZN=uzn,uNN=unn,uNE=une,uEN=uen,uEE=uee))
}

