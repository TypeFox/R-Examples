
Ellipsoidal.Distance<-function(olat, olon, tlat, tlon, a=6378137.0 , b=6356752.314, tol=10^(-12))
{

######   Vincenty's formulae are two related iterative methods used in
  #######             geodesy to calculate the distance
######    between two points on the surface of an spheroid, developed by
  ######  Thaddeus Vincenty in 1975. They are based on the assumption that
  ##########              the figure of the Earth is an oblate
######    spheroid, and hence are more accurate than methods such as
  ###########     great-circle distance which assume a spherical Earth.

######   The first (direct) method computes the location of a point which is a given
  ######  distance and azimuth (direction) from another point. The second (inverse) method
  ######  computes the geographical distance and azimuth between two given points.
  ######  They have been widely used in geodesy because they are
  ######  accurate to within 0.5 mm (0.020) on the Earth ellipsoid

  #####  default is WGS-84, these are in meters
  
if(missing(a)) a = 6378137.0 
if(missing(b)) b = 6356752.314
if(missing(tol))  tol = 10^(-12)

##   daz = distaz(12, 23, 32,    65)
##   ed = Ellipsoidal.Distance(12, 23, 32,    65)

###  olat= 12; olon=23; tlat=-32;    tlon=-65

###    R.MAPK = 6378.2064;
#   ed =  Ellipsoidal.Distance(12, 23, 32,    65, a=R.MAPK*1000, b=R.MAPK*1000)

  ######  http://en.wikipedia.org/wiki/Vincenty%27s_formulae

err=0
GIVE =  list(dist=NA, az=NA, revaz=NA, err=err)

if(is.na(olat) | is.na(tlat) | is.na(tlon) | is.na(olon))
  {
    return(GIVE)
  }

if(olat < -90 | olat > 90){  return(GIVE)  }
if(tlat < -90 | tlat > 90){  return(GIVE)  }

if( olat == tlat & olon == tlon )
  {
    GIVE =  list(dist=0, az=NA, revaz=NA, err=err)
    return(GIVE)
  }


f = (a-b)/a

##############  latitudes
phi1 = olat*pi/180
phi2 = tlat*pi/180

U1 = atan((1-f)*tan(phi1))
U2 = atan((1-f)*tan(phi2))

cU2 = cos(U2)
cU1 = cos(U1)

sU2 = sin(U2)
sU1 = sin(U1)

###########   longitudes 
lam1 = olon*pi/180
lam2 = tlon*pi/180

L = lam2-lam1

lam = L

####   tolerance
K = tol+1

while(K>tol)
  {

slam = sin(lam)
clam = cos(lam)
    
sinsig = sqrt((cU2*slam)^2 + (cU1*sU2 - sU1*cU2*clam)^2)
if(sinsig==0) {
  print("1 aborting Ellipsoidal.Distance")
  return(GIVE) }
cossig = sU1*sU2 + cU1*cU2*clam

sig = atan2(sinsig, cossig)
sinalpha = (cU1*cU2*slam)/sinsig

cossqalpha = (1-sinalpha^2)

if(cossqalpha==0) {
   print("2 aborting Ellipsoidal.Distance")
  return(GIVE) }

cos2sigm = cossig - (2*sU1*sU2)/cossqalpha

C = (f/16)*cossqalpha*(4+f*(4-3*cossqalpha))

lam2 = L+(1-C)*f*sinalpha*(sig+C*sinsig*(cos2sigm+C*cossig*(-1+2*cos2sigm^2)))

K = abs(lam2-lam)

lam=lam2

  }

usq = cossqalpha * (a^2 - b^2)/b^2

A = 1 + (usq/16384)*(4096+usq*(-768+usq*(320-175*usq)))

B = usq*(256 +usq*(-128+usq*(74-47*usq)))/1024

delsig = B*sinsig*(cos2sigm+0.25*B*(cossig*(-1+2*cos2sigm^2)-(1/6)*B*cos2sigm*(-3+4*sinsig^2)*(-3+4*cos2sigm^2)))

s = b*A*(sig-delsig)

alpha1 = atan2(cU2*slam, cU1*sU2 - sU1*cU2*clam)
alpha2 = atan2(cU1*slam, (-sU1*cU2 + cU1*sU2*clam))

err=1
GIVE =  list(dist=s/1000, az=alpha1*180/pi, revaz=alpha2*180/pi, err=err)


return(GIVE)

}









