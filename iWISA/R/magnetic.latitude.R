#'@title Magnetic latitude
#'
#'@description Returns cos of magnetic latitude used in Dst. This function is used in SAIndex to adjust for the magnetic latitude.
#'
#'@param col Colatitude of a station
#'@param lam Longitude of a station
#'@seealso \code{\link{SAIndex}}
#'@export
#'
magnetic.latitude <-
function(col=68.68,lam=202.00)
{
# returns cos of magnetic latitude used in Dst
# (phi,lam)-geographic coord. of station 
# (Phi,Lam)-geomagnetic coord. of station 
# first convert to radii
phi0<-79*2*pi/360
lam0<-290*2*pi/360
phi<-(90-col)*2*pi/360
lam<-lam*2*pi/360
Phi<-asin(sin(phi)*sin(phi0)+cos(phi)*cos(phi0)*cos(lam-lam0))
Lambda<-asin(cos(phi)*sin(lam-lam0)/cos(Phi))
list(Phi=Phi*360/(2*pi),Lambda=Lambda*360/(2*pi))
#returns geomagnetic coordinates in degrees!!!
}
