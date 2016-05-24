leafsfirst.new<-function(pcf=NULL, lev=NULL, refe=NULL, type="lst",
levmet="radius", ordmet="etaisrec", ngrid=NULL,
dendat=NULL, rho=0, propor=NULL, dist.type="euclid")
{
# pcf is a piecewise constant object
# type= "lst"/"shape"
# levmet= "radius"/"proba"

if ((!is.null(lev)) || (!is.null(propor))) type<-"shape"
if (!is.null(dendat)) type<-"tail"

if (type=="tail") 

lst<-leafsfirst.tail(dendat=dendat, rho=rho, refe=refe, dist.type=dist.type)


return(lst)
}


