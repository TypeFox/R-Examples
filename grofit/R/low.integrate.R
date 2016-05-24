low.integrate <-
function (x,y)
{

if(is.vector(x)==FALSE || is.vector(y)==FALSE)
     stop("low.integrate: two vectors x and y are needed !")
if(length(x)!=length(y))
     stop("low.integrate: x and y have to be of same length !")


#spline fit and computation of the maximum derivative
y.spl <- NULL
try(y.spl <- smooth.spline(x,y))
if (is.null(y.spl)==TRUE){
    warning("Spline could not be fitted to data!")
    stop("Error in low.integrate")	
}

f     <- function(t){
          p<-predict(y.spl,t)
          f<-p$y
         }

low.integrate <- integrate(f,min(x),max(x))$value

}

