"storemapattr" <- function(x)
{
    ## Verifications
    if ((!inherits(x,"asc"))&(!inherits(x,"kasc")))
        stop("x should be a map of class asc or kasc")

    ## creates an object containing only the attributes of the arguments
    toto<-0
    if (inherits(x, "asc"))
        x<-as.kasc(list(x=x))
    toto<-getkascattr(x,toto)
    class(toto)<-"mapattr"

    ## output
    return(toto)
  }

