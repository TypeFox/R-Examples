# Version: 23-12-2012, Daniel Fischer

`plot.mdr` <- function(x, which=NULL, ...){
 
 if(is.null(which)) which <- x$fold

 if ((which==1) && (which<=x$fold)) plotThis <- "evalOne"
 if ((which==2) && (which<=x$fold)) plotThis <- "evalTwo"
 if ((which==3) && (which<=x$fold)) plotThis <- "evalThree"
 if ((which==4) && (which<=x$fold)) plotThis <- "evalfour"

 temp <- x[["mdr"]][[plotThis]][,7]
 temp <- temp[temp!=0]

 plot(density(temp,na.rm=T),main="Precision Density Plot")

 invisible()
} 