seqeconstraint <- function(maxGap=-1, windowSize=-1, ageMin=-1, ageMax=-1, ageMaxEnd=-1, countMethod=1){
	## check that all constraints are coherent
	if(ageMaxEnd != -1 && ageMax == -1){
		ageMax <- ageMaxEnd
	}
	if(maxGap!= -1 && windowSize!=-1 && maxGap>windowSize){
		stop(" [!] maxGap is greater than windowSize")
	}
	if(ageMin!= -1 && ageMax!=-1 && ageMin>ageMax){
		stop(" [!] ageMin is greater than ageMax or ageMaxEnd")
	}
        if(!countMethod%in%seq(1,6,1)&
           !countMethod%in%c("COBJ","CDIST_O","CWIN","CMINWIN",
                             "CDIST"))
          {
            stop(" [!] unknown countMethod input")
          }
	ret <- list()
	ret$maxGap <- maxGap
	ret$windowSize <- windowSize
	ret$ageMin <- ageMin
	ret$ageMax <- ageMax
	ret$ageMaxEnd <- ageMaxEnd
        if (is.character(countMethod))
          {
            ret$countMethod <- switch(countMethod,
                                      COBJ=1,CDIST_O=2,CWIN=3,
                                      CMINWIN=4,CDIST=5)
          } else {
            ret$countMethod <- countMethod
          }
	class(ret) <- "seqeconstraint"
	return(ret)
}

print.seqeconstraint<-function(x, ...){
	z<-data.frame(Constraint=names(x),Value=as.numeric(x))
	z <- z[z$"Value"!=-1, ]
        if (z[z$"Constraint"=="countMethod","Value"] == 1) {
          z[z$"Constraint"=="countMethod","Value"] <-
            "COBJ"
	}
	if (z[z$"Constraint"=="countMethod","Value"] == 2) {
          z[z$"Constraint"=="countMethod","Value"] <-
            "CDIST_0"
	}
        if (z[z$"Constraint"=="countMethod","Value"] == 3) {
          z[z$"Constraint"=="countMethod","Value"] <-
            "CWIN"
	}
        if (z[z$"Constraint"=="countMethod","Value"] == 4) {
          z[z$"Constraint"=="countMethod","Value"] <-
            "CMINWIN"
	}
        if (z[z$"Constraint"=="countMethod","Value"] == 5) {
          z[z$"Constraint"=="countMethod","Value"] <-
            "CDIST"
	}
	if(nrow(z) > 0) { 
		print(z, row.names=FALSE,...) 
	}
}
