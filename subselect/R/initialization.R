initialization <- function(mat,criterion,r){
 
#####################################
#  set parameters to default values #
#####################################

  
	if (!is.null(attr(mat,"FisherI")) && attr(mat,"FisherI")==TRUE)  {
		if (r!=1) stop("In generalized linear model problems the argument r should always be set to 1.\n")
		if (criterion=="default")  criterion <- "WALD"
		else if (criterion!="WALD" && criterion!="Wald" && criterion!="wald") 
			stop("In generalized linear model problems only the Wald criterion is currently available.\n")
	}
	else {
		if (r==0 && criterion=="default")  criterion <- "RM"
        	if (r>0  && criterion=="default")  criterion <- "TAU_2"
		if (criterion=="WALD" || criterion=="Wald" || criterion=="wald") 
			warning("Since the Wald criterion was chosen, it was assumed that the mat matrix represented\n the Fisher information resulting from the estimation of a generalized linear model.\n")
	}

        assign("criterion",criterion,pos=parent.frame())

}

