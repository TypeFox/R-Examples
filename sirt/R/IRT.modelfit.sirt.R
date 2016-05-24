

###########################################################
# general model fit function for sirt objects
IRT.modelfit.sirt <- function( object , mod ){
    res <- modelfit.sirt(object)
	res$IRT.IC <- CDM::IRT.IC(object)
	res$objname <- mod
	class(res) <- paste0("IRT.modelfit." , class(object) )
	return(res)
		}
###########################################################
# summary IRT.modelfit.xxx
summary.IRT.modelfit.sirt <- function( object , ... ){
	class(object) <- "modelfit.sirt"
	summary(object)
			}
###########################################################			

###########################################################
# IRT.modelfit for objects of class rasch.mml
IRT.modelfit.rasch.mml <- function( object , ... ){
    cl <- paste(match.call())[2]	
    res <- IRT.modelfit.sirt( object , mod=cl )
	return(res)	
						}									
summary.IRT.modelfit.rasch.mml <- summary.IRT.modelfit.sirt						
###########################################################						
# IRT.modelfit smirt
IRT.modelfit.smirt <- IRT.modelfit.rasch.mml
summary.IRT.modelfit.smirt <- summary.IRT.modelfit.sirt
###########################################################
# IRT.modelfit rasch.mirtlc
IRT.modelfit.rasch.mirtlc <- IRT.modelfit.rasch.mml
summary.IRT.modelfit.rasch.mirtlc <- summary.IRT.modelfit.sirt
###########################################################
# IRT.modelfit gom
IRT.modelfit.gom <- IRT.modelfit.rasch.mml
summary.IRT.modelfit.gom <- summary.IRT.modelfit.sirt
###########################################################


###########################################################
# summary IRT.modelfit2
summary.IRT.modelfit.sirt2 <- function( object , ... ){
	class(object) <- "tam.modelfit"
	summary(object)
			}
###########################################################	

###########################################################
# modelfit rm.facets
IRT.modelfit.rm.facets <- function(object, ... ){
	mod <- paste(match.call())[2]
	data <- object$procdata$dat2.NA
	probs <- object$probs
	theta.k <- object$theta.k
	f.qk.yi <- object$f.qk.yi
	res <- modelfit.cor.poly( data , probs, theta.k , f.qk.yi )
	res$IRT.IC <- CDM::IRT.IC(object)
	res$objname <- mod
	class(res) <- paste0("IRT.modelfit." , class(object) )
	return(res)		
		}
summary.IRT.modelfit.rm.facets <- summary.IRT.modelfit.sirt2		
############################################################
# model fit rm.sdt		
IRT.modelfit.rm.sdt <- IRT.modelfit.rm.facets
summary.IRT.modelfit.rm.sdt <- summary.IRT.modelfit.sirt2		
##########################################################
