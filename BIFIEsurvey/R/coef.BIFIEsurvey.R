
###########################################
# general BIFIE method coef
coef.BIFIEsurvey <- function( object , type = NULL , ... ){	
	parsres <- extract.replicated.pars( BIFIE.method=object , type=type )
	pars <- parsres$parsM
	pars <- rowMeans(pars)
	names(pars) <- parsres$parnames
	return(pars)
		}
###############################################
# BIFIE.correl
coef.BIFIE.correl <- function( object , type = NULL , ... ){
	pars <- coef.BIFIEsurvey( object=object , type=type )
	return(pars)
		}
# further BIFIE functions		
coef.BIFIE.by <- function( object , ... ){
	pars <- coef.BIFIEsurvey( object=object , type= NULL )
	return(pars)
		}
coef.BIFIE.crosstab <- function( object , ... ){
	pars <- coef.BIFIEsurvey( object=object , type= NULL )
	return(pars)
		}		
coef.BIFIE.derivedParameters <- function( object , ... ){
	pars <- coef.BIFIEsurvey( object=object , type= NULL )
	return(pars)
		}				
coef.BIFIE.freq <- function( object , ... ){
	pars <- coef.BIFIEsurvey( object=object , type= NULL )
	return(pars)
		}				
coef.BIFIE.linreg <- function( object , ... ){
	pars <- coef.BIFIEsurvey( object=object , type= NULL )
	return(pars)
		}		
coef.BIFIE.logistreg <- function( object , ... ){
	pars <- coef.BIFIEsurvey( object=object , type= NULL )
	return(pars)
		}				
coef.BIFIE.univar <- function( object , ... ){
	pars <- coef.BIFIEsurvey( object=object , type= NULL )
	return(pars)
		}	
coef.BIFIE.twolevelreg <- function( object , ... ){
    if (object$se){
		pars <- coef.BIFIEsurvey( object=object , type= NULL )
				} else {
		pars <- coef(object$micombs)										
				}
	return(pars)
		}		
coef.BIFIE.pathmodel <- function( object , ... ){
	pars <- coef.BIFIEsurvey( object=object , type= NULL )
	return(pars)
		}			