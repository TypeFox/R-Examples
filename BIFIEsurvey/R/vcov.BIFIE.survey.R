
####################################################
# vcov.BIFIEsurvey
vcov.BIFIEsurvey <- function( object, type=NULL , eps=1E-10 , avoid.singul = FALSE ){
	# extract replicated parameters		
	parsres <- extract.replicated.pars( BIFIE.method = object , type=type )
	res1 <- object
	#*****			
	parsM <- parsres$parsM
	parsrepM <- parsres$parsrepM	
	parnames <- parsres$parnames
	RR <- object$RR	
	# avoid.singul <- FALSE
	if ( ( class(object)=="BIFIE.correl" ) & is.null(type) ){ 
				avoid.singul <- TRUE 
						}
	if ( ( class(object) %in% c("BIFIE.freq" , "BIFIE.crosstab") )  ){ 
				avoid.singul <- TRUE 
						}		
	if ( avoid.singul ){			
#		eps <- 1E-10
		parsrepM <- parsrepM * ( 1 + stats::runif(prod(dim(parsrepM)) , 0 , eps ) )
				}
	
	fayfac <- res1$fayfac
	NP <- nrow(parsM)
	Cdes <- matrix( 1 , ncol=NP , nrow=1 )
	Ccols <- which( colSums( abs( Cdes) ) > 0 )
	rdes <- c(0)
	# compute covariance matrices
	res0 <- .Call("bifie_comp_vcov" ,  parsM = parsM , parsrepM = parsrepM , 
						Cdes , rdes , Ccols - 1 , fayfac=fayfac ,
						PACKAGE="BIFIEsurvey")
	var_w <- res0$var_w
	var_b <- res0$var_b
	Nimp <- res1$Nimp
	
	# total variance
	var_tot <- var_w  + ( 1 + 1/Nimp ) * var_b 
	rownames(var_tot) <- colnames(var_tot) <- parnames
	if (object$NMI){	
		var_tot <- BIFIE_NMI_inference_parameters( parsM , parsrepM , fayfac ,
				RR , Nimp , object$Nimp_NMI , comp_cov = TRUE )$Tm	
						}	
	return(var_tot)
	}
#########################################################
vcov.BIFIE.correl <- function( object , type = NULL , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type=type, ... )
	return(pars)
		}
# further BIFIE functions		
vcov.BIFIE.by <- function( object , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type= NULL , ...)
	return(pars)
		}
vcov.BIFIE.derivedParameters <- function( object , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type= NULL , ... )
	return(pars)
		}
vcov.BIFIE.crosstab <- function( object , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type= NULL , ... )
	return(pars)
		}
vcov.BIFIE.freq <- function( object , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type= NULL , ...)
	return(pars)
		}
vcov.BIFIE.linreg <- function( object , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type= NULL , ... )
	return(pars)
		}
vcov.BIFIE.logistreg <- function( object , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type= NULL , ...)
	return(pars)
		}
vcov.BIFIE.univar <- function( object , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type= NULL , ...)
	return(pars)
		}
vcov.BIFIE.twolevelreg <- function( object , ... ){
    if (object$se){
		pars <- vcov.BIFIEsurvey( object=object , type= NULL , ... )
				} else {
		pars <- vcov( object$micombs )		
				}
	return(pars)
		}		
vcov.BIFIE.pathmodel <- function( object , ... ){
	pars <- vcov.BIFIEsurvey( object=object , type= NULL , ...)
	return(pars)
		}
		