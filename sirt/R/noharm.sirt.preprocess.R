

# e1 <- environment()	
# .attach.environment( res , envir=e1 )


############################################################
# data preprocessing noharm.sirt
.noharm.sirt.preproc <- function( dat , weights , Fpatt , Fval , 
	Ppatt , Pval , Psipatt , Psival , wgtm , dimensions ){

	res <- NULL
	res$dat0 <- dat

	#*****************
	# data processing
	res$N <- N <- nrow(dat)
	res$I <- I <- ncol(dat)
	if (is.null(weights)){ weights <- rep(1,nrow(dat) ) }
	res$weights <- weights
	dat.resp <- 1-is.na(dat)
	dat0 <- dat
	dat[ dat.resp==0] <- 0
	res$dat <- dat
	res$dat.resp <- dat.resp
	# calculate (weighted) product moment correlation
	ss <- as.matrix( crossprod( weights*dat.resp ) )
	eps <- .000001
	pm <- as.matrix( crossprod( as.matrix(dat*weights*dat.resp) ) / ( ss+eps ) )
	res$pm <- pm
	res$pm0 <- pm
	res$ss <- ss
	# CFA or EFA?
	if ( is.null(dimensions) ){ 
		model.type <- "CFA"
		modtype <- 3	# 3 - multidimensional CFA
			} else {
		model.type <- "EFA"
		modtype <- 2	# 2 - multidimensional EFA
		D <- dimensions
		Pval <- diag(D)
		Ppatt <- 0*diag(D)
		Fpatt <- matrix(1,nrow=I,ncol=D)
		if (D>1){		
			for (dd in 2:D){  Fpatt[dd,1:(dd-1)] <- 0 }
					}
		Fval <- .5*(Fpatt>0)		
		if ( D == 1 ){ 	# 1 dimension
				model.type <- "CFA" 
				modtype <- 3
						}		
					}
					
	res$model.type <- model.type
	res$modtype  <- modtype
	# initial values if they are not provided
	if ( is.null(Psival) ){ Psival <- 0*diag(res$I) }
	if ( is.null(Psipatt) ){ Psipatt <- 0*diag(res$I) }
	if ( is.null(Fval) ){ Fval <- .5*(Fpatt>0) }	
	if ( is.null(Pval) ){ 
		Pval <- diag( ncol(Ppatt) )
			}		
	# weight matrix
	wgtm.default <- FALSE
	if ( is.null(wgtm) ){ 
		wgtm <- matrix(1,I,I) 
		wgtm.default <- TRUE
					}
	diag(wgtm) <- 0					
	wgtm <- wgtm * ( ss > 0 )
	res$wgtm <- wgtm		
	res$sumwgtm <- ( sum( wgtm > 0 ) - sum( diag(wgtm) > 0 ) ) / 2	
	#***
	# column names
	D <- ncol(Ppatt)
	cn <- paste0("F",1:D)
	if (is.null(colnames(Fpatt) ) ){
		colnames(Fpatt) <- cn }
	if (is.null(colnames(Fval) ) ){
		colnames(Fval) <- colnames(Fpatt) 
					}		
	
	if (is.null(colnames(Pval))){ 
			colnames(Pval) <- colnames(Ppatt) 
					}
	if (is.null(colnames(Pval))){ 
			colnames(Pval) <- colnames(Fval) 
					}
	rownames(Pval) <- colnames(Pval)
	
	#*****
	# matrix conversion
	#****
	res$Fpatt <- as.matrix(Fpatt) 
	res$Fval <- as.matrix(Fval)
	res$Ppatt <- as.matrix(Ppatt)
	res$Pval <- as.matrix(Pval)
	res$Psipatt <- as.matrix(Psipatt)
	res$Psival <- as.matrix(Psival)
	res$Psipatt <- ( Psipatt + t(Psipatt) )/2
	res$Ppatt <- ( Ppatt + t(Ppatt) )/2
	res$D <- ncol(Ppatt)
	res$Fpatt <- 1*(res$Fpatt>0)
	res$Ppatt <- 1*(res$Ppatt>0)
	res$Psipatt <- 1*(res$Psipatt>0)
	res$wgtm.default <- wgtm.default
	return(res)
	}

######################################################################		