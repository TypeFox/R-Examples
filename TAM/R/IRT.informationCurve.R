


#####################################################
# S3 method
IRT.informationCurves <- function (object, ...) {
    UseMethod("IRT.informationCurves")
       }
#####################################################


######################################################
informationCurves_mml <- function( object , h=.0001 , 
	iIndex = NULL , theta = NULL , ... ){

	eps <- 1E-20
	I <- ncol(object$resp)
	if ( is.null(iIndex) ){
		iIndex <- 1:I
				}
	resp <- object$resp
	maxcat <- apply( resp , 2 , max , na.rm=TRUE )
	A <- object$A
	AXsi <- object$AXsi
	B <- object$B
	if ( dim(B)[3] > 1 ){
		stop(paste0("Information curves can only be calculated for \n",
				"    unidimensional models.\n") )
					}
	xsi <- object$xsi$xsi
	maxK <- object$maxK
	if ( is.null( theta ) ){
		theta <- object$theta		
		theta0 <- seq( min(theta) , max(theta) , length=100 )
		theta <- matrix( theta0 , ncol=1)
						} else {
		theta <- matrix( theta , ncol=1 )
						}
						
	nnodes <- nrow(theta)
	guess <- object$item$guess

	#****
	# calculate probabilities
	if ( class(object) %in% c("tam.mml","tam.mml.2pl" ,
				"tam.mml.mfr") ){		
		p0 <- calc_prob.v5(iIndex, A, AXsi, B, xsi, theta, nnodes, 
					maxK, recalc=TRUE)$rprobs
		p1 <- calc_prob.v5(iIndex, A, AXsi, B, xsi, theta + h, 
				nnodes, maxK, recalc=TRUE)$rprobs
		p2 <- calc_prob.v5(iIndex, A, AXsi, B, xsi, theta - h, 
				nnodes, maxK, recalc=TRUE)$rprobs
								}
	if ( class(object) %in% c("tam.mml.3pl" ) ){	
		p0 <- .mml.3pl.calc_prob.v5(iIndex, A, AXsi, B, xsi, theta, 
			nnodes, maxK, recalc=TRUE , guess)
		p1 <- .mml.3pl.calc_prob.v5(iIndex, A, AXsi, B, xsi, theta - h, 
			nnodes, maxK, recalc=TRUE , guess)
		p2 <- .mml.3pl.calc_prob.v5(iIndex, A, AXsi, B, xsi, theta + h, 
			nnodes, maxK, recalc=TRUE , guess)
						}
			
	p0a <- p0
	p0[ is.na(p0) ] <- 0
	p1[ is.na(p1) ] <- 0
	p2[ is.na(p2) ] <- 0
	p0 <- p0 + eps
	
	# 1st derivative
	d1 <- ( p1 - p2  ) / ( 2 * h )    # negative sign?
	# second derivative
	# f(x+h)+f(x-h) = 2*f(x) + f''(x)*h^2
	d2 <- ( p1 + p2 - 2*p0 ) / h^2

	# calculate information curves
	IP <- dim(p0)[1]
	TP <- nnodes
	info_curves_item <- matrix( 0 , nrow= IP , ncol=TP)
	rownames(info_curves_item) <- colnames(object$resp)[ iIndex	]
	
	
	for (kk in 1:maxK){
		# kk <- 1
		info_curves_item <- info_curves_item + 
					p0[ , kk , ] * ( ( d1[,kk,] / p0[,kk,] )^2 - 
										d2[,kk,] / p0[,kk,] )
						}
	
	# test information curves
	test_info_curve <- colSums( info_curves_item )
	# standard error curve
	se_curve <- sqrt( 1 / test_info_curve )
	# category information curves
	info_curves_categories <- array( 0 ,dim=c(IP , maxK ,TP) )
	dimnames(info_curves_categories)[[1]] <- rownames(info_curves_item) 
	for (kk in 1:maxK){
		info_curves_categories[,kk,] <- info_curves_item * p0a[,kk,]
					}
	
	
	res <- list( se_curve = se_curve , test_info_curve=test_info_curve ,
				 info_curves_item = info_curves_item ,
				 info_curves_categories = info_curves_categories , 
	             theta = theta , maxK = maxK , maxcat = maxcat)
	class(res) <- "IRT.informationCurves"
	return(res)
	
		}
######################################################################
# S3 methods
IRT.informationCurves.tam.mml <- informationCurves_mml
IRT.informationCurves.tam.mml.2pl <- informationCurves_mml
IRT.informationCurves.tam.mml.3pl <- informationCurves_mml
IRT.informationCurves.tam.mml.mfr <- informationCurves_mml
######################################################################