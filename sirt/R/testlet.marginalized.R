testlet.marginalized <-
function( tam.fa.obj=NULL , a1=NULL, d1=NULL, testlet=NULL, a.testlet=NULL, var.testlet=NULL){
    
	if ( ! is.null(a1) ){
		itemlabel <- rep(1:length(a1))
					}
	if (! is.null(tam.fa.obj) ){
		mod1 <- tam.fa.obj
		a1 <- mod1$B[,2,1]
		d1 <- mod1$AXsi_[,2]
		testlets <- apply( mod1$B[,2,-1] , 1 , FUN =function(zz){ 
			l1 <- which(zz > 0 ) 
			if (length(l1) == 0 ){ l1 <- NA }
			l1
				} 
				)
		a.testlet <- mod1$B[ , 2 ,  - 1]
		a.testlet <- a.testlet[ cbind( seq(1,length(testlets) ) , testlets ) ]
		var.testlet <- diag( mod1$variance )[-1]
		testlet <- testlets
		itemlabel <- colnames(mod1$resp)
						}
    # compute marginalized item intercepts and item slopes
    k <- 16*sqrt(3) / ( 15*pi )     # multiplication constant
            # 1 / k = 1.700    
    # compute lambda_logit
    multfac <- 1 - is.na(testlet)
	testlet1 <- testlet
	testlet1[ is.na(testlet1) ] <- 1
	a.testlet[ is.na(testlet) ] <- 0
    lambda.logit <- 1 / sqrt( 1 + k^2 * multfac * a.testlet[ testlet1 ]^2 * var.testlet[ testlet1 ] )
    # compute item parameters
    dfr <- data.frame("item" = itemlabel )
    dfr$testlet <- testlet
    dfr$a1 <- a1
    dfr$d1 <- d1
    dfr$a.testlet <- a.testlet
    dfr$var.testlet <- var.testlet[ testlet ]
    # marginal parameters
    dfr$a1_marg <- dfr$a1 * lambda.logit
    dfr$d1_marg <- dfr$d1 * lambda.logit
    # output
    return(dfr)
        }
