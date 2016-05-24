
##########################################################
# convert a fitted tam object into a mirt object
tam2mirt <- function( tamobj ){
    est.mirt <- FALSE
	# extract intercept
	AXsi <- tamobj$AXsi
	# extract loadings
	B <- tamobj$B
	# number of dimensions
	D <- dim(B)[3]
	# extract trait distribution
	mean.trait <- tamobj$beta
	cov.trait <- tamobj$variance
	# extract data
	dat <- tamobj$resp
	# factors
	if (D==1){ factors <- "F1" }
	if (D>1){ factors <- dimnames(tamobj$B)[[3]] }	
	# lavaan syntax with fixed values
	lavsyn <- tam2mirt_fix( D , factors , B , dat , AXsi ,
		   mean.trait , cov.trait , tamobj )	
	# lavaan syntax with freed values
	lavsyn.freed <- tam2mirt_freed( D , factors , B , dat , AXsi ,
		   mean.trait , cov.trait , tamobj )		
	# pseudo-estimate model in mirt: just create mirt object structure
	res <- lavaan2mirt( dat , lavsyn , est.mirt=TRUE )	
	res$mirt@nest <- as.integer(tamobj$ic$np ) # number of estimated parameters	
	# recalculate AIC, BIC, AICc and SABIC
	res$mirt@AIC <- tamobj$ic$AIC
	res$mirt@BIC <- tamobj$ic$BIC
	res$mirt@AICc <- tamobj$ic$AICc
	res$mirt@SABIC <- tamobj$ic$aBIC
	# use theta grid from estimation in TAM
	res$mirt@Theta <- tamobj$theta	
	# output
	res$lavaan.syntax.fixed <- lavsyn
	res$lavaan.syntax.freed <- lavsyn.freed
	# res$tamobj <- tamobj
	return(res)
		}
################################################################