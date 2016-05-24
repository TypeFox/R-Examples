

#############################################
# imputation for grouped data
mice.impute.grouped <- function (y, ry, x, low=NULL , upp=NULL ,  ...){
    x <- cbind(1, as.matrix(x))
    newstate <- get( "newstate" , pos = parent.frame() )  
    vname <- get("vname", pos = parent.frame())
    Y <- cbind( low[[vname]] , upp[[vname]] )        
    # draw bootstrap sample
    N <- nrow(Y)
    ind <- base::sample( 1:N , replace=TRUE )
    X1 <- x[ ind , -1]
    Y1 <- Y[ ind , ]
    # fit grouped model to bootstrap sample
	mod <- grouped::grouped( Y1 ~ X1 )
    beta <- coef(mod)
	sigma <- summary(mod)$sigma	
	# calculate predicted values
	ypred <- x %*% beta
	# compute quantiles
	qlow <- stats::pnorm( Y[,1] , mean=ypred , sd = sigma )
	qupp <- stats::pnorm( Y[,2] , mean=ypred , sd = sigma )
    # draw randomly a quantile
	samp_quant <- qlow + ( qupp - qlow)* stats::runif(N)
	# draw imputed value
	yimp <- stats::qnorm( samp_quant , mean=ypred , sd = sigma )
# dfr <- cbind( ypred , Y , qlow , qupp , samp_quant , yimp)
# Revalpr("round(dfr[1:30,],4)")	
    return(yimp)    
    }