#############################################
# Rcpp version of kernel PLS regression
kernelpls.fit2 <- function(X, Y, ncomp ){
	e1 <- environment()
	tsqs <- NULL
    ## Save dimnames:
    dnX <- dimnames(X)
    dnY <- dimnames(Y)
    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]
	if (nresp>1){
	    stop("PLS regression is only provided for one-dimensional responses.")
				}
    ## Center variables:
    Xmeans <- colMeans(X)
    X <- X - rep(Xmeans, each = nobj)
    Ymeans <- colMeans(Y)
    Y <- Y - rep(Ymeans, each = nobj)
	# apply Rcpp function
	res <- kernelpls_1dim(Y,X , comp=ncomp)
	.attach.environment( res=res , envir=e1 )	
	#****
	# output management copied from kernelpls.fit function
	# from the pls package
	    residuals <- - fitted + c(Y)
        fitted <- fitted + rep(Ymeans, each = nobj) # Add mean
        ## Add dimnames:
        objnames <- dnX[[1]]
        if (is.null(objnames)) objnames <- dnY[[1]]
        prednames <- dnX[[2]]
        respnames <- dnY[[2]]
        compnames <- paste("Comp", 1:ncomp)
        nCompnames <- paste(1:ncomp, "comps")
        dimnames(TT) <- dimnames(U) <- list(objnames, compnames)
        dimnames(R) <- dimnames(W) <- dimnames(P) <-
            list(prednames, compnames)
        dimnames(tQ) <- list(compnames, respnames)
        dimnames(B) <- list(prednames , nCompnames)
        dimnames(fitted) <- dimnames(residuals) <-
            list(objnames,  nCompnames)
        class(TT) <- class(U) <- "scores"
        class(P) <- class(W) <- class(tQ) <- "loadings"

        res <- list(coefficients = B,
             scores = TT, loadings = P,
             loading.weights = W,
             Yscores = U, Yloadings = t(tQ),
             projection = R,
             Xmeans = Xmeans, Ymeans = Ymeans,
             fitted.values = fitted, residuals = residuals,
             Xvar = colSums(P * P) * tsqs,
             Xtotvar = sum(X * X) )
		# R^2 measures
		R2 <- cumsum(res$Xvar) / res$Xtotvar
		R21 <- sapply( 1:ncomp , FUN = function(cc){
             1 - stats::var( Y[,1] -  res$fitted.values[,cc] ) / stats::var( Y[,1] )
            } )
		R2 <- rbind( R2 , R21)
		rownames(R2) <- c("R2(X)" , "R2(Y)")
		colnames(R2) <- compnames
		res$R2 <- R2
		class(res) <- "kernelpls.fit2"
		return(res)
		}
#######################################################
# attach all elements in a list in a local environment
.attach.environment <- function( res , envir ){
	CC <- length(res)
	for (cc in 1:CC){
		assign( names(res)[cc] , res[[cc]] , envir=envir )		
					}
					}
########################################################
# Call to Rcpp function
kernelpls_1dim <- function (Y,X,comp){ 
	res <- .Call("kernelpls_1dim_C", 
			Y,X,comp, 
			PACKAGE = "miceadds")
	return(res)
					}
##########################################################