print.BayesSAE <- 
function(x, digits = max(3, getOption("digits") - 3), ...)
{
   	cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

    m <- x$m
    p <- x$p
    beta <- as.matrix(x[[1]][,(m+1):(m+p)])
    beta <- colMeans(beta)	
    if(length(beta)) {
   	   	cat(paste("Regression coefficients in the linking model:\n", sep = ""))
   	   	print.default(format(beta, digits = digits), print.gap = 2, quote = FALSE)
   	   	cat("\n")
   	} 
   	else cat("No coefficients in the linking model \n\n")
   	invisible(x)
}

summary.BayesSAE <- 
function(object, HB = TRUE, ...)
{
   	## prediction
   	m <- object$m
    innov <- object$innov
    subset <- object$subset
   	if ((HB) && object$type != "UFH" && object$type != "UYC"){
        theta <- object$HB
        object$HB = TRUE
    }
    else{
        theta <- colMeans(as.matrix(object[[1]][,1:m]))
        object$HB <- FALSE
    }
    object$theta <- theta
	
   	## extend coefficient table
    p <- object$p
   	cf <- cbind(summary(object[[1]])[[1]][(m+1):(m+p), 1:2], summary(object[[1]])[[2]][(m+1):(m+p), c(1,5)])
   	object$cf <- cf

    ## sampling variance 	
    if (innov == "t"){     
        sig2 <- cbind(summary(object[[1]])[[1]][(m+p+2):(2*m+p+1), 1:2], 
            summary(object[[1]])[[2]][(m+p+2):(2*m+p+1), c(1,5)])
        object$sig2 <- sig2
    }
	
    ## residual standard error
    sigv <- c(summary(object[[1]])[[1]][m+p+1,1:2], summary(object[[1]])[[2]][m+p+1, c(1,5)])    
    object$sigv <- sigv
	
   	## delete some slots
    object$mcmc <- NULL    
	
   	## return
   	class(object) <- "summary.BayesSAE"
   	object
}

print.summary.BayesSAE <- 
function(x, digits = max(3, getOption("digits") - 3), ...)
{
   	cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "", sep = "\n")

	type <- x$type
	if (type == "FH")
        Type <- "Basic Fay-Herriort Model"
    else if (type == "YC")
        Type <- "You-Chanpman Model"
    else if (type == "UFH")
        Type <- "Unmatched Fay-Herriort Model"
    else if (type == "UYC")
        Type <- "Unmatched You-Chapman Model"
    else if (type == "SFH")
        Type <- "CAR Area-Level Fay-Herriot Model"
    else
        Type <- "CAY Area-Level You-Chapman Model"
	cat(Type, "\n")
    
    name <- names(x$mf[1,])[1]
    p <- x$p
    m <- x$m
	cat(paste("Sampling model: ", name, " ~ theta\n", sep = ""))

    name <- names(x$mf[1,])[2:(p)]
    if (x$tran == "F")
        cat("Linking model: theta ~ ")
    else if (x$tran == "log")
        cat("Linking model: log(theta) ~ ")
    else 
        cat("Linking model: logit(theta) ~ ")
    cat(sprintf("%s +%s", name[1:(p-2)], ""), name[p-1]) 
    if (x$spatial) 
        cat(" with spatial random effect\n")
    else
        cat("\n")

    if (x$HB) 
        cat("\nRao-Blackwellization of theta's based on the simulation:\n") 
    else
        cat("\nPosterior mean of theta's based on the simulation:\n")
    print(structure(round(as.vector(x$theta), digits = digits), .Names = 1:m))

    cat(paste("\nCoefficients in the linking model:\n", sep = ""))
    printCoefmat(x$cf, digits = digits, signif.legend = FALSE)
	
   	if (x$innov == "t") {
   	   	cat(paste("\nSampling variance in the sampling model:\n", sep = ""))
   	   	printCoefmat(x$sig2, digits = digits, signif.legend = FALSE)
   	} 

    cat(paste("\nVariance of residual in the linking model:\n", sep = ""))
    sigv <- matrix(x$sigv, c(4, 1))
    row.names(sigv) <- c("Mean", "SD", "2.5%", "97.5")	
    printCoefmat(t(sigv), digits = digits, signif.legend = FALSE)	
	
   	cat("\nDIC:", x$DIC, "\n")
     
   	invisible(x)
}

MCMC <-
function(object, ...){
    result <- object$mcmc
    return(result)
}