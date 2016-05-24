glmnet.cr <-
function (x, y, method = "backward", weights, offset=NULL,
alpha = 1, nlambda = 100, lambda.min.ratio=NULL, lambda=NULL,
standardize = TRUE, thresh = 1e-04, 
exclude, penalty.factor = NULL, maxit=100) 
{
    n <- nobs <- dim(x)[1]
    p <- m <- nvars <- dim(x)[2]
    k <- length(unique(y))
    x <- as.matrix(x)
	if (is.null(penalty.factor)) penalty.factor<-rep(1, nvars) else penalty.factor<-penalty.factor
	if (is.null(lambda.min.ratio)) lambda.min.ratio<-ifelse(nobs<nvars,0.01,0.0001)
    if (c("backward", "forward")[charmatch(method, c("backward", 
													 "forward"))] == "backward") {
        restructure <- cr.backward(x = x, y = y)
    }
    if (c("backward", "forward")[charmatch(method, c("backward", 
													 "forward"))] == "forward") {
        restructure <- cr.forward(x = x, y = y)
    }
    glmnet.data <- list(x = restructure[, -1], y = restructure[, 
						"y"])
    object <- glmnet(glmnet.data$x, glmnet.data$y, family = "binomial", weights=weights, offset=offset, alpha = alpha, nlambda = nlambda, lambda.min.ratio = lambda.min.ratio, lambda=lambda,
					 standardize = standardize, thresh = thresh, exclude=exclude, penalty.factor = c(penalty.factor,rep(0,k-1)), maxit=maxit, type.gaussian=ifelse(nvars<500,"covariance","naive"))
	object$x<-x
	object$y<-y
	object$method<-method
    class(object) <- "glmnet.cr"
    object
}

