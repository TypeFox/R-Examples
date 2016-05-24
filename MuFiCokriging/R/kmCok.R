##****************************************************************************
## Code by Loic Le Gratiet and Olivier Roustant, David Ginsbourger....
##
## slightly re-factored by Yves to illustrate the S4 mechanism.
## 
##
##***************************************************************************** 


##*****************************************************************************
## Define the new class 'kmCok' extending 'km'
## by simply giving a representation for the NEW slots
##
## The validation function can (must)  be enhanced
## 
##***************************************************************************** 

setClass("kmCok", 		
	representation( 
		d = "integer",						   ## spatial dimension
		n = "integer",						   ## observations number
			## data
		X = "matrix",				      	   ## the design of experiments, size nxd
		y = "matrix", 					      ## the observations, size nx1
			## trend information
		p = "integer",						   ## 1+number of trend basis functions
		F = "matrix",						   ## the experimental matrix, size nxp
		trend.formula = "formula",		   ## trend form
		trend.coef = "numeric",			   ## trend coefficients, size px1
			## covariance
		covariance = "covKernel",  ## covariance structure (new S4 class, see covStruct.R)
			## noisy observations
		noise.flag = "logical",  	   	   ## Are observations noisy ? 
		noise.var = "numeric",		      ## vector of length n
			## model information
		known.param = "character",		   ## known parameters: "None", "All" or "Trend"
		case = "character",				   ## "NoNugget" : deterministic observations, no nugget effect
												   ## "1Nugget"  : homogenous nugget effect, to be estimated
												   ## "Nuggets"  : known nugget or (exclusive) noisy observations 
			## optimisation 
		param.estim = "logical",			   ## Are parameter estimated ??
		method = "character", 			   ## Statistical criterion : "MLE" or "PMLE"
		penalty = "list",		   			   ## fun ("SCAD"), fun.derivative, value
		optim.method = "character",		   ## Optimisation algorithm : "BFGS" or "gen"
		lower = "numeric",					   ##
		upper = "numeric",					   ## boudaries for parameter estimation. Length covariance@param.n
		control = "list",					   ## pop.size, wait.generations, max.generations, BFGSburnin 
		gr = "logical",						   ## is analytical gradient to be used ?
		call = "language", 				   ## user call
		parinit = "numeric", 				   ## initial values used (given by user or computed)
		logLik = "numeric",				   ## objective function value (concentrated log-Likelihood)
			## auxilliary variables 
		T = "matrix",                     ## Upper triang. factor of the Choleski dec. 
												   ## of the cov. matrix : t(T)%*%T = C. Size nxn
		z = "numeric",						   ## t(T)^(-1)*(y - F%*%beta). Size nx1
		M = "matrix",							   ## t(T)^(-1)%*%F. Size nxp
		AR.p = "integer",
		AR.formula = "formula",
		AR.Z = "numeric",
		AR.F = "matrix"
	), 
	validity = function(object) {
		if (object@n <= object@d) stop("the number of experiments must be larger than the spatial dimension")
		
		if (ncol(object@y) != 1) stop("the response must be 1-dimensional")
		
		if (!identical(nrow(object@X), nrow(object@y))) stop("the number of observations is not equal to the number of experiments")
		
		if (object@covariance@nugget.flag & object@noise.flag) stop("'nugget' and 'noise' cannot be specified together")
		
		if (object@noise.flag) {	
			if (!identical(length(object@noise.var), object@n)) stop("the length of the vector 'noise.var' must be equal to the number of experiments")
		}
		
		if (!identical(object@trend.coef, numeric(0))) {
			if (!identical(length(object@trend.coef), (object@p+object@AR.p))) stop("the number of trend coefficients is not compatible with the trend formula")
		}
		if ( object@AR.p > 0 ) {
             if ( length(object@AR.Z) != object@n ) {
               return("'AR.Z' must be numeric of length 'n'")
             }
             if ( nrow(object@AR.F) != object@n ) {
               return("'AR.F' must have  'n' rows")
             }
             if ( ncol(object@AR.F) != object@AR.p ) {
               return("'AR.F' must have  'AR.p' cols")
             }
           }
	}
)



##*****************************************************************************
## Construct a 'kmCok' object as in 'km' but with a few more slots AR (4 new
## slots), filled from two new args.
## 
##
##*****************************************************************************

`kmCok` <-
 function(formula = ~1,
           design,
           response,
           formula.rho = ~1,      
           Z = NULL,              
           covtype = "matern5_2",  
           coef.trend = NULL,
           coef.cov = NULL,
           coef.var = NULL,
           nugget = NULL, nugget.estim = FALSE,
           noise.var = NULL,
	     estim.method = "MLE",
           penalty = NULL, 
           optim.method = "BFGS",
           lower = NULL, upper = NULL,
           parinit = NULL, control = NULL, gr = TRUE, iso = FALSE,
           scaling = FALSE, knots = NULL) {
    
    model <- new("kmCok")	
    model@call <- match.call()
      
    data <- data.frame(design)
    model@trend.formula <- formula <- drop.response(formula, data = data)
    F <- model.matrix(formula, data = data)
    
    X <- as.matrix(design)
    y <- as.matrix(response)
    
    model@X <- X
    model@y <- y
    model@d <- ncol(X)
    model@n <- nrow(X)
    
    ## Is there an autoregressive part? then parse the formula
    
    if (length(Z) > 0L) {
      model@AR.Z <- Z
      model@AR.formula <- drop.response(formula.rho, data = data)
      AR.F <- model.matrix(model@AR.formula, data = data)
      ## mutiply row i by Z[i] for i in 1:n
      model@AR.F <- sweep(x = AR.F, MARGIN = 1, STATS = Z, FUN = "*")
      model@AR.p <- ncol(model@AR.F)
      model@F <- cbind(model@AR.F, F)
    } else {
      if (!identical(formula.rho, ~1)) {
        warning("'formula.rho' ignored since 'Z' has length 0")
      }
      model@AR.formula <- ~1    ##
      model@AR.Z <- numeric(0)
      model@AR.F <- matrix(NA, nrow = model@n , ncol = 0L)
      model@AR.p <- 0L
      model@F <- F
    }
    model@p <- ncol(F)

    
    model@noise.flag <- (length(noise.var) != 0)
    model@noise.var <- as.numeric(noise.var)
    
    isTrend <- length(coef.trend) != 0
    isCov <- length(coef.cov) != 0
    isVar <- length(coef.var) != 0

    known.param <- ( (length(coef.trend) != 0) && (length(coef.cov) != 0) && (length(coef.var) != 0) )
    if ((isTrend && isCov && isVar) || (covtype == "covUser")) {
        known.param <- "All"
        nugget.estim <- FALSE
    }else if ((isTrend) && ((!isCov) || (!isVar))) {
        known.param <- "Trend"
    }else if ((!isTrend) && isCov && isVar) {
        known.param <- "CovAndVar"
        nugget.estim <- FALSE
    }else {
        known.param <- "None"
        coef.var <- coef.cov <- NULL
    }
    if (isCov) {
        known.covparam <- "All"
    }else {
        known.covparam <- "None"
    }
    
    model@covariance <-
      covStruct.create(covtype = covtype,
                       d = model@d, 
                       known.covparam = known.covparam,
                       var.names = colnames(X), 
                       coef.cov = coef.cov,
                       coef.var = coef.var,
                       nugget = nugget, 
                       nugget.estim = nugget.estim,
                       nugget.flag = ((length(nugget) != 0) || nugget.estim),  
                       iso = iso,
                       scaling = scaling,
                       knots = knots)
    
    ## Now, at least some parameters are unknown
    model@known.param <- known.param
    if (known.param == "All") {
        model@trend.coef <- as.numeric(coef.trend)
        model@param.estim <- FALSE
        validObject(model)
        model <- computeAuxVariables(model)
        return(model)
    }
    if (known.param == "CovAndVar") {
        model@param.estim <- TRUE
        validObject(model)
        model <- computeAuxVariables(model)
        x <- backsolve(t(model@T), model@y, upper.tri = FALSE)
        beta <- compute.beta.hat(x = x, M = model@M)
        z <- compute.z(x = x, M = model@M, beta = beta)
        model@z <- z
        model@trend.coef <- beta
        return(model)
    }
    if (known.param == "Trend") {
        model@trend.coef <- as.numeric(coef.trend)
    }    
    if (length(penalty) == 0) {
        if (is.element(estim.method, c("MLE", "LOO"))) {
            model@method <- estim.method
        }
        else {
            stop("estim.method must be: 'MLE' or 'LOO'")
        }
    } else {
        if (covtype != "gauss") {
            stop("At this stage, Penalized Maximum Likelihood is coded only for Gaussian covariance")
        }
        penalty.set <- c("SCAD")
        if (!is.element(penalty$fun, penalty.set)) {
            stop("At this stage, the penalty #function has to be one of : SCAD")
        }
        if (length(penalty$value) == 0) {
            penalty$value <- sqrt(2 * log(model@n)/model@n) * 
                seq(from = 1, by = 0.5, length = 15)
        }
        penalty$fun.derivative <- paste(penalty$fun, ".derivative", 
            sep = "")
        model@penalty <- penalty
        model@method <- "PMLE"
    }
    model@param.estim <- TRUE
    model@optim.method <- as.character(optim.method)
    if ((length(lower) == 0) || (length(upper) == 0)) {
        bounds <- covParametersBounds(model@covariance, design)
        if (length(lower) == 0) 
            lower <- bounds$lower
        if (length(upper) == 0) 
            upper <- bounds$upper
    }
    model@lower <- as.numeric(lower)
    model@upper <- as.numeric(upper)
    model@parinit <- as.numeric(parinit)
    if (optim.method == "BFGS") {
        if (length(control$pop.size) == 0) 
            control$pop.size <- 20
        if (identical(control$trace, FALSE)) 
            control$trace <- 0
        if ((length(control$trace) == 0) || (identical(control$trace, 
            TRUE))) {
            control$trace <- 3
        }
    }
    if (optim.method == "gen") {
        d <- ncol(design)
        if (length(control$pop.size) == 0) 
            control$pop.size <- min(20, floor(4 + 3 * log(d)))
        if (length(control$max.generations) == 0) 
            control$max.generations <- 5
        if (length(control$wait.generations) == 0) 
            control$wait.generations <- 2
        if (length(control$BFGSburnin) == 0) 
            control$BFGSburnin <- 0
        if (identical(control$trace, FALSE)) {
            control$trace <- 0
        }
        else control$trace <- 1
    }
    upper.alpha <- control$upper.alpha
    if (length(upper.alpha) == 0) {
        control$upper.alpha <- 1 - 1e-08
    }else if ((upper.alpha < 0) || (upper.alpha > 1)) {
        control$upper.alpha <- 1 - 1e-08
    }
    model@control <- control
    model@gr <- as.logical(gr)
    envir.logLik <- new.env()
    validObject(model, complete = TRUE)
    if ((length(noise.var) != 0) || ((length(nugget) != 0) && 
        (!nugget.estim))) {
        model@case <- "Nuggets"
    }
    if ((length(nugget) == 0) && (!nugget.estim) && (length(noise.var) == 
        0)) {
        model@case <- "NoNugget"
    }
    if ((length(noise.var) == 0) && (nugget.estim)) {
        model@case <- "1Nugget"
    }
    if ((model@method == "LOO") & (model@case != "NoNugget")) {
        stop("leave-One-Out is not available for this model")
    }
    f <- kmEstimate
    if (identical(model@method, "PMLE")) {
        cv <- function(lambda, object, f) {
            object@penalty$value <- lambda
            object@control$trace <- 0
            object <- f(object, envir = envir.logLik)
            criterion <- sum((object@y - leaveOneOut.km(object, 
                type = "UK")$mean)^2)
            return(criterion)
        }
        lambda.val <- model@penalty$value
        nval <- length(lambda.val)
        u <- rep(0, nval)
        for (i in 1L:nval) {
            u[i] <- cv(lambda.val[i], object = model, f)
        }
        plot(lambda.val, u)
        lambda <- lambda.val[which.min(u)]
        model@penalty$value <- lambda
        model <- f(model, envir = envir.logLik)
    }else {
        model <- f(model, envir = envir.logLik)
    }        
    return(model)
    
  }

##*****************************************************************************
## Prediction function that will be registred as S4 method
##
## Loic's code with a few changes
##
##*****************************************************************************

`predict.kmCok` <-
  function (object,
            newdata,
            newZ,          
            type,
            se.compute = TRUE,
            cov.compute = FALSE, 
            checkNames = FALSE,
            ...) {

    ## Write some checks here: 'newZ' should be a vector
    ## with length in accordance with 'model', ...
    ## 
    
    X <- object@X
    y <- object@y
    
    if (checkNames) {
      newdata <- checkNames(X1 = X, X2 = newdata, X1.name = "the design", 
                            X2.name = "newdata")
    } else {
      newdata <- as.matrix(newdata)
      d.newdata <- ncol(newdata)
      if (d.newdata != object@d) 
        stop("newdata must have the same numbers of columns than the experimental design")
      if (!identical(colnames(newdata), colnames(X))) {
        colnames(newdata) <- colnames(X)
      }
    }
    
    T <- object@T
    z <- object@z
    M <- object@M
    beta <- object@trend.coef
    F.newdata <- model.matrix(object@trend.formula, data = data.frame(newdata))

    Fr.newdata <- model.matrix(object@AR.formula, data = data.frame(newdata))
    dr <- dim(Fr.newdata)[2]
    Fr.newdata <- sweep(x = Fr.newdata, MARGIN = 1, STATS = newZ, FUN = "*")

    
    F.newdata <- cbind(Fr.newdata, F.newdata)
    y.predict.trend <- F.newdata %*% beta
    
    c.newdata <- covMat1Mat2(object@covariance,
                             X1 = X, X2 = newdata, 
                             nugget.flag = object@covariance@nugget.flag)
    
    Tinv.c.newdata <- backsolve(t(T), c.newdata, upper.tri = FALSE)
    y.predict.complement <- t(Tinv.c.newdata) %*% z
    y.predict <- y.predict.trend + y.predict.complement
    y.predict <- as.numeric(y.predict)
    output.list <- list(mean = y.predict, c = c.newdata, Tinv.c = Tinv.c.newdata)

    if (se.compute) {     
      s2.predict.1 <- apply(Tinv.c.newdata, 2, crossprod)

      if (object@covariance@nugget.flag) {
        total.sd2 <- object@covariance@sd2 + object@covariance@nugget
      } else total.sd2 <- object@covariance@sd2

      if (type == "SK") {
        s2.predict <- pmax(total.sd2 - s2.predict.1, 0)
        s2.predict <- as.numeric(s2.predict)
        q95 <- qnorm(0.975)
      } else if (type == "UK") {
        T.M <- chol(t(M) %*% M)
        s2.predict.mat <- backsolve(t(T.M), t(F.newdata - 
                                              t(Tinv.c.newdata) %*% M), upper.tri = FALSE)
        s2.predict.2 <- apply(s2.predict.mat, 2, crossprod)
        s2.predict <- pmax(total.sd2 - s2.predict.1 + s2.predict.2, 
                           0)
            s2.predict <- as.numeric(s2.predict)
        s2.predict <- s2.predict * object@n/(object@n - object@p)
        q95 <- qt(0.975, object@n - object@p)
      }
      lower95 <- y.predict - q95 * sqrt(s2.predict)
      upper95 <- y.predict + q95 * sqrt(s2.predict)
      output.list$sd <- sqrt(s2.predict)
      output.list$lower95 <- lower95
      output.list$upper95 <- upper95
    }
    
    if (cov.compute) {     
      if (object@covariance@nugget.flag) {
        total.sd2 <- object@covariance@sd2 + object@covariance@nugget
      }else total.sd2 <- object@covariance@sd2
      C.newdata <- covMatrix(object@covariance, newdata)[[1]]
      cond.cov <- C.newdata - crossprod(Tinv.c.newdata)
      
      if (type == "UK") {
        T.M <- chol(t(M) %*% M)
        s2.predict.mat <-
          backsolve(t(T.M),
                    t(F.newdata -  t(Tinv.c.newdata) %*% M), upper.tri = FALSE)
        cond.cov <- cond.cov + crossprod(s2.predict.mat)
        cond.cov <- cond.cov * object@n/(object@n - object@p)
      }
      
      output.list$cov <- cond.cov

    }
    
    return(output.list)
    
}

##*****************************************************************************
## Register the 'predict.km' function as a S4 'predict' method for the class
## "kmCok".
##
## See section 10.1 to 10.4 of
##
## John Chambers
## "Software for Data Analysis, programming with R".
## Springer 2008. Statistics and Computing series.
##
##
##*****************************************************************************

#if(!isGeneric("predict")) {
  setGeneric(name = "predict",
             def = function(object, ...) standardGeneric("predict")
             )
#}

## in 'signature' you must give a named character
## vector with the class of the generic function args
##

setMethod(f = "predict",
          signature = signature(object = "kmCok"), 
          definition = function(object,
            newdata,
            newZ,             
            type,
            se.compute = TRUE,
            cov.compute = FALSE,
            checkNames = FALSE,
            ... ) {
            predict.kmCok(object = object,
                          newdata = newdata,
                          newZ = newZ,
                          type = type,
                          se.compute = se.compute,
                          cov.compute = cov.compute,
                          checkNames = checkNames,
                          ...)
          })
