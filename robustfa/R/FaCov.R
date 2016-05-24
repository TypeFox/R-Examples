##  The S3 version
FaCov <- function (x, ...)
    UseMethod("FaCov")

# setMethod("getQuan", "FaCov", function(obj) obj@n.obs)

FaCov.formula <- function (formula, data = NULL, factors = 2, cor = FALSE, method = "mle", scoresMethod = "none", ...)
## formula and data arguments are used to construct x
## method = c("mle", "pca", "pfa"), scoresMethod = c("none", "regression", "Bartlett")
## subset, na.action ignored
{
    cl <- match.call()

    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0)
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE)
    mf$... <- NULL
    if (!is.null(factors)) mf$factors <- NULL
    if (!is.null(method)) mf$method <- NULL
    if (!is.null(scoresMethod)) mf$scoresMethod <- NULL
	if (!is.null(cor)) mf$cor <- NULL
    ## so that mf has ONLY 3 elements to be compatible with the following code "mf <- eval.parent(mf)"

    mf[[1]] <- as.name("model.frame")
    ## mf = model.frame(formula = ~., data = list(x1 = c(1531125205, 106581996.9,...
    ## mf[[1]] = model.frame of class "name"; mf[[2]] = ~. of class "formula"; mf[[3]] = data of class "data.frame"

    mf <- eval.parent(mf)
    ## this is not a 'standard' model-fitting function,
    ## so no need to consider contrasts or levels

	## if (rrcov:::.check_vars_numeric(mf)) # Unexported object imported by a ¡®:::¡¯ call
    ##     stop("Fa applies only to numerical variables")

    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    attr(mt, "intercept") <- 0
    x <- model.matrix(mt, mf)

    res <- FaCov.default(x, factors = factors, method = method, scoresMethod = scoresMethod, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("FaCov")
    res@call <- cl
    res
}

FaCov.default <- function(x, factors = 2, cor = FALSE, cov.control = CovControlMcd(), method = c("mle", "pca", "pfa"), scoresMethod = c("none", "regression", "Bartlett"), ...) 
# na.action = na.fail, trace = FALSE
{

    cl <- match.call()

    if(missing(x)){
        stop("You have to provide at least some data")
    }
    data <- as.matrix(x)
    n <- nrow(data)
    p <- ncol(data)

    if(n < p)
        stop("'FaCov' can only be used with more units than variables")
    
    if (missing(method)) method = "mle"
    if (missing(scoresMethod)) scoresMethod = "none"

    ##
    ## verify and set the input parameters: k and kmax
    ##
    # kmax <- max(min(floor(kmax), floor(n/2), rankMM(x)),1)
    # if((k <- floor(k)) < 0)
    #     k <- 0
    # else if(k > kmax) {
    #     warning(paste("The number of principal components k = ", k, " is larger than kmax = ", kmax, "; k is set to ", kmax,".", sep = ""))
    #     k <- kmax
    # }
    # if(k ! =  0)
    #     k <- min(k, p)
    # else {
    #     k <- min(kmax, p)
    #     if(trace)
    #         cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep = "")
    # }

######################################################################

    ## add the option for classic covariance estimates - if cov.control = NULL
    covx <- if(!is.null(cov.control)) restimate(cov.control, data) else Cov(data)
    covmat <- list(cov = getCov(covx), center = getCenter(covx), n.obs = covx@n.obs)

    out <- switch(method, 
                  pca = factorScorePca(x = data, factors = factors, covmat = covmat, cor = cor, scoresMethod = scoresMethod), 
                  pfa = factorScorePfa(x = data, factors = factors, covmat = covmat, cor = cor, scoresMethod = scoresMethod),
                  mle = factanal(factors = factors, covmat = covmat)) # "factanal" does not have "cor" argument
	
	if(scoresMethod !=  "none" && method ==  "mle")
	out <- computeScores(out, x = data, covmat = covmat, cor = cor, scoresMethod = scoresMethod) # "computeScores" is defined in "utils.R"

    # scores   <- predict(out, newdata = data)
    # center   <- getCenter(covx)
    # sdev     <- out$sdev
    # scores   <- scores[, 1:k, drop = FALSE]
    # loadings <- out$loadings # [, 1:k, drop = FALSE]
    # eigenvalues  <- (sdev^2)[1:k]

######################################################################
    # names(eigenvalues) <- NULL
    # if(is.list(dimnames(data)))
    # {
    #     rownames(scores) <- rownames(data)  # dimnames(scores)[[1]] <- dimnames(data)[[1]]
    # }
    # dimnames(scores)[[2]] <- paste("PC", seq_len(ncol(scores)), sep = "")
    # dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("FaCov")
    res <- new("FaCov", call = cl,
			converged = out$converged,
			loadings = out$loadings[], # class(out$loadings) == "loadings", class(out$loadings[]) == "matrix"
			communality = out$communality,
			uniquenesses = out$uniquenesses,
			cor = cor,
			covariance = out$covariance,
			correlation = out$correlation,
			usedMatrix = out$usedMatrix,
			reducedCorrelation = out$reducedCorrelation,
			criteria = out$criteria,
			factors = out$factors,
			dof = out$dof,
			method = out$method,
			scores = out$scores,
			scoresMethod = scoresMethod,
			scoringCoef = out$scoringCoef,
			meanF = out$meanF,
			corF = out$corF,
			STATISTIC = out$STATISTIC,
			PVAL = out$PVAL,
			n.obs = n,
			center = getCenter(covx),
			eigenvalues = out$eigenvalues,
			cov.control = cov.control)

    ## Compute distances and flags
    # res <- rrcov:::.distances(x, p, res) # see the function ".distances.R"
    return(res)
}
