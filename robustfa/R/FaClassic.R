## setMethod("getQuan", "FaClassic", function(obj) obj@n.obs)

##  The S3 version
FaClassic <- function (x, ...) UseMethod("FaClassic")

FaClassic.formula <- function (formula, data = NULL, factors = 2, cor = FALSE, method = "mle", scoresMethod = "none", ...)
## formula and data arguments are used to construct x
## method = c("mle", "pca", "pfa"), scoresMethod = c("none", "regression", "Bartlett")
## subset, na.action ignored
{
    cl <- match.call() # class(cl) = "call"

    mt <- terms(formula, data = data)
    if (attr(mt, "response") > 0)
        stop("response not allowed in formula")
    mf <- match.call(expand.dots = FALSE) # class(mf) = "call"
    mf$... <- NULL
    if (!is.null(factors)) mf$factors <- NULL
    if (!is.null(method)) mf$method <- NULL
    if (!is.null(scoresMethod)) mf$scoresMethod <- NULL
	if (!is.null(cor)) mf$cor <- NULL
    ## so that mf has ONLY 3 elements to be compatible with the following code "mf <- eval.parent(mf)"

    mf[[1]] <- as.name("model.frame") 
    ## mf[[1]] is its name, mf is not a list
    ## mf is still a call, however, its name has changed to "model.frame"
    ## mf = model.frame(formula = ~., data = list(x1 = c(1531125205, 106581996.9,...
    ## mf[[1]] = model.frame of class "name"; mf[[2]] = ~. of class "formula"; mf[[3]] = data of class "data.frame"

    mf <- eval.parent(mf)
    ## mf = data now is a data.frame, mf = eval.parent(mf) = eval(mf, parent.frame(1)) = eval(mf) = data
    ## parent.frame(1) = R_GlobalEnv
    ## Now mf is an object of class "data.frame" with attributes:
    ## names, terms (with attributes: variables, factors, term.labels, order, intercept, response, .Environment, predvars, dataClasses), row.names, class
    ## this is not a `standard' model-fitting function,
    ## so no need to consider contrasts or levels
    
	## if (rrcov:::.check_vars_numeric(mf)) # Unexported object imported by a ¡®:::¡¯ call
    ##     stop("Fa applies only to numerical variables")

    na.act <- attr(mf, "na.action")
    mt <- attr(mf, "terms")
    ## mt is an object of class c("terms", "formula") with attributes: variables, factors, term.labels, order, intercept, response, .Environment, predvars, dataClasses
    attr(mt, "intercept") <- 0

    x <- model.matrix(mt, mf)
    ## x = as.matrix(data)
    ## x is not an object with class "matrix", and its attributes are:
    ## dim, dimnames (a list of dimnames[[1]] and dimnames[[2]]), assign

    res <- FaClassic.default(x, factors = factors, method = method, scoresMethod = scoresMethod, ...)

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("FaClassic") 
    ## cl is a call, its name is changed from "FaClassic.formula" to "FaClassic"
    res@call <- cl
    res
}

FaClassic.default <- function(x, factors = 2, cor = FALSE, method = c("mle", "pca", "pfa"), scoresMethod = c("none", "regression", "Bartlett"), ...)
{
    cl <- match.call()

    if(missing(x)){
        stop("You have to provide at least some data")
    }
    data <- as.matrix(x)
    n <- nrow(data)
    p <- ncol(data)

    if (missing(method)) method = "mle"
    if (missing(scoresMethod)) scoresMethod = "none"

#     Xsvd <- kernelEVD(data) # kernelEVD is defined in "Fa.R"
#     if(Xsvd$rank = =  0) {
#         stop("All data points collapse!")
#     }

    ##
    ## verify and set the input parameters: k and kmax
    ##
#     kmax <- max(min(floor(kmax), floor(n/2), Xsvd$rank),1)
#     if((k <- floor(k)) < 0)
#         k <- 0
#     else if(k > kmax) {
#         warning(paste("The number of principal components k = ", k, " is larger than kmax = ", kmax, "; k is set to ", kmax,".", sep = ""))
#         k <- kmax
#     }
#     if(k ! =  0)
#         k <- min(k, ncol(data))
#     else {
#         k <- min(kmax,ncol(data))
#         if(trace)
#             cat("The number of principal components is defined by the algorithm. It is set to ", k,".\n", sep = "")
#     }

#     loadings    <- Xsvd$loadings[, 1:k, drop = FALSE]
#     eigenvalues <- as.vector(Xsvd$eigenvalues[1:k])
#     center      <- as.vector(Xsvd$center)
#     scores      <- Xsvd$scores[, 1:k, drop = FALSE]

#     if(is.list(dimnames(data))) {
#         dimnames(scores)[[1]] <- dimnames(data)[[1]]
#    } else {
#         dimnames(scores)[[1]] <- 1:n
#     }
#     dimnames(scores)[[2]] <- paste("PC", seq_len(ncol(scores)), sep = "")
#     dimnames(loadings) <- list(colnames(data), paste("PC", seq_len(ncol(loadings)), sep = ""))

    ## use Cov (or CovClassic) to calculate the covariance matrix
    covx <- Cov(data)
    covmat <- list(cov = getCov(covx), center = getCenter(covx), n.obs = covx@n.obs)

    out <- switch(method, 
                  pca = factorScorePca(x = data, factors = factors, covmat = covmat, cor = cor, scoresMethod = scoresMethod),
                  pfa = factorScorePfa(x = data, factors = factors, covmat = covmat, cor = cor, scoresMethod = scoresMethod), 
                  mle = factanal(factors = factors, covmat = covmat)) # "factanal" does not have "cor" argument

    if(scoresMethod != "none" && method == "mle")
	out <- computeScores(out, x = data, covmat = covmat, cor = cor, scoresMethod = scoresMethod) # "computeScores" is defined in "utils.R"

    ## fix up call to refer to the generic, but leave arg name as `formula'
    cl[[1]] <- as.name("FaClassic")
    ## cl is a call, its name is changed from "FaClassic.default" to "FaClassic"
    res <- new("FaClassic", call = cl,
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
			cov.control = NULL)

    ## Compute distances and flags
    # res <- .distances(data, Xsvd$rank, res) # .distances is defined in "Fa.R"
    return(res)
}
