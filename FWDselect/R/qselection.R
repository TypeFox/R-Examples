#'Selecting variables for several subset sizes
#'
#'@description Function that enables to obtain the best variables for more than
#'  one size of subset. Returns a table with the chosen covariates to be
#'  introduced into the models and their information criteria. Additionally, an
#'   asterisk is shown next to the size of subset which minimizes the
#'   information criterion.
#'@param x A data frame containing all the covariates.
#'@param y A vector with the response values.
#'@param qvector A vector with more than one variable-subset size to be
#'  selected.
#'@param criterion The information criterion to be used.
#'  Default is the deviance. Other functions provided
#'  are the coefficient of determination (\code{"R2"}), the residual
#'  variance (\code{"variance"}), the Akaike information criterion (\code{"aic"}),
#'  AIC with a correction for finite sample sizes (\code{"aicc"})
#'  and the Bayesian information criterion (\code{"bic"}). The deviance,
#'  coefficient of determination and variance are calculated by cross-validation.
#'@param method A character string specifying which regression method is used,
#'  i.e., linear models (\code{"lm"}), generalized additive models
#'  (\code{"glm"}) or generalized additive models (\code{"gam"}).
#'@param family A description of the error distribution and link function to be
#'  used in the model: (\code{"gaussian"}), (\code{"binomial"}) or
#'  (\code{"poisson"}).
#'@param nfolds Number of folds for the cross-validation procedure, for
#'\code{deviance}, \code{R2} or \code{variance} criterion.
#'@param cluster A logical value. If  \code{TRUE} (default), the
#'  procedure is  parallelized. Note that there are cases without enough
#'  repetitions (e.g., a low number of initial variables) that R will gain in
#'  performance through serial computation. R takes time to distribute tasks
#'  across the processors also it will need time for binding them all together
#'  later on. Therefore, if the time for distributing and gathering pieces
#'  together is greater than the time need for single-thread computing, it does
#'  not worth parallelize.
#'@param ncores An integer value specifying the number of cores to be used
#' in the parallelized procedure. If \code{NULL} (default), the number of cores to be used
#' is equal to the number of cores of the machine - 1.
#'@return
#'\item{q}{A vector of subset sizes.}
#'\item{criterion}{A vector of Information criterion values.}
#'\item{selection}{Selected variables for each size.}
#'@author Marta Sestelo, Nora M. Villanueva and Javier Roca-Pardinas.
#'@seealso \code{\link{selection}} \code{\link{plot.qselection}}.
#' @examples
#' library(FWDselect)
#' data(diabetes)
#' x = diabetes[ ,2:11]
#' y = diabetes[ ,1]
#' obj2 = qselection(x, y, qvector = c(1:9), method = "lm", criterion = "variance", cluster = FALSE)
#' obj2
#'
#'@export

qselection = function(x, y, qvector, criterion = "deviance",
    method = "lm", family = "gaussian", nfolds = 5, cluster = TRUE, ncores = NULL) {

    if (missing(x)) {
        stop("Argument \"x\" is missing, with no default")
    }
    if (missing(y)) {
        stop("Argument \"y\" is missing, with no default")
    }
    if (missing(qvector)) {
        stop("Argument \"qvector\" is missing, with no default")
    }

    in_c <- c()
    var <- c()
    qq <- c()
    cont <- 0
    res <- c()
    for (q in qvector) {
        cont <- cont + 1
        if (q == qvector[1]) {
          prevar <- NULL
        }else{
          prevar <- aux$Variable_numbers
        }
        aux <- selection(x = x, y = y, q = q, prevar = prevar, criterion = criterion,
            method = method, family = family, seconds = F, cluster = cluster, ncores = ncores)
        in_c[cont] <- round(aux$Information_Criterion,
            3)
        var[cont] <- toString(aux$Variable_names)
        qq[cont] <- q
        print(paste("Selecting subset of size",
            q, "...", sep = " "))
    }
    res <- data.frame(qq, in_c, var)
    colnames(res) <- c("q", paste(criterion), "selection")
    class(res) <- "qselection"
    return(res)
}
