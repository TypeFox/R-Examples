#' Post-inference for clustered data
#' 
#' Correct standard-errors to account for clustered data, doing either a degrees of freedom correction or using a heteroskedasticidty-cluster robust covariance matrix
#' possibly on the range specified by bandwidth
#' @param object Object of class lm, from which rdd_reg also inherits.
#' @param clusterVar The variable containing the cluster attributions. 
#' @param vcov. Specific covariance function to pass to coeftest. See help of sandwich
#' @param type The type of cluster correction to use: either the degrees of freedom, or a HC matrix. 
#' @param \ldots Further arguments passed to coeftest
#' @return The output of the coeftest function, which is itself of class \code{coeftest}
#' @seealso \code{\link{vcovCluster}}, which implements the cluster-robust covariance matrix estimator used by \code{cluserInf}
#' @references Wooldridge (2003) Cluster-sample methods in applied econometrics. 
#' \emph{AmericanEconomic Review}, 93, p. 133-138
#' @export
#' @import sandwich
#' @import lmtest
#' @examples
#' data(house)
#' house_rdd <- rdd_data(y=house$y, x=house$x, cutpoint=0)
#' reg_para <- rdd_reg_lm(rdd_object=house_rdd)
#' 
#' # here we just generate randomly a cluster variable:
#' nlet <- sort(c(outer(letters, letters, paste, sep='')))
#' clusRandom <- sample(nlet[1:60], size=nrow(house_rdd), replace=TRUE)
#'
#' # now do post-inference:
#' clusterInf(reg_para, clusterVar=clusRandom)
#' clusterInf(reg_para, clusterVar=clusRandom, type='HC')


clusterInf <- function(object, clusterVar, vcov. = NULL, type = c("df-adj", "HC"), ...) {
    
    if (is.null(clusterVar)) 
        stop("clusterVar seems to be NULL?")
    type <- match.arg(type)
    
    if (type == "df-adj") {
        nClus <- if (is.factor(clusterVar)) 
            nlevels(clusterVar) else length(unique(clusterVar))
        res <- coeftest(object, vcov. = vcov., df = nClus, ...)
    } else {
        if (!is.null(vcov.)) 
            warning("arg 'vcov.' not used when 'type=HC' (default vcovCluster used)")
        res <- coeftest(object, vcov. = function(x) vcovCluster(x, clusterVar = clusterVar), ...)
    }
    
    return(res)
}

#' @export
estfun.rdd_reg_np <- function(x, ...) {
    inf_met <- infType(x)  ## def in Misc.R
    if (inf_met == "se") 
        stop("No 'vcovHC', 'vcovCluster', 'estfun' etc can be applied to RDDrg_np with non-parametric inference estimators")
    estfun(x$RDDslot$model)
}

#' @export
bread.rdd_reg_np <- function(x, ...) {
    inf_met <- infType(x)  ## def in Misc.R
    if (inf_met == "se") 
        stop("No 'vcovHC', 'vcovCluster', 'estfun' etc can be applied to RDDrg_np with non-parametric inference estimators")
    bread(x$RDDslot$model)
}


# sandwich.rdd_reg_np <- function (x, bread. = bread, meat. = meat, ...){ inf_met <- infType(x) ## def in Misc.R
# if(inf_met=='se') stop('No 'vcovHC', 'vcovCluster', 'estfun' etc can be applied to RDDrg_np with non-parametric inference
# estimators') sandwich(x$RDDslot$model, bread.=bread., meat.=meat., ...)  }

#' @export
model.frame.rdd_reg_np <- function(formula, ...) model.frame(formula$RDDslot$model)

#' Cluster Heteroskedasticity-consistent estimation of the covariance matrix. 
#' 
#' Offer a cluster variant of the usual Heteroskedasticity-consistent 
#' @param object Object of class lm, from which rdd_reg also inherits.
#' @param clusterVar The variable containing the cluster attributions. 
#' @return A matrix containing the covariance matrix estimate.
#' @author Mahmood Arai, see \url{http://people.su.se/~ma/econometrics.html}
#' @references Cameron, C.,  Gelbach, J. and Miller, D. (2011) Robust Inference With Multiway Clustering,
#' \emph{Journal of Business and Economic Statistics},  vol. 29(2), pages 238-249.
#' #' @references Wooldridge (2003) Cluster-sample methods in applied econometrics. 
#' \emph{American Economic Review}, 93, p. 133-138
#' @references Arai, M. (2011) Cluster-robust standard errors using R, Note available \url{http://people.su.se/~ma/clustering.pdf}. 
#' @export
#' @seealso \code{\link{clusterInf}} for a direct function, allowing also alternative cluster inference methods. 
#' @examples
#' data(STAR_MHE)
#' if(all(c(require(sandwich), require(lmtest)))){
#' 
#' # Run simple regression:
#' reg_krug <- lm(pscore~cs, data=STAR_MHE)
#' 
#' # Row 1 of Table 8.2.1, inference with standard vcovHC:
#' coeftest(reg_krug,vcov.=vcovHC(reg_krug, 'HC1'))[2,2]
#' 
#' # Row 4 of Table 8.2.1, inference with cluster vcovHC:
#' coeftest(reg_krug,vcov.=vcovCluster(reg_krug, clusterVar=STAR_MHE$classid))[2,2]
#' }

vcovCluster <- function(object, clusterVar) {
    M <- length(unique(clusterVar))
    N <- length(clusterVar)
    K <- getModelRank(object)
    dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
    uj <- apply(estfun(object), 2, function(x) tapply(x, clusterVar, sum))
    # require('sandwich')
    dfc * sandwich::sandwich(object, meat. = crossprod(uj)/N)
}

#' @rdname vcovCluster
#' @param clusterVar1,clusterVar2 The two cluster variables for the 2-cluster case.
#' @export
vcovCluster2 <- function(object, clusterVar1, clusterVar2) {
    # R-codes (www.r-project.org) for computing multi-way clustered-standard errors. Mahmood Arai, Jan 26, 2008.  See: Thompson
    # (2006), Cameron, Gelbach and Miller (2006) and Petersen (2006).  reweighting the var-cov matrix for the within model
    
    K <- getModelRank(object)
    estF <- estfun(object)
    
    clusterVar12 <- paste(clusterVar1, clusterVar2, sep = "")
    M1 <- length(unique(clusterVar1))
    M2 <- length(unique(clusterVar2))
    M12 <- length(unique(clusterVar12))
    N <- length(clusterVar1)
    
    dfc1 <- (M1/(M1 - 1)) * ((N - 1)/(N - K))
    dfc2 <- (M2/(M2 - 1)) * ((N - 1)/(N - K))
    dfc12 <- (M12/(M12 - 1)) * ((N - 1)/(N - K))
    
    u1j <- apply(estF, 2, function(x) tapply(x, clusterVar1, sum))
    u2j <- apply(estF, 2, function(x) tapply(x, clusterVar2, sum))
    u12j <- apply(estF, 2, function(x) tapply(x, clusterVar12, sum))
    vc1 <- dfc1 * sandwich(object, meat. = crossprod(u1j)/N)
    vc2 <- dfc2 * sandwich(object, meat. = crossprod(u2j)/N)
    vc12 <- dfc12 * sandwich(object, meat. = crossprod(u12j)/N)
    vcovMCL <- vc1 + vc2 - vc12
    vcovMCL
}

getModelRank <- function(object, ...) UseMethod("getModelRank")

getModelRank.default <- function(object, ...) object$rank

getModelRank.rdd_reg_np <- function(object, ...) getModelRank.default(object$RDDslot$model) 
