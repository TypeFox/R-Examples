##' Extract subset of latent variable model
##'
##' Extract measurement models or user-specified subset of model
##'
##'
##' @aliases measurement
##' @param x \code{lvm}-object.
##' @param vars Character vector or formula specifying variables to include in
##' subset.
##' @param \dots Additional arguments to be passed to the low level functions
##' @return A \code{lvm}-object.
##' @author Klaus K. Holst
##' @keywords models regression
##' @examples
##'
##' m <- lvm(c(y1,y2)~x1+x2)
##' subset(m,~y1+x1)
##'
##' @export
##' @method subset lvm
subset.lvm <- function(x, vars, ...) {
    if (missing(vars)) return(x)
    if (inherits(vars,"formula")) vars <- all.vars(vars)
    if (!all(vars%in%vars(x))) stop("Not a subset of model")
    latentvars <- intersect(vars,latent(x))
    ##  g0 <- subGraph(vars, Graph(x))
    ##  res <- graph2lvm(g0)
    res <- lvm(vars)
    M <- t(x$M[vars,vars,drop=FALSE])
    for (i in seq_len(nrow(M))) {
        if (any(M[,i]==1)) {
            res <- regression(res, y=rownames(M)[M[,i]==1], x=rownames(M)[i], ...)
        }
    }
    if (length(latentvars)>0)
        latent(res) <- latentvars
    res$cov[vars,vars] <- x$cov[vars,vars]
    ## Fixed parameters:
    res$par[vars,vars] <- x$par[vars,vars]
    res$fix[vars,vars] <- x$fix[vars,vars]
    res$covpar[vars,vars] <- x$covpar[vars,vars]
    res$covfix[vars,vars] <- x$covfix[vars,vars]
    res$mean[vars] <- x$mean[vars]
    res$attributes <- x$attributes
    for (i in seq_along(x$attributes)) {
        val <- x$attributes[[i]]
        if (length(val)>0) {
            val <- val[intersect(vars,names(val))]
            res$attributes[[i]] <- val
        }
    }
    index(res) <- reindex(res)
    return(res)
}
