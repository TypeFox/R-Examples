##' Missing value generator
##'
##' This function adds a binary variable to a given \code{lvm} model
##' and also a variable which is equal to the original variable where
##' the binary variable is equal to zero
##'
##' @title Missing value generator
##' @param object  \code{lvm}-object.
##' @param formula The right hand side specifies the name of a latent
##' variable which is not always observed. The left hand side
##' specifies the name of a new variable which is equal to the latent
##' variable but has missing values.  If given as a string then this
##' is used as the name of the latent (full-data) name, and the
##' observed data name is 'missing.data'
##' @param Rformula Missing data mechanism with left hand side
##' specifying the name of the observed data indicator (may also just
##' be given as a character instead of a formula)
##' @param missing.name Name of observed data variable (only used if
##' 'formula' was given as a character specifying the name of the
##' full-data variable)
##' @param suffix If missing.name is missing, then the name of the
##' oberved data variable will be the name of the full-data variable +
##' the suffix
##' @param ... Passed to binomial.lvm.
##' @return lvm object
##' @aliases Missing, Missing<-
##' @examples
##' library(lava)
##' set.seed(17)
##' m <- lvm(y0~x01+x02+x03)
##' m <- Missing(m,formula=x1~x01,Rformula=R1~0.3*x02+-0.7*x01,p=0.4)
##' sim(m,10)
##'
##'
##' m <- lvm(y~1)
##' m <- Missing(m,"y","r")
##' ## same as
##' ## m <- Missing(m,y~1,r~1)
##' sim(m,10)
##'
##' ## same as
##' m <- lvm(y~1)
##' Missing(m,"y") <- r~x
##' sim(m,10)
##'
##' m <- lvm(y~1)
##' m <- Missing(m,"y","r",suffix=".")
##' ## same as
##' ## m <- Missing(m,"y","r",missing.name="y.")
##' ## same as
##' ## m <- Missing(m,y.~y,"r")
##' sim(m,10)
##'
##' @export
##' @author Thomas A. Gerds <tag@@biostat.ku.dk>
Missing <- function(object,formula,Rformula,missing.name,suffix="0",...){
    if (is.character(Rformula)) {
        indicatorname <- Rformula
        Rformula <- toformula(Rformula,1)
    } else {
        indicatorname <- all.vars(Rformula)[1]
    }
    if (length(all.vars(formula))==1) formula <- all.vars(formula)
    if (is.character(formula)) {
        if (missing(missing.name)) missing.name <- paste0(formula,suffix)
        formula <- toformula(missing.name,formula)
    }
    newf <- update(formula,paste(".~.+",indicatorname))
    if (is.null(distribution(object,indicatorname)[[1]]) || length(list(...))>0) {
        distribution(object,indicatorname) <- binomial.lvm(...)
    }
    transform(object,newf) <- function(u){
        out <- u[,1]
        out[u[,2]==0] <- NA
        out
    }
    regression(object) <- Rformula
    object
}

##' @export
"Missing<-" <- function(object,formula,...,value) {
    Missing(object,formula,value,...)
}
