##' Extract i.i.d. decomposition (influence function) from model object
##'
##' Extract i.i.d. decomposition (influence function) from model object
##' @export
##' @usage
##'
##' iid(x,...)
##'
##' \method{iid}{default}(x,bread,id=NULL,folds=0,maxsize=(folds>0)*1e6,...)
##'
##' @aliases iid.default
##' @param x model object
##' @param id (optional) id/cluster variable
##' @param bread (optional) Inverse of derivative of mean score function
##' @param folds (optional) Calculate aggregated iid decomposition (0:=disabled)
##' @param maxsize (optional) Data is split in groups of size up to 'maxsize' (0:=disabled)
##' @param ... additional arguments
##' @examples
##' m <- lvm(y~x+z)
##' distribution(m, ~y+z) <- binomial.lvm("logit")
##' d <- sim(m,1e3)
##' g <- glm(y~x+z,data=d,family=binomial)
##' crossprod(iid(g))
##'
iid <- function(x,...) UseMethod("iid")

##' @export
iid.default <- function(x,bread,id=NULL,folds=0,maxsize=(folds>0)*1e6,...) {
    if (!any(paste("score",class(x),sep=".") %in% methods("score"))) {
        warning("Not available for this class")
        return(NULL)
    }
    if (folds>0) {
        if (!requireNamespace("mets",quietly=TRUE)) stop("Requires 'mets'")
        U <- Reduce("rbind",mets::divide.conquer(function(data) score(x,data=data,...),
                                                 id=id,
                                                 data=data,size=round(nrow(data)/folds)))
    } else {
        U <- score(x,indiv=TRUE,...)
    }
    n <- NROW(U)
    pp <- pars(x)
    if (!missing(bread) && is.null(bread)) {
        bread <- vcov(x)
    }
    if (missing(bread)) bread <- attributes(U)$bread
    if (is.null(bread)) {
        bread <- attributes(x)$bread
        if (is.null(bread)) bread <- x$bread
        if (is.null(bread)) {
            if (maxsize>0) {
                ff <- function(p) colSums(Reduce("rbind",mets::divide.conquer(function(data) score(x,data=data,p=p,...),
                                                                              data=data,size=maxsize)))
                I <- -numDeriv::jacobian(ff,pp,method=lava.options()$Dmethod)
            } else {
                I <- -numDeriv::jacobian(function(p) score(x,p=p,indiv=FALSE,...),pp,method=lava.options()$Dmethod)
            }
            bread <- Inverse(I)
        }
    }
    iid0 <- U%*%bread
    if (!missing(id)) {
        N <- nrow(iid0)
        if (!lava.options()$cluster.index) {
            iid0 <- matrix(unlist(by(iid0,id,colSums)),byrow=TRUE,ncol=ncol(bread))
        } else {
            iid0 <- mets::cluster.index(id,mat=iid0,return.all=FALSE)
        }
        attributes(iid0)$N <- N
    }
    colnames(iid0) <- colnames(U)
  return(structure(iid0,bread=bread))
}


##' @export
iid.multigroupfit <- function(x,...) iid.default(x,combine=TRUE,...)

##' @export
iid.matrix <- function(x,...) {
    p <- ncol(x); n <- nrow(x)
    mu <- colMeans(x,na.rm=TRUE); S <- var(x,use="pairwise.complete.obs")*(n-1)/n
    iid1 <- t(t(x)-mu)
    iid2 <- matrix(ncol=(p+1)*p/2,nrow=n)
    pos <- 0
    nn <- c()
    cc <- mu
    for (i in seq(p))
        for (j in seq(i,p)) {
            pos <- pos+1
            cc <- c(cc,S[i,j])
            iid2[,pos] <- (iid1[,i]*iid1[,j])-cc[length(cc)]
            nn <- c(nn,paste(colnames(x)[c(i,j)],collapse=lava.options()$symbols[2]))
        }
    colnames(iid1) <- colnames(x); colnames(iid2) <- nn
    names(cc) <- c(colnames(iid1),colnames(iid2))
    iid1[is.na(iid1)] <- 0
    iid2[is.na(iid2)] <- 0
    structure(cbind(iid1/n,iid2/n),
              coef=cc,
              mean=mu, var=S)
}

##' @export
iid.numeric <- function(x,...) {
    n <- length(x)
    mu <- mean(x); S <- var(x)*(n-1)/n
    iid1 <- t(t(x)-mu)
    structure(cbind(mean=iid1/n,var=(iid1^2-S)/n),coef=c(mean=mu,var=S),mean=mu,var=S)
}


##' @export
iid.data.frame <- function(x,...) {
    if (!all(apply(x[1,,drop=FALSE],2,function(x) inherits(x,c("numeric","integer")))))
        stop("Don't know how to handle data.frames of this type")
    iid(as.matrix(x))
}
