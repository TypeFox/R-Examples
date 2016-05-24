pzmax <- function(alpha,S) {
    ##P(Zmax > z) Family wise error rate, Zmax = max |Z_i|
    k <- nrow(S)
    z <- qnorm(1-alpha/2)
    1-mets::pmvn(lower=rep(-z,k),upper=rep(z,k),sigma=cov2cor(S))
}


##' @export
p.correct <- function(object,idx,alpha=0.05) {
    S <- vcov(object); if (!missing(idx)) S <- S[idx,idx,drop=FALSE]
    f <- function(a) pzmax(a,S)-alpha
    uniroot(f,lower=0,upper=0.05)$root
}

##' Closed testing procedure
##'
##' Closed testing procedure
##' @aliases closed.testing p.correct
##' @param object estimate object
##' @param idx Index of parameters to adjust for multiple testing
##' @param null Null hypothesis value
##' @param ... Additional arguments
##' @export
##' @examples
##' m <- lvm()
##' regression(m, c(y1,y2,y3,y4,y5,y6,y7)~x) <- c(0,0.25,0,0.25,0.25,0,0)
##' regression(m, to=endogenous(m), from="u") <- 1
##' variance(m,endogenous(m)) <- 1
##' set.seed(2)
##' d <- sim(m,200)
##' l1 <- lm(y1~x,d)
##' l2 <- lm(y2~x,d)
##' l3 <- lm(y3~x,d)
##' l4 <- lm(y4~x,d)
##' l5 <- lm(y5~x,d)
##' l6 <- lm(y6~x,d)
##' l7 <- lm(y7~x,d)
##'
##' (a <- merge(l1,l2,l3,l4,l5,l6,l7,subset=2))
##' p.correct(a)
##' as.vector(closed.testing(a))
##'
closed.testing <- function(object,idx=seq_along(coef(object)),null=rep(0,length(idx)),...) {
    B <- diag(nrow=length(idx))
    e <- estimate(object,keep=idx)
    combs <- pvals <- c()
    for (i in seq_along(idx)) {
        co <- combn(length(idx),i)
        pp <- numeric(ncol(co))
        for (j in seq_along(pp)) {
            pp[j] <- compare(e,contrast=B[co[,j],,drop=FALSE],null=null[co[,j]],...)$p.value
        }
        combs <- c(combs,list(co))
        pvals <- c(pvals,list(pp))
    }
    pmax <- c()
    for (k in seq_along(idx)) {
        pk <- c()
        for (i in seq_along(idx)) {
            cols <- apply(combs[[i]],2,function(x) k%in%x)
            pk <- c(pk,pvals[[i]][which(cols)])
        }
        pmax <- c(pmax,max(pk))
    }
    return(structure(pmax,comb=combs,pval=pvals))
}
