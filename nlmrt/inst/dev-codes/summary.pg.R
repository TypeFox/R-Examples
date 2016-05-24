summary.philipp<-function (object, ...) 
{
    smalltol <- .Machine$double.eps * 1000
    options(digits = 5)
    resname <- deparse(substitute(object))
    JJ <- object$jacobian
    res <- object$resid
    coef <- object$coeffs
    resname <- deparse(substitute(x))
    pname <- names(coef)
    npar <- length(coef)
    lo <- object$lower
    if (is.null(lo)) 
        lo <- rep(-Inf, npar)
    up <- object$upper
    if (is.null(up)) 
        up <- rep(Inf, npar)
    mi <- object$maskidx
    mt <- rep(" ", npar)
    mt[mi] <- "M"
    bdmsk <- rep(1, npar)
    bdmsk[mi] <- 0
    ct <- rep(" ", npar)
    for (i in seq_along(coef)) {
        if (lo[[i]] - coef[[i]] > 0) {
            ct[[i]] <- "-"
            if (bdmsk[[i]] == 1) 
                bdmsk[[i]] <- -3
        }
        else {
            if (coef[[i]] - lo[[i]] < smalltol * (abs(coef[[i]]) + 
                smalltol)) {
                ct[[i]] <- "L"
                if (bdmsk[[i]] != 0) 
                  bdmsk[[i]] <- -3
            }
        }
        if (coef[[i]] - up[[i]] > 0) {
            ct[[i]] <- "+"
            if (bdmsk[[i]] == 1) 
                bdmsk[[i]] <- -1
        }
        else {
            if (up[[i]] - coef[[i]] < smalltol * (abs(coef[[i]]) + 
                smalltol)) {
                ct[[i]] <- "U"
                if (bdmsk[[i]] != 0) 
                  bdmsk[[i]] <- -1
            }
        }
    }
    ss <- object$ssquares
    nobs <- length(res)
    ndof <- nobs - npar
    if (ndof <= 0) 
        stop(paste("Inadmissible degrees of freedom =", ndof, 
            sep = ""))
    sighat2 <- ss/(ndof)
    dec <- svd(JJ)
    U <- dec$u
    V <- dec$v
    Sd <- dec$d
    Sinv <- 1/Sd
    Sinv[which(bdmsk != 1)] <- 0
    VS <- crossprod(t(V), diag(Sinv))
    Jinv <- crossprod(t(VS))
    var <- Jinv * sighat2
    SEs <- sqrt(diag(var))
    gr <- crossprod(JJ, res)
    tstat <- coef/SEs
    pval <- 2 * (1 - pt(tstat, df = ndof))
    output<-data.frame(coef=coef,SE=SEs,tstat=tstat,pval=pval,stringsAsFactors=FALSE)
rownames(output)<-pname
    output[,4]<-pt(abs(output[,3]),df=length(res)-nrow(output),lower.tail=FALSE)*2
    output
}