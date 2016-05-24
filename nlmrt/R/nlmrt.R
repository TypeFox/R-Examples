# nlmrt.R -- print and summary methods
#   result <- list(resid = resbest, jacobian = Jac, feval = feval, 
#        jeval = jeval, coefficients = pnum, ssquares = ssbest)
#    class(result) <- "nlmrt"


summary.nlmrt <- function(object, ...) {
    sumnlmrt<-list() # set up the stub
    smalltol <- .Machine$double.eps * 1000
    options(digits = 5) # 7 is default
    resname <- deparse(substitute(object))
    JJ <- object$jacobian
    res <- object$resid
    coeff <- object$coefficients
    pname<-names(coeff)
    npar <- length(coeff)
    lo <- object$lower
    if (is.null(lo)) lo <- rep( -Inf, npar)
    up <- object$upper
    if (is.null(up)) up <- rep( Inf, npar)
    mi <- object$maskidx
    mt <- rep(" ",npar) # start with all "unmasked"
    mt[mi] <- "M" # Put in the masks
    bdmsk <- rep(1, npar) # bounds and masks indicator ?? should it be 1L
    bdmsk[mi] <- 0 # masked
    ct <- rep(" ",npar) # start with all "free"
    for (i in seq_along(coeff)){
       if (lo[[i]] - coeff[[i]] > 0) {
          ct[[i]] <- "-" # lower bound violation
          if (bdmsk[[i]] == 1) bdmsk[[i]] <- -3
       } else { 
          if (coeff[[i]] - lo[[i]] < smalltol*(abs(coeff[[i]])+smalltol) ) {
             ct[[i]] <- "L" # "at" lower bound
             if (bdmsk[[i]] != 0) bdmsk[[i]] <- -3 # leave mask indication intact
          }
       }
       if (coeff[[i]] - up[[i]] > 0) {
          ct[[i]] <- "+" # upper bound violation
          if (bdmsk[[i]] == 1) bdmsk[[i]] <- -1
       } else { 
          if (up[[i]] - coeff[[i]] < smalltol*(abs(coeff[[i]])+smalltol) ) {
             ct[[i]] <- "U" # "at" upper bound
             if (bdmsk[[i]] != 0) bdmsk[[i]] <- -1 # leave mask indication intact
          }
       }
    }
    ss <- object$ssquares
    nobs <- length(res)
    ndof <- nobs - npar
    if (ndof <= 0) {
          if (ndof < 0) { stop(paste("Inadmissible degrees of freedom =",ndof,sep='')) }
          else { sighat2 <- Inf }
    } else {
       sighat2 <- as.numeric(ss)/(ndof) # 160302 NOT an array
    }
    dec <- svd(JJ)
    U <- dec$u
    V <- dec$v
    Sd <- dec$d
    if (min(Sd) <= smalltol * max(Sd)) { # singular
       SEs <- rep(NA, npar) # ?? Inf or NA
    } else {
       Sinv <- 1/Sd
       Sinv[which(bdmsk != 1)] <- 0
       if (npar > 1) {
           VS <- crossprod(t(V), diag(Sinv))
       } else {
           VS <- V*Sinv
       }
       Jinv <- crossprod(t(VS))
       var <- Jinv * sighat2
       SEs <- sqrt(diag(var))
    }
    gr <- crossprod(JJ, res)
    if (any(is.na(SEs))) {
        tstat<-rep(NA, npar)
    } else {
        if (any(SEs == 0)){
           tstat <- rep(0, npar)
        } else {
           tstat <- coeff/SEs
        }
    }
    pval<-2*(1-pt(abs(tstat), df=ndof)) # This will carry through NAs
    object<-list(resname=resname, ssquares=ss, nobs=nobs, coeff=coeff, ct=ct, mt=mt, 
           SEs=SEs, tstat=tstat, pval=pval, Sd=Sd, gr=gr, jeval=object$jeval,
           feval=object$feval)
##? LEAVE OUT: JJ  res 
##? Sd
##? gr
    invisible(object)
}

# ?? coef() function
coef.nlmrt <- function(object, ...) {
       out <- object$coefficients
       # print(object$coefficients)
       attr(out,"pkgname")<-"nlmrt"
       invisible(out)
}

print.nlmrt <- function(x, ...) {
    xx<-summary(x)
    with(xx, { 
	cat("nlmrt class object:",resname,"\n")
	pname<-names(coeff)
	npar <- length(coeff)
        cat("residual sumsquares = ",ssquares," on ",nobs,"observations\n")
        cat("    after ",jeval,"   Jacobian and ",feval,"function evaluations\n")
        cat("  name     ","      coeff    ","     SE   ","   tstat  ",
             "   pval  ","   gradient  "," JSingval  ","\n")
        for (i in seq_along(coeff)){
            cat(format(pname[i], width=10)," ")
            cat(format(coeff[[i]], digits=6, width=12))
            cat(ct[[i]],mt[[i]]," ")
            cat(format(SEs[[i]], digits=4, width=9)," ")
            cat(format(tstat[[i]], digits=4, width=9)," ")
            cat(format(pval[[i]], digits=4, width=9)," ")
            cat(format(gr[[i]], digits=4, width=10)," ")
            cat(format(Sd[[i]], digits=4, width=10)," ")
            cat("\n")
        }
    }) # remember to close with()
    invisible(x)
}

