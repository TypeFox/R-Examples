.packageName <- "dr"
#####################################################################
#
#     Dimension reduction methods for R and Splus
#     Revised in July, 2004 by Sandy Weisberg and Jorge de la Vega
#     tests for SAVE by Yongwu Shao added April 2006
#     'ire/pire' methods added by Sanford Weisberg Summer 2006
#     Version 3.0.0 August 2007
#     copyright 2001, 2004, 2006, 2007, Sanford Weisberg
#     sandy@stat.umn.edu
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#####################################################################
#     dr is the primary function
#####################################################################
dr <- 
function(formula, data, subset, group=NULL, na.action=na.fail, weights,...)
{      
    mf <- match.call(expand.dots=FALSE)
    mf$na.action <- na.action
    if(!is.null(mf$group)) {
     mf$group <- eval(parse(text=as.character(group)[2]), 
                     data,environment(group)) }
    mf$... <- NULL # ignore ...
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    Y <- model.extract(mf,"response")
    X <- model.matrix(mt, mf)
    y.name <- if(is.matrix(Y)) colnames(Y) else
        as.character(attr(mt, "variables")[2])
    int <- match("(Intercept)", dimnames(X)[[2]], nomatch=0)
    if (int > 0) X <- X[, -int, drop=FALSE] # drop the intercept from X
    weights <- mf$"(weights)"               
    if(is.null(weights)){                   # create weights if not present
         weights <- rep(1,dim(X)[1])} else
        {
         if(any(weights<0))stop("Negative weights")
         pos.weights <- weights > 1.e-30
         weights <- weights[pos.weights]
         weights <- dim(X)[1]*weights/sum(weights)
         X <- X[pos.weights,]
         Y <- if(is.matrix(Y)) Y[pos.weights,] else Y[pos.weights]}
    ans <- dr.compute(X,Y,weights=weights,group=mf$"(group)",...)
    ans$call <- match.call()
    ans$y.name <- y.name
    ans$terms <- mt
    ans
}

dr.compute <-
function(x,y,weights,group=NULL,method="sir",chi2approx="bx",...)
{ 
    if (NROW(y) != nrow(x))     #check lengths
        stop("The response and predictors have different number of observations")
    if (NROW(y) < ncol(x))
        stop("The methods in dr require more observations than predictors")
    #set the class name
    classname<- if (is.matrix(y)) c(paste("m",method,sep=""),method) else 
        if(!is.null(group)) c(paste("p",method,sep=""),method) else
        method
    genclassname<-"dr"
    sweights <- sqrt(weights)
    qrz <- qr(scale(apply(x,2,function(a, sweights) a  * sweights, sweights), 
                     center=TRUE, scale=FALSE))  # QR decomp of WEIGHTED x matrix
#initialize the object and then call the fit method
    ans <- list(x=x,y=y,weights=weights,method=method,cases=NROW(y),qr=qrz,
                group=group,chi2approx=chi2approx)
    class(ans) <-  c(classname,genclassname)   #set the class
    ans <- dr.fit(object=ans,...)   # ... contains args for the method
    ans$x <- ans$x[,qrz$pivot[1:qrz$rank]] # reorder x and reduce to full rank
    ans #return the object
}
    
#######################################################
#    accessor functions
#######################################################

dr.x <- function(object) {object$x}
dr.wts <- function(object) {object$weights}
dr.qr <- function(object) {object$qr}
dr.Q <- function(object){UseMethod("dr.Q")}
dr.Q.default <- function(object){ qr.Q(dr.qr(object))[,1:object$qr$rank] }
dr.R <- function(object){UseMethod("dr.R")}
dr.R.default <- function(object){ qr.R(dr.qr(object))[1:object$qr$rank,1:object$qr$rank]}
dr.z <- function(object) { sqrt(object$cases) * dr.Q(object) }
dr.y.name <- function(object) {object$y.name}
dr.basis <- function(object,numdir) {UseMethod("dr.basis")}
dr.basis.default <- function(object,numdir=object$numdir){
 object$evectors[,1:numdir]}
dr.evalues <- function(object) {UseMethod("dr.evalues")}
dr.evalues.default <- function(object) object$evalues

#####################################################################
#     Fitting function
#####################################################################
dr.fit <- function(object,numdir=4,...) UseMethod("dr.fit")

dr.fit.default <-function(object,numdir=4,...){ 
    M <- dr.M(object,...)  # Get the kernel matrix M, method specific
    D <- if(dim(M$M)[1] == dim(M$M)[2]) eigen(M$M) else {
          if(ncol(M$M)==1) eigen(M$M %*% t(M$M)) else eigen(t(M$M) %*% M$M)} 
    or <- rev(order(abs(D$values)))
    evalues <- D$values[or]
    raw.evectors <- D$vectors[,or]  
    evectors <- backsolve(sqrt(object$cases)*dr.R(object),raw.evectors)
    evectors <- if (is.matrix(evectors)) evectors else matrix(evectors,ncol=1)
    evectors <- apply(evectors,2,function(x) x/sqrt(sum(x^2)))
    names <- colnames(dr.x(object))[1:object$qr$rank]
    dimnames(evectors)<-
         list(names, paste("Dir", 1:NCOL(evectors), sep=""))
    aa<-c( object, list(evectors=evectors,evalues=evalues, 
                numdir=min(numdir,dim(evectors)[2],dim(M$M)[1]),
                raw.evectors=raw.evectors), M)
    class(aa) <- class(object)
    return(aa)
}


#####################################################################
###
###  dr methods. Each method REQUIRES a dr.M function, and
###  may also have a dr.y function and an function method.
###
#####################################################################

#####################################################################
###  generic functions
#####################################################################

dr.M <- function(object, ...){UseMethod("dr.M")}
dr.y <- function(object) {UseMethod("dr.y")}
dr.y.default <- function(object){object$y}
dr.test <- function(object, numdir,...){  UseMethod("dr.test")}
dr.test.default <-function(object, numdir,...) {NULL}
dr.coordinate.test <- function(object,hypothesis,d,chi2approx,...)
                     {  UseMethod("dr.coordinate.test")  }
dr.coordinate.test.default <-
       function(object, hypothesis,d,chi2approx,...) {NULL}
dr.joint.test <- function(object,...){ UseMethod("dr.joint.test")}
dr.joint.test.default <- function(object,...){NULL}

#####################################################################
#     OLS
#####################################################################

dr.M.ols <- function(object,...) {
 ols <- t(dr.z(object)) %*% (sqrt(dr.wts(object)) * dr.y(object))
 return(list( M= ols %*% t(ols), numdir=1))}
       

#####################################################################
#     Sliced Inverse Regression and multivariate SIR
#####################################################################

dr.M.sir <-function(object,nslices=NULL,slice.function=dr.slices,sel=NULL,...) {
    sel <- if(is.null(sel))1:dim(dr.z(object))[1] else sel
    z <- dr.z(object)[sel,]
    y <- dr.y(object)
    y <- if (is.matrix(y)) y[sel,] else y[sel]
# get slice information    
    h <- if (!is.null(nslices)) nslices else max(8, NCOL(z)+3)
    slices <- slice.function(y,h)
    #slices<- if(is.null(slice.info)) dr.slices(y,h) else slice.info
# initialize slice means matrix
    zmeans <- matrix(0,slices$nslices,NCOL(z))
    slice.weight <- rep(0,slices$nslices)  # NOT rep(0,NCOL(z))
    wts <- dr.wts(object)
# compute weighted means within slice 
    wmean <- function (x, wts) { sum(x * wts) / sum (wts) }
    for (j in 1:slices$nslices){
      sel <- slices$slice.indicator==j
      zmeans[j,]<- apply(z[sel,,drop=FALSE],2,wmean,wts[sel])
      slice.weight[j]<-sum(wts[sel])}
# get M matrix for sir
    M <- t(zmeans) %*% apply(zmeans,2,"*",slice.weight)/ sum(slice.weight)
    return (list (M=M,slice.info=slices))
}

dr.M.msir <-function(...) {dr.M.sir(...)}

dr.test.sir<-function(object,numdir=object$numdir,...) {
#compute the sir test statistic for the first numdir directions
    e<-sort(object$evalues)
    p<-length(object$evalues)
    n<-object$cases
    st<-df<-pv<-0
    nt <- min(p,numdir)
    for (i in 0:(nt-1))
      {st[i+1]<-n*(p-i)*mean(e[seq(1,p-i)])
       df[i+1]<-(p-i)*sum(object$slice.info$nslices-i-1)
       pv[i+1]<-1-pchisq(st[i+1],df[i+1])
      }
    z<-data.frame(cbind(st,df,pv))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    dimnames(z)<-list(rr,c("Stat","df","p.value")) 
    z
}

#  Written by Yongwu Shao, May 2006
dr.coordinate.test.sir<-function(object,hypothesis,d=NULL,
    chi2approx=object$chi2approx,pval="general",...){
    gamma <- if (class(hypothesis) == "formula")
        coord.hyp.basis(object, hypothesis)
        else as.matrix(hypothesis)
    p<-length(object$evalues)
    n<-object$cases
    z <- dr.z(object)
    ev <- object$evalues
    slices<-object$slice.info
    h<-slices$nslices
    d <- if(is.null(d)) min(h,length(ev)) else d
    M<-object$M
    r<-p-dim(gamma)[2]
    H<- (dr.R(object)) %*% gamma  
    H <- qr.Q(qr(H),complete=TRUE) # a p times p matrix
    QH<- H[,1:(p-r)] # first p-r columns
    H <- H[,(p-r+1):p] # last r columns
    st<-n*sum(ev[1:d])-n*sum(sort(eigen(t(QH)%*%M%*%QH)$values,
               decreasing=TRUE)[1:min(d,p-r)])
    wts <- 1-ev[1:min(d,h-1)]
# each eigenvalue occurs r times.  
    testr <- dr.pvalue(rep(wts,r),st,chi2approx=chi2approx) 
# general test
    epsilon<-array(0,c(n,h))
    zmeans<-array(0,c(p,h)) 
    for (i in 1:h) {
        sel<-(slices$slice.indicator==i)
        f_k<-sum(sel)/n
        zmeans[,i]<-apply(z[sel,],2,mean)
        epsilon[,i]<-(sel-f_k-z%*%zmeans[,i]*f_k)/sqrt(f_k)
        }
    HZ<- z%*%H
    Phi<-svd(zmeans,nv=h)$v[,1:min(d,h)]
    epsilonHZ<-array(0,c(n,r*d))
    for (j in 1:n) epsilonHZ[j,]<-t(Phi)%*%t(t(epsilon[j,]))%*%t(HZ[j,])
    wts <- eigen(((n-1)/n)*cov(epsilonHZ))$values
    testg<-dr.pvalue(wts[wts>0],st,chi2approx=chi2approx)
    pv <- if (pval=="restricted") testg$pval.adj else testr$pval.adj
    z <- data.frame(cbind(st, pv))
    dimnames(z) <- list("Test", c("Statistic", "P.value"))
    z
}

#####################################################################
#     Sliced Average Variance Estimation
#  Original by S. Weisberg; modified by Yongwu Shao 4/27/2006
#  to compute the A array needed for dimension tests
#####################################################################

dr.M.save<-function(object,nslices=NULL,slice.function=dr.slices,sel=NULL,...) {
    sel <- if(is.null(sel)) 1:dim(dr.z(object))[1] else sel
    z <- dr.z(object)[sel,]
    y <- dr.y(object)
    y <- if (is.matrix(y)) y[sel,] else y[sel]
    wts <- dr.wts(object)[sel]
    h <- if (!is.null(nslices)) nslices else max(8, ncol(z) + 3)
    slices <- slice.function(y,h)
    #slices <- if (is.null(slice.info)) dr.slices(y, h) else slice.info
# initialize M
    M <- matrix(0, NCOL(z), NCOL(z))
# A is a new 3D array, needed to compute tests.
    A <- array(0,c(slices$nslices,NCOL(z),NCOL(z)))
    ws <- rep(0, slices$nslices)
# Compute weighted within-slice covariance matrices, skipping any slice with
# total weight smaller than 1
    wvar <- function(x, w) {
        (if (sum(w) > 1) {(length(w) - 1)/(sum(w) - 0)} else {0}) * 
                 var(sweep(x, 1, sqrt(w), "*"))}
    for (j in 1:slices$nslices) {
        ind <- slices$slice.indicator == j
        IminusC <- diag(rep(1, NCOL(z))) - wvar(z[ind, ], wts[ind])
        ws[j] <- sum(wts[ind])
        A[j,,] <- sqrt(ws[j])*IminusC  # new
        M <- M + ws[j] * IminusC %*% IminusC
    }
    M <- M/sum(ws)
    A <- A/sqrt(sum(ws)) # new
    return(list(M = M, A = A, slice.info = slices))
}

# Written by Yongwu Shao, 4/27/2006
dr.test.save<-function(object,numdir=object$numdir,...) {
    p<-length(object$evalues)
    n<-object$cases
    h<-object$slice.info$nslices
    st.normal<-df.normal<-st.general<-df.general<-0
    pv.normal<-pv.general<-0
    nt<-numdir
    A<-object$A
    M<-object$M
    D<-eigen(M)
    or<-rev(order(abs(D$values)))
    evectors<-D$vectors[,or]
    for (i in 1:nt-1) {
        theta<-evectors[,(i+1):p]
        st.normal[i+1]<-0
        for (j in 1:h) {
            st.normal[i+1]<-st.normal[i+1] + 
                            sum((t(theta) %*% A[j,,] %*% theta)^2)*n/2
        }
        df.normal[i+1]<-(h-1)*(p-i)*(p-i+1)/2
        pv.normal[i+1]<-1-pchisq(st.normal[i+1],df.normal[i+1])
# general test
        HZ<- dr.z(object) %*% theta
        ZZ<-array(0,c(n,(p-i)*(p-i)))
        for (j in 1:n) ZZ[j,]<- t(t(HZ[j,])) %*% t(HZ[j,])
        Sigma<-cov(ZZ)/2
        df.general[i+1]<-sum(diag(Sigma))^2/sum(Sigma^2)*(h-1)
        st.general[i+1]<-st.normal[i+1] *sum(diag(Sigma))/sum(Sigma^2)
        pv.general[i+1]<-1-pchisq(st.general[i+1],df.general[i+1])
    }
    z<-data.frame(cbind(st.normal,df.normal,pv.normal,pv.general))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    dimnames(z)<-list(rr,c("Stat","df(Nor)","p.value(Nor)","p.value(Gen)"))
    z
}

# Written by Yongwu Shao, 4/27/2006
dr.coordinate.test.save <-
function (object, hypothesis, d=NULL, chi2approx=object$chi2approx,...) 
{
    gamma <- if (class(hypothesis) == "formula") 
        coord.hyp.basis(object, hypothesis)
        else as.matrix(hypothesis)
    p <- length(object$evalues)
    n <- object$cases
    h <- object$slice.info$nslices
    st <- df <- pv <- 0
    gamma <- (dr.R(object)) %*% gamma
    r <- p - dim(gamma)[2]
    H <- as.matrix(qr.Q(qr(gamma), complete = TRUE)[, (p - r + 
        1):p])
    A <- object$A
    st <- 0
    for (j in 1:h) {
        st <- st + sum((t(H) %*% A[j, , ] %*% H)^2) * n/2
    }
# Normal predictors
    df.normal <- (h - 1) * r * (r + 1)/2
    pv.normal  <- 1 - pchisq(st, df.normal)
# General predictors
    {
        HZ <- dr.z(object) %*% H
        ZZ <- array(0, c(n, r^2))
        for (j in 1:n) {
            ZZ[j, ] <- t(t(HZ[j, ])) %*% t(HZ[j, ])
        }
        wts <- rep(eigen( ((n-1)/n)*cov(ZZ)/2)$values,h-1)
        testg <- dr.pvalue(wts[wts>0],st,a=chi2approx)
    }  
    z <- data.frame(cbind(st, df.normal, pv.normal, testg$pval.adj))
    dimnames(z) <- list("Test", c("Statistic", "df(Nor)", "p.val(Nor)", 
        "p.val(Gen)"))
    z
}

#####################################################################
#     pHd, pHdy and pHdres
#####################################################################

dr.M.phdy <- function(...) {dr.M.phd(...)}
dr.M.mphd <- function(...) stop("Multivariate pHd not implemented!")

dr.M.phdres <- function(...) {dr.M.phd(...)}
dr.M.mphdres <- function(...) stop("Multivariate pHd not implemented!")
dr.M.mphy <- function(...) stop("Multivariate pHd not implemented!")

dr.M.phd <-function(object,...) {
    wts <- dr.wts(object)
    z <- dr.z(object)
    y <- dr.y(object)
    M<- (t(apply(z,2,"*",wts*y)) %*% z) / sum(wts)
    return(list(M=M))
}

dr.M.mphd <- function(...) stop("Multivariate pHd not implemented!")


dr.y.phdy <- function(object) {y <- object$y ; y - mean(y)}
dr.y.phdres <- function(object) { 
   y <- object$y
   sw <- sqrt(object$weights)
   qr.resid(object$qr,sw*(y-mean(y)))}
dr.y.phd <- function(object) {dr.y.phdres(object)}

dr.test.phd<-function(object,numdir=object$numdir,...) {
# Modified by Jorge de la Vega, February, 2001
#compute the phd asymptotic test statitics under restrictions for the
#first numdir directions, based on response = OLS residuals
# order the absolute eigenvalues
    e<-sort(abs(object$evalues))
    p<-length(object$evalues)
# get the response
    resi<-dr.y(object)
    varres<-2*var(resi)
    n<-object$cases
    st<-df<-pv<-0
    nt<-min(p,numdir)
    for (i in 0:(nt-1))
# the statistic is partial sums of squared absolute eigenvalues
      {st[i+1]<-n*(p-i)*mean(e[seq(1,p-i)]^2)/varres
# compute the degrees of freedom
       df[i+1]<-(p-i)*(p-i+1)/2
# use asymptotic chi-square distrbution for p.values.  Requires normality.
       pv[i+1]<-1-pchisq(st[i+1],df[i+1])
      }
# compute the test for complete independence
    indep <- dr.indep.test.phdres(object,st[1])
# compute tests that do not require normal theory (linear combination of 
# chi-squares):
    lc <- dr.test2.phdres(object,st)
# report results
    z<-data.frame(cbind(st,df,pv,c(indep[[2]],rep(NA,length(st)-1)),lc[,2]))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    cc<-c("Stat","df","Normal theory","Indep. test","General theory")
    dimnames(z)<-list(rr,cc)
    z
}

dr.test.phdres<-function(object,numdir,...){dr.test.phd(object,numdir)}

dr.test.phdy<-function(object,numdir,...) {
#compute the phd test statitics for the first numdir directions
#based on response = y.  According to Li (1992), this requires normal
#predictors.
# order the absolute eigenvalues
    e<-sort(abs(object$evalues))
    p<-length(object$evalues)
# get the response
    resi<-dr.y(object)
    varres<-2*var(resi)
    n<-object$cases
    st<-df<-pv<-0
    nt<-min(p,numdir)
    for (i in 0:(nt-1))
# the statistic is partial sums of squared absolute eigenvalues
      {st[i+1]<-n*(p-i)*mean(e[seq(1,p-i)]^2)/varres
# compute the degrees of freedom
       df[i+1]<-(p-i)*(p-i+1)/2
# use asymptotic chi-square distrbution for p.values.  Requires normality.
       pv[i+1]<-1-pchisq(st[i+1],df[i+1])
      }
# report results
    z<-data.frame(cbind(st,df,pv))
    rr<-paste(0:(nt-1),"D vs >= ",1:nt,"D",sep="")
    cc<-c("Stat","df","p.value")
    dimnames(z)<-list(rr,cc)
    z}

#####################################################################
# phdq, Reference:  Li (1992, JASA)
# Function to build quadratic form in the full quadratic fit.
# Corrected phdq by Jorge de la Vega 7/10/01
#####################################################################

dr.M.phdq<-function(object,...){
 pars <- fullquad.fit(object)$coef
 k<-length(pars)
 p<-(-3+sqrt(9+8*(k-1)))/2 #always k=1+2p+p(p-1)/2
 mymatrix <- diag(pars[round((p+2):(2*p+1))]) # doesn't work in S without 'round'
 pars <- pars[-(1:(2*p+1))]
 for (i in 1:(p - 1)) {
      mymatrix[i,(i+1):p] <- pars[1:(p - i)]/2
      mymatrix[(i+1):p,i] <- pars[1:(p - i)]/2
      pars <- pars[-(1:(p - i))]
 }
 return(list(M=mymatrix))
 }

fullquad.fit <-function(object) {
 x <- dr.z(object)
 y <- object$y
 w <- dr.wts(object)
 z <- cbind(x,x^2)
 p <- NCOL(x)
 for (j in 1:(p-1)){
   for (k in (j+1):p) { z <- cbind(z,matrix(x[,j]*x[,k],ncol=1))}}
 lm(y~z,weights=w)
 }

 
dr.y.phdq<-function(object){residuals(fullquad.fit(object),type="pearson")}
dr.test.phdq<-function(object,numdir,...){dr.test.phd(object,numdir)}
dr.M.mphdq <- function(...) stop("Multivariate pHd not implemented!")




#####################################################################
#     Helper methods
#####################################################################
#####################################################################
#  Auxillary function to find matrix H used in marginal coordinate tests
#  We test Y indep X2 | X1, and H is a basis in R^p for the column
#  space of X2.
#####################################################################

coord.hyp.basis <- function(object,spec,which=1) {
 mod2 <- update(object$terms,spec)
 Base <- attr(object$terms,"term.labels")
 New  <- attr(terms(mod2), "term.labels")
 cols <-na.omit(match(New,Base))
 if(length(cols) != length(New)) stop("Error---bad value of 'spec'")
 as.matrix(diag(rep(1,length(Base)))[,which*cols])
 }

#####################################################################
#     Recover the direction vectors
#####################################################################

dr.directions <- function(object, which, x) {UseMethod("dr.direction")}
dr.direction  <- function(object, which, x) {UseMethod("dr.direction")}

# standardize so that the first coordinates are always positive.
dr.direction.default <- 
  function(object, which=NULL,x=dr.x(object)) {  
    ans <- (apply(x,2,function(x){x-mean(x)}) %*% object$evectors)
    which <- if (is.null(which)) seq(dim(ans)[2]) else which 
    ans <- apply(ans,2,function(x) if(x[1] <= 0)  -x else x)[,which]
    if (length(which) > 1) dimnames(ans) <-
                       list ( attr(x,"dimnames")[[1]],paste("Dir",which,sep=""))
    ans
  }

  
#####################################################################
#     Plotting methods
#####################################################################

plot.dr <- function
      (x,which=1:x$numdir,mark.by.y=FALSE,plot.method=pairs,...) {
 d <- dr.direction(x,which)
 if (mark.by.y == FALSE) {
    plot.method(cbind(dr.y(x),d),labels=c(dr.y.name(x),colnames(d)),...)
    }
 else {plot.method(d,labels=colnames(d),col=markby(dr.y(x)),...)}
 }


markby <-
function (z, use = "color", values = NULL, color.fn = rainbow,
    na.action = "na.use")
{
    u <- unique(z)
    lu <- length(u)
    ans <- 0
    vals <- if (use == "color") {
        if (!is.null(values) && length(values) == lu)
            values
        else color.fn(lu)
    }
    else {
        if (!is.null(values) && length(values) == lu)
            values
        else 1:lu
    }
    for (j in 1:lu) if (is.na(u[j])) {
        ans[which(is.na(z))] <- if (na.action == "na.use")
            vals[j]
        else NA
    }
    else {
        ans[z == u[j]] <- vals[j]
    }
    ans
}
   

###################################################################
#
#  basic print method for dimension reduction
#
###################################################################
"print.dr" <-
function(x, digits = max(3, getOption("digits") - 3), width=50, 
    numdir=x$numdir, ...)
{   
    fout <- deparse(x$call,width.cutoff=width)
    for (f in fout) cat("\n",f)
    cat("\n")
    cat("Estimated Basis Vectors for Central Subspace:\n")
    evectors<-x$evectors
    print.default(evectors[,1:min(numdir,dim(evectors)[2])])
    cat("Eigenvalues:\n")
    print.default(x$evalues[1:min(numdir,dim(evectors)[2])])
    cat("\n")
    invisible(x)
}

###################################################################
#  basic summary method for dimension reduction
###################################################################
"summary.dr" <- function (object, ...)
{   z <- object
    ans <- z[c("call")]
    nd <- min(z$numdir,length(which(abs(z$evalues)>1.e-15)))
    ans$evectors <- z$evectors[,1:nd]
    ans$method <- z$method
    ans$nslices <- z$slice.info$nslices
    ans$sizes <- z$slice.info$slice.sizes
    ans$weights <- dr.wts(z)
    sw <- sqrt(ans$weights)
    y <- z$y #dr.y(z)
    ans$n <- z$cases #NROW(z$model)
    ols.fit <- qr.fitted(object$qr,sw*(y-mean(y)))
    angles <- cosangle(dr.direction(object),ols.fit)
    angles <- if (is.matrix(angles)) angles[,1:nd] else angles[1:nd]
    if (is.matrix(angles)) dimnames(angles)[[1]] <- z$y.name
    angle.names <- if (!is.matrix(angles)) "R^2(OLS|dr)" else 
          paste("R^2(OLS for ",dimnames(angles)[[1]],"|dr)",sep="")
    ans$evalues <-rbind (z$evalues[1:nd],angles)
    dimnames(ans$evalues)<-
     list(c("Eigenvalues",angle.names),
          paste("Dir", 1:NCOL(ans$evalues), sep=""))
    ans$test <- dr.test(object,nd)
    class(ans) <- "summary.dr"
    ans
}

###################################################################
#
# basic print.summary method for dimension reduction
#
###################################################################
"print.summary.dr" <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    cat("Method:\n")#S: ' ' instead of '\n'
    if(is.null(x$nslices)){
       cat(paste(x$method, ", n = ", x$n,sep=""))
       if(diff(range(x$weights)) > 0)
                  cat(", using weights.\n") else cat(".\n")}
       else {
         cat(paste(x$method," with ",x$nslices, " slices, n = ",
                   x$n,sep=""))
         if(diff(range(x$weights)) > 0) 
              cat(", using weights.\n") else cat(".\n")
         cat("\nSlice Sizes:\n")#S: ' ' instead of '\n'
         cat(x$sizes,"\n")}
    cat("\nEstimated Basis Vectors for Central Subspace:\n")
    print(x$evectors,digits=digits)
    cat("\n")
    print(x$evalues,digits=digits)
    if (length(x$omitted) > 0){
      cat("\nModel matrix is not full rank.  Deleted columns:\n")
      cat(x$omitted,"\n")}
    if (!is.null(x$test)){
      cat("\nLarge-sample Marginal Dimension Tests:\n")
      #cat("\nAsymp. Chi-square tests for dimension:\n")
      print(as.matrix(x$test),digits=digits)}
    invisible(x)
}



###################################################################3
##
##  Translation of methods in Arc for testing with pHd to R
##  Original lisp functions were mostly written by R. D. Cook
##  Translation to R by S. Weisberg, February, 2001
##
###################################################################3
# this function is a translation from Arc.  It computes the matrices W and
# eW described in Sec. 12.3.1 of Cook (1998), Regression Graphics.
# There are separate versions of this function for R and for Splus because
cov.ew.matrix <- function(object,scaled=FALSE) {
  mat.normalize <- function(a){apply(a,2,function(x){x/(sqrt(sum(x^2)))})}
  n <- dim(dr.x(object))[1]
  TEMPwts <- object$weights
  sTEMPwts <- sqrt(TEMPwts)
  v <- sqrt(n)* mat.normalize(
         apply(scale(dr.x(object),center=TRUE,scale=FALSE),2,"*",sTEMPwts) %*% 
               object$evectors)
  y <- dr.y(object) # get the response
  y <- if (scaled) y-mean(sTEMPwts*y) else 1 # a multiplier in the matrix
  p <- dim(v)[2]
  ew0 <- NULL
  for (i in 1:p){
   for (j in i:p){
    ew0 <- cbind(ew0, if (i ==j) y*(v[,j]^2-1) else y*sqrt(2)*v[,i]*v[,j])}}
  wmean <- function (x,w)  { sum(x * w) / sum (w) }
  tmp <- apply(ew0,2,function(x,w,wmean){sqrt(w)*(x-wmean(x,w))},TEMPwts,wmean)
  ans<-(1/sum(TEMPwts)) * t(tmp) %*% tmp
  ans} 

#translation of :general-pvalues method for phd in Arc
dr.test2.phdres <- function(object,stats){
  covew <- cov.ew.matrix(object,scaled=TRUE)
  C <- .5/var(dr.y(object))
  p <- length(stats)
  pval <- NULL
  d2 <-dim(dr.x(object))[2]
  start <- -d2
  end <- dim(covew)[2]
  for (i in 1:p) {
   start <- start + d2-i+2
   evals <- 
     eigen(as.numeric(C)*covew[start:end,start:end],only.values=TRUE)$values
   pval<-c(pval,
     dr.pvalue(evals,stats[i],chi2approx=object$chi2approx)$pval.adj)
#   pval<-c(pval,wood.pvalue(evals,stats[i]))
   }
# report results
    z<-data.frame(cbind(stats,pval))
    rr<-paste(0:(p-1),"D vs >= ",1:p,"D",sep="")
    dimnames(z)<-list(rr,c("Stat","p.value"))
    z}
   
dr.indep.test.phdres <- function(object,stat) {
  eval <- eigen(cov.ew.matrix(object,scaled=FALSE),only.values=TRUE)
  pval<-dr.pvalue(.5*eval$values,stat,
           chi2approx=object$chi2approx)$pval.adj
# report results
    z<-data.frame(cbind(stat,pval))
    dimnames(z)<-list(c("Test of independence"),c("Stat","p.value"))
    z}

dr.pvalue <- function(coef,f,chi2approx=c("bx","wood"),...){
 method <- match.arg(chi2approx)
 if (method == "bx") {
   bentlerxie.pvalue(coef,f)} else{
   wood.pvalue(coef,f,...)}}


bentlerxie.pvalue <- function(coef,f) {
# Returns the Bentler-Xie approximation to P(coef'X > f), where
# X is a vector of iid Chi-square(1) r.v.'s, coef is a vector of weights,
# and f is the observed value.
  trace1 <- sum(coef)                   
  trace2 <- sum(coef^2)
  df.adj <- trace1^2/trace2   
  stat.adj <- f *  df.adj/trace1
  bxadjusted <- 1 - pchisq(stat.adj, df.adj)  
  ans <- data.frame(test=f,test.adj=stat.adj,
                    df.adj=df.adj,pval.adj=bxadjusted)
  ans
}

wood.pvalue <- function (coef, f, tol=0.0, print=FALSE){
#Returns an approximation to P(coef'X > f) for X=(X1,...,Xk)', a vector of iid
#one df chi-squared rvs.  coef is a list of positive coefficients. tol is used
#to check for near-zero conditions.
#See Wood (1989), Communications in Statistics, Simulation 1439-1456.
#Translated from Arc function wood-pvalue.
#  error checking
  if (min(coef) < 0) stop("negative eigenvalues")
  if (length(coef) == 1)
     {pval <- 1-pchisq(f/coef,1)} else
     {k1 <-     sum(coef)
      k2 <- 2 * sum(coef^2)
      k3 <- 8 * sum(coef^3)
      t1 <- 4*k1*k2^2 + k3*(k2-k1^2)
      t2 <- k1*k3 - 2*k2^2
    if ((t2 <= tol) && (tol < t2) ){
        a1 <- k1^2/k2
    b  <- k1/k2
    pval <- 1 - pgamma(b*f,a1)
        if (print) 
      print(paste("Satterthwaite-Welsh Approximation =", pval))}
      else if( (t1 <= tol) && (tol < t2)){
        a1 <-2 + (k1/k2)^2
    b  <- (k1*(k1^2+k2))/k2
    pval <- if (f < tol) 1.0 else 1 - pgamma(b/f,a1)
        if (print) print(paste("Inverse gamma approximation =",pval))}
      else if ((t1 > tol) && (t2 > tol)) {
        a1 <- (2*k1*(k1*k3 + k2*k1^2 -k2^2))/t1
     b <- t1/t2
    a2 <- 3 + 2*k2*(k2+k1^2)/t2
    pval <- 1-pbeta(f/(f+b),a1,a2)
        if (print) print(paste(
          "Three parameter F(Pearson Type VI) approximation =", pval))}
      else {
        pval <- NULL
        if (print) print("Wood's Approximation failed")}}
   data.frame(test=f,test.adj=NA,df.adj=NA,pval.adj=pval)}


#########################################################################
#
# permutation tests for dimenison reduction
#
#########################################################################

dr.permutation.test <- function(object,npermute=50,numdir=object$numdir) {
 if (inherits(object,"ire")) stop("Permutation tests not implemented for ire")
 else{
   permute.weights <- TRUE
   call <- object$call
   call[[1]] <- as.name("dr.compute")
   call$formula <- call$data <- call$subset <- call$na.action <- NULL
   x <- dr.directions(object) # predictors with a 'better' basis
   call$y <- object$y
   weights <- object$weights
# nd is the number of dimensions to test
   nd <- min(numdir,length(which(abs(object$evalues)>1.e-8))-1)
   nt <- nd + 1
# observed value of the test statistics = obstest
   obstest<-dr.permutation.test.statistic(object,nt)
# count and val keep track of the results and are initialized to zero
   count<-rep(0,nt)
   val<-rep(0,nt)
# main loop
   for (j in 1:npermute){   #repeat npermute times
    perm<-sample(1:object$cases)
    call$weights<- if (permute.weights) weights[perm] else weights    
# inner loop
    for (col in 0:nd){
# z gives a permutation of x.  For a test of dim = col versus dim >= col+1,
# all columns of x are permuted EXCEPT for the first col columns.
        call$x<-if (col==0) x[perm,] else   cbind(x[,(1:col)],x[perm,-(1:col)])
        iperm <- eval(call)
        val[col+1]<- dr.permutation.test.statistic(iperm,col+1)[col+1]
        }  # end of inner loop
# add to counter if the permuted test exceeds the obs. value
     count[val>obstest]<-count[val>obstest]+1
   }
# summarize
   pval<-(count)/(npermute+1)
   ans1 <- data.frame(cbind(obstest,pval))
   dimnames(ans1) <-list(paste(0:(nt-1),"D vs >= ",1:nt,"D",sep=""),
                         c("Stat","p.value"))
   ans<-list(summary=ans1,npermute=npermute)
   class(ans) <- "dr.permutation.test"
   ans
   }}

"print.dr.permutation.test" <-
function(x, digits = max(3, getOption("digits") - 3), ...)
{
   cat("\nPermutation tests\nNumber of permutations:\n")
   print.default(x$npermute)
   cat("\nTest results:\n")
   print(x$summary,digits=digits) 
   invisible(x)
}

"summary.dr.permutation.test" <- function(...)
              {print.dr.permutation.test(...)}


#########################################################################
#
# dr.permutation.test.statistic method
#
#########################################################################

dr.permutation.test.statistic <- function(object,numdir)
  {UseMethod("dr.permutation.test.statistic")}

dr.permutation.test.statistic.default <- function(object,numdir){
   object$cases*rev(cumsum(rev(object$evalues)))[1:numdir]}

dr.permutation.test.statistic.phdy <- function(object,numdir){
       dr.permutation.test.statistic.phd(object,numdir)}
dr.permutation.test.statistic.phdres <- function(object,numdir){
       dr.permutation.test.statistic.phd(object,numdir)}
dr.permutation.test.statistic.phd <- function(object,numdir){
   (.5*object$cases*rev(cumsum(rev(object$evalues^2)))/var(dr.y(object)))[1:numdir]}

#####################################################################
#
#     dr.slices returns non-overlapping slices based on y
#     y is either a list of n numbers or an n by p matrix
#     nslices is either the total number of slices, or else a
#     list of the number of slices in each dimension
#
#####################################################################
dr.slices <- function(y,nslices) {
  dr.slice.1d <- function(y,h) {
      z<-unique(y)
      if (length(z) > h) dr.slice2(y,h) else dr.slice1(y,sort(z))}
  dr.slice1 <- function(y,u){
      z <- sizes <- 0
      for (j in 1:length(u)) {
          temp <- which(y==u[j])
          z[temp] <- j
          sizes[j] <- length(temp) } 
      list(slice.indicator=z, nslices=length(u), slice.sizes=sizes)}
  dr.slice2 <- function(y,h){
       myfind <- function(x,cty){
          ans<-which(x <= cty)
          if (length(ans)==0) length(cty) else ans[1]} 
       or <- order(y)     # y[or] would return ordered y
       cty <- cumsum(table(y))  # cumulative sums of counts of y
       names(cty) <- NULL # drop class names
       n <- length(y)     # length of y
       m<-floor(n/h)      # nominal number of obs per slice
       sp <- end <- 0     # initialize
       j <- 0             # slice counter will end up <= h
       ans <- rep(1,n)    # initialize slice indicator to all slice 1
       while(end < n-2) { # find slice boundaries:  all slices have at least 2 obs
          end <- end+m
          j <- j+1       
          sp[j] <- myfind(end,cty) 
          end <- cty[sp[j]]}
       sp[j] <- length(cty)
       for (j in 2:length(sp)){ # build slice indicator
         firstobs <- cty[sp[j-1]]+1
         lastobs <- cty[sp[j]]
         ans[or[firstobs:lastobs]] <- j}
       list(slice.indicator=ans, nslices=length(sp),
            slice.sizes=c(cty[sp[1]],diff(cty[sp]))) }
  p <- if (is.matrix(y)) dim(y)[2] else 1
  h <- if (length(nslices) == p) nslices else rep(ceiling(nslices^(1/p)),p)
  a <- dr.slice.1d( if(is.matrix(y)) y[,1] else y, h[1])
  if (p > 1){
    for (col in 2:p) {
       ns <- 0
       for (j in unique(a$slice.indicator)) {
         b <- dr.slice.1d(y[a$slice.indicator==j,col],h[col])
         a$slice.indicator[a$slice.indicator==j] <- 
                a$slice.indicator[a$slice.indicator==j] + 10^(col-1)*b$slice.indicator
         ns <- ns + b$nslices}
       a$nslices <- ns }
#recode unique values to 1:nslices and fix up slice sizes
    v <- unique(a$slice.indicator)
    L <- slice.indicator <- NULL
    for (i in 1:length(v)) {
       sel <- a$slice.indicator==v[i]
       slice.indicator[sel] <- i
       L <- c(L,length(a$slice.indicator[sel]))}
    a$slice.indicator <- slice.indicator
    a$slice.sizes <- L }
  a}

dr.slices.arc<-function(y,nslices)  # matches slices produced by Arc
{
  if(is.matrix(y))  stop("dr.slices.arc is used for univariate y only.  Use dr.slices")
  h <- nslices
  or <- order(y)
  n <- length(y)
  m<-floor(n/h)
  r<-n-m*h
  start<-sp<-ans<-0
  j<-1
  while((start+m)<n)
    { if (r==0)
        start<-start
      else 
        {start<-start+1
         r<-r-1
        }
       while (y[or][start+m]==y[or][start+m+1])
          start<-start+1
       sp[j]<-start+m
       start<-sp[j]
       j<-j+1
     }
# next line added 6/17/02 to assure that the last slice has at least 2 obs.
  if (sp[j-1] == n-1) j <- j-1
  sp[j]<-n
  ans[or[1:sp[1]]] <- 1
  for (k in 2:j){ans[ or[(sp[k-1]+1):sp[k] ] ] <- k}
  list(slice.indicator=ans, nslices=j, slice.sizes=c(sp[1],diff(sp)))
}
     
#####################################################################
#
#     Misc. Auxillary functions: cosine of the angle between two 
#     vectors and between a vector and a subspace.
#
#####################################################################

#
# angle between a vector vec and a subspace span(mat)
#
cosangle <- function(mat,vecs){
 ans <-NULL
 if (!is.matrix(vecs)) ans<-cosangle1(mat,vecs) else {
   for (i in 1:dim(vecs)[2]){ans <-rbind(ans,cosangle1(mat,vecs[,i]))}
   dimnames(ans) <- list(colnames(vecs),NULL) }
 ans}

cosangle1 <- function(mat,vec) {
  ans <- 
   cumsum((t(qr.Q(qr(mat))) %*% scale(vec)/sqrt(length(vec)-1))^2) 
# R returns a row vector, but Splus returns a column vector.  The next line
# fixes this difference
  if (version$language == "R") ans else t(ans)
}



#####################################################################
# R Functions for reweighting for elliptical symmetry
# modified from reweight.lsp for Arc
# Sanford Weisberg, sandy@stat.umn.edu
# March, 2001, rev June 2004

# Here is an outline of the function:
#   1.  Estimates of the mean m and covariance matrix S are obtained.  The
#       function cov.rob in the MASS package is used for this purpose.  
#   2.  The matrix X is replaced by Z = (X - 1t(m))S^{-1/2}.  If the columns
#       of X were from a multivariate normal N(m,S), then the rows of Z are
#       multivariate normal N(0, I).
#   3.  Randomly sample from N(0, sigma*I), where sigma is a tuning
#       parameter.  Find the row of Z closest to the random data, and increase
#       its weight by 1
#   4.  Return the weights divided by nsamples and multiplied by n.
#
#     dr.weights
#
#####################################################################
dr.weights <-
function (formula, data = list(), subset, na.action=na.fail, 
   sigma=1,nsamples=NULL,...)
{
# Create design matrix from the formula and call dr.estimate.weights
    mf1 <- match.call(expand.dots=FALSE)
    mf1$... <- NULL # ignore ...
    mf1$covmethod  <- mf1$nsamples <- NULL
    mf1[[1]] <- as.name("model.frame")
    mf <- eval(mf1, sys.frame(sys.parent()))
    mt <- attr(mf,"terms")
    x <- model.matrix(mt, mf)
    int <- match("(Intercept)", dimnames(x)[[2]], nomatch=0)
    if (int > 0) x <- x[, -int, drop=FALSE] # drop the intercept from X
    ans <- cov.rob(x, ...) 
    m <- ans$center
    s<-svd(ans$cov)
    z<-sweep(x,2,m) %*% s$u %*% diag(1/sqrt(s$d))
    n <- dim(z)[1]   # number of obs
    p <- dim(z)[2]   # number of predictors
    ns <- if (is.null(nsamples)) 10*n else nsamples
    dist <- wts <- rep(0,n)  # initialize distances and weights
    for (i in 1:ns) {
       point <- rnorm(p) * sigma      # Random sample from a normal N(0, sigma^2I)
       dist <- apply(z,1,function(x,point){sum((point-x)^2)},point) 
#Dist to each point 
       sel <- dist == min(dist)               # Find closest point(s)
       wts[sel]<-wts[sel]+1/length(wts[sel])} # Increase weights for those points
    w <- n*wts/ns
    if (missing(subset)) return(w)
    if (is.null(subset)) return(w) else {
# find dimension of mf without subset specified
    mf1$subset <- NULL
    w1 <- rep(NA,length=dim(eval(mf1))[1])
    w1[subset] <- w
    return(w1)}}
 
###################################################################
## drop1 methods 
###################################################################
drop1.dr <-
function (object, scope = NULL, update = TRUE, test = "general", 
    trace=1, ...) 
{               
    keep <- if (is.null(scope)) 
        NULL
    else attr(terms(update(object$terms, scope)), "term.labels")
    all <- attr(object$terms, "term.labels")
    candidates <- setdiff(all, keep)
    if (length(candidates) == 0) 
        stop("Error---nothing to drop")
    ans <- NULL
    for (label in candidates) {
        ans <- rbind(ans, dr.coordinate.test(object, as.formula(paste("~.-", 
            label, sep = "")), ...))
    }
    row.names(ans) <- paste("-", candidates)
    ncols <- ncol(ans)
    or <- order(-ans[, if (test == "general") 
        ncols
    else (ncols - 1)])
    form <- formula(object)
    attributes(form) <- NULL
    fout <- deparse(form, width.cutoff = 50)
    if(trace > 0) {
      for (f in fout) cat("\n", f)
      cat("\n")
      print(ans[or, ]) }
    if (is.null(object$stop)) {
        object$stop <- 0
    }
    stopp <- if (ans[or[1], if (test == "general") 
        ncols
    else (ncols - 1)] < object$stop) 
        TRUE
    else FALSE
    if (stopp == TRUE) {
        if(trace > 0) cat("\nStopping Criterion Met\n")
        object$stop <- TRUE
        object
    }
    else if (update == TRUE) {
        update(object, as.formula(paste("~.", row.names(ans)[or][1], 
            sep = "")))
    }
    else invisible(ans)
}
      
dr.step <-
function (object, scope = NULL, d = NULL, minsize = 2, stop = 0, 
     trace=1,...) 
{
    if (is.null(object$stop)) {
        object$stop <- stop
    }
    if (object$stop == TRUE) 
        object
    else {
        minsize <- max(2, minsize)
        keep <- if (is.null(scope)) 
            NULL
        else attr(terms(update(object$terms, scope)), "term.labels")
        all <- attr(object$terms, "term.labels")
        if (length(keep) >= length(all)) {
            if(trace > 0) cat("\nNo more variables to remove\n")
            object
        }
        else if (length(all) <= minsize) {
            if (trace > 0) cat("\nMinimum size reached\n")
            object$numdir <- minsize
            object
        }
        else {
            if (dim(object$x)[2] <= minsize) {
                if (trace > 0) cat("\nMinimum size reached\n")
                object}
            else {   
                obj1 <- drop1(object, scope = scope, d = d, trace=trace, ...)
                dr.step(obj1, scope = scope, d = d, stop = stop,
                       trace=trace, ...)
            }
        }
    }
}

 
#####################################################################
##
##  Add functions to Splus that are built-in to R
##
#####################################################################

if (version$language != "R") {

"is.empty.model" <- function (x)
{
    tt <- terms(x)
    (length(attr(tt, "factors")) == 0) & (attr(tt, "intercept")==0)
}
"NROW" <-
function(x) if(is.array(x)||is.data.frame(x)) nrow(x) else length(x)
"NCOL" <-
function(x) if(is.array(x)||is.data.frame(x)) ncol(x) else as.integer(1)

"colnames" <-
function(x, do.NULL = TRUE, prefix = "col")
{
    dn <- dimnames(x)
    if(!is.null(dn[[2]]))
    dn[[2]]
    else {
    if(do.NULL) NULL else paste(prefix, seq(length=NCOL(x)), sep="")
    }
}

"getOption" <- function(x) options(x)[[1]]
# end of special functions
}

#####################################################################
#     ire:  Cook, RD and Ni, L (2005). Sufficient Dimension reduction via
# inverse regression:  A minimum discrepancy approach.  JASA 410-428
#  This code does all computations by computing the QR factorization of the
#  centered X matrix, and then using sqrt(n) Q in place of X in the
#  computations.  Back transformation to X scale is done only when direction
#  vectors are needed.
#####################################################################
dr.fit.ire <-function(object,numdir=4,nslices=NULL,slice.function=dr.slices,
    tests=TRUE,...){ 
    z <- dr.z(object)    
    h <- if (!is.null(nslices)) nslices else max(8, NCOL(z)+3)
    slices <- slice.function(dr.y(object),h)
    object$slice.info <- slices
    f <- slices$slice.sizes/sum(slices$slice.sizes)
    n <- object$cases
    weights <- object$weights
    p <- dim(z)[2]
    numdir <- min(numdir,p-1)               
    h <- slices$nslices
    xi <- matrix(0,nrow=p,ncol=slices$nslices)
    for (j in 1:h){
     sel <- slices$slice.indicator==j
     xi[,j]<-apply(z[sel,],2,function(a,w) sum(a*w)/sum(w),weights[sel])
    } 
    object$sir.raw.evectors <- eigen(xi %*% diag(f) %*% t(xi))$vectors
    object$slice.means <- xi # matrix of slice means
    An <- qr.Q(qr(contr.helmert(slices$nslices)))
    object$zeta <- xi %*% diag(f) %*% An
    rownames(object$zeta) <- paste("Q",1:dim(z)[2],sep="")
    ans <- NULL

    if (tests == TRUE) {    
      object$indep.test <- dr.test(object,numdir=0,...)
      Gz <- Gzcomp(object,numdir)  # This is the same for all numdir > 0
      for (d in 1:numdir){
         ans[[d]] <- dr.test(object,numdir=d,Gz,...)
         colnames(ans[[d]]$B) <- paste("Dir",1:d,sep="")     
      }
# This is different from Ni's lsp code.  Gamma_zeta is computed from the
# fit of the largest dimension and stored.
    object$Gz <- Gzcomp(object,d,span=ans[[numdir]]$B)
    }
    aa<-c(object, list(result=ans,numdir=numdir,n=n,f=f))
    class(aa) <- class(object)
    return(aa)
}
     
dr.test.ire <- function(object,numdir,Gz=Gzcomp(object,numdir),steps=1,...){       
   ans <- dr.iteration(object,Gz,d=numdir,
             B=object$sir.raw.evectors[,1:numdir,drop=FALSE],...)
   if (steps > 0 & numdir > 0){ # Reestimate if steps > 0.
    for (st in 1:steps){
     Gz <- Gzcomp(object,numdir,ans$B[,1:numdir])
     ans <- dr.iteration(object,Gz,d=numdir,
             B=ans$B[,1:numdir,drop=FALSE],...)}}
# reorder B matrix according to importance (col. 2 of p. 414)
   if (numdir > 1){
     ans0 <- dr.iteration(object,Gz,d=1,T=ans$B)
     sumry <- ans0$summary
     B0 <- matrix(ans0$B/sqrt(sum((ans0$B)^2)),ncol=1)
     C <- ans$B
     for (d in 2:numdir) {
        common <- B0[,d-1]
        for (j in 1:dim(C)[2]){
            C[,j] <- C[,j] - sum(common*(C[,j]))*B0[,d-1]}
        C <- qr.Q(qr(C))[,-dim(C)[2],drop=FALSE]
        ans0 <- dr.iteration(object,Gz,d=1,T=C)
        B0 <- cbind(B0,ans0$B/sqrt(sum((ans0$B)^2)))
        sumry <- rbind(sumry,dr.iteration(object,Gz,d=d,T=B0)$summary)
        }
# scale to length 1 and make first element always positive
    ans$B <- apply(B0,2,function(x){
        b <- x/sqrt(sum(x^2))
        if (b[1] < 0) - b else b})
    ans$sequential <- sumry
    }
   if (numdir > 0) {colnames(ans$B) <- paste("Dir",1:numdir,sep="")}
   if (numdir > 1) rownames(ans$sequential) <- paste(1:numdir)
   ans}
   
Gzcomp <- function(object,numdir,span){UseMethod("Gzcomp")}
Gzcomp.ire <- function(object,numdir,span=NULL){
    slices <- object$slice.info
    An <- qr.Q(qr(contr.helmert(slices$nslices)))
    n <- object$cases
    weights <- object$weights
    z <- dr.z(object)  # z is centered with correct weights.
    p <- dim(object$slice.means)[1]
    h <- slices$nslices
    f <- slices$slice.sizes/sum(slices$slice.sizes)
# page 411, eq (1) except projecting to d dimensions using span.
    span <- if (!is.null(span)) span else diag(rep(1,p))
    xi <- if (numdir > 0){qr.fitted(qr(span),object$slice.means)} else {
          matrix(0,nrow=p,ncol=h)}        
# We next compute the inner product matrix Vn = inverse(Gamma_zeta)
# We compute only the Cholesky decomposition of Gamma_zeta, as that
# is all that is needed.  The name is reused for intermediate quantities
# to save space (it is a p*(h-1) by p*(h-1) matrix).
# First, compute Gamma, defined on line 2 p. 414
# make use of vec(A %*% B) = kronecker(t(B), A)
# and also that the mean of vec, as defined below, is zero.
    Gz <- array(0,c(p*h,p*h)) 
    Gmat <- NULL
    for (i in 1:n){
     epsilon_i <- -f -z[i,,drop=FALSE] %*% xi %*% diag(f)
     epsilon_i[slices$slice.indicator[i]]<-1+epsilon_i[slices$slice.indicator[i]]
     vec <- as.vector( z[i,] %*% epsilon_i)
     Gmat <- rbind(Gmat,vec)}
     Gz <- chol(kronecker(t(An),diag(rep(1,p))) %*% (((n-1)/n)*cov(Gmat)) %*%
                 kronecker(An,diag(rep(1,p))))
     Gz}
   
#####################################################################
##  ire iteration
##  B, if set, is a p by d matrix whose columns are a starting value
##     for a basis for the central subspace.  The default is the first
##     d columns of I_p
##  R, if set, is a restriction matrix, such that the estimated
##     spanning vectors are RB, not B itself.
##  Returned matrix B is really R %*% B
##  --fn is the objective function, eq (5) of Cook & Ni (2004)
##  --updateC is step 2 of the algorithm on p. 414
##  --updateB is step 3 of the algorithm on p. 414
##  --PBk in the paper is the projection on an orthogonal completment
##      of the operator QBk.  It is the QR decomposition of the projection,
##      not its complement (hence the use of qr.resid, not qr.fitted)
######################################################################
dr.iteration <- function(object,Gz,d=2,B,T,eps,itmax,verbose){UseMethod("dr.iteration")}
dr.iteration.ire <- function(object,Gz,d=2,B=NULL,T=NULL,eps=1.e-6,itmax=200,
   verbose=FALSE){
   n <- object$cases
   zeta <- object$zeta
   p <- dim(zeta)[1]
   h1 <- dim(zeta)[2]  # zeta is p by (h-1) for ire and sum(h-1) for pire
   if (d == 0){ 
     err <- n*sum(forwardsolve(t(Gz),as.vector(zeta))^2)
     data.frame(Test=err,df=(p-d)*(h1-d),
               p.value=pchisq(err,(p-d)*(h1-d),lower.tail=FALSE),iter=0)} else {
   T <- if(is.null(T)) diag(rep(1,p)) else T
   B <- if(is.null(B)) diag(rep(1,ncol(T)))[,1:d,drop=FALSE] else B
   fn <- function(B,C){ 
           n * sum( forwardsolve(t(Gz),as.vector(zeta)-as.vector(T%*%B%*%C))^2 ) }
   updateC <- function() {
        matrix( qr.coef(qr(forwardsolve(t(Gz),kronecker(diag(rep(1,h1)),T%*%B))),
                           forwardsolve(t(Gz),as.vector(zeta))), nrow=d)}
   updateB <- function() { 
      for (k in 1:d) { 
         alphak <- as.vector(zeta - T %*% B[,-k,drop=FALSE] %*% C[-k,])
         PBk <- qr(B[,-k])  
         bk <- qr.coef(
                qr(forwardsolve(t(Gz),
                    t(qr.resid(PBk,t(kronecker(C[k,],T)))))),
                forwardsolve(t(Gz),as.vector(alphak)))
         bk[is.na(bk)] <- 0  # can use any OLS estimate; eg, set NA to 0
         bk <- qr.resid(PBk,bk)
         B[,k] <- bk/sqrt(sum(bk^2))}
         B}
     C <- updateC()  # starting values
     err <- fn(B,C)
     iter <- 0
     repeat{
       iter <- iter+1
       B <- updateB()
       C <- updateC()
       errold <- err
       err <- fn(B,C)
       if(verbose==TRUE) print(paste("Iter =",iter,"Fn =",err),quote=FALSE)
       if ( abs(err-errold)/errold < eps || iter > itmax ) break
       }  
     B <- T %*% B
     rownames(B) <- rownames(zeta)
     list(B=B,summary=data.frame(Test=err,df=(p-d)*(h1-d),
               p.value=pchisq(err,(p-d)*(h1-d),lower.tail=FALSE),iter=iter))
     }}

# Equation (12), p. 415 of Cook and Ni (2004) for d=NULL
# Equation (14), p. 415 for d > 0
# Remember...all calculations are in the Q scale, where X=QR...
dr.coordinate.test.ire<-function(object,hypothesis,d=NULL,...){
    gamma <- if (class(hypothesis) == "formula")
        coord.hyp.basis(object, hypothesis)
        else as.matrix(hypothesis)
    gamma <- dr.R(object)%*%gamma  # Rotate to Q-coordinates:
    p <- dim(object$x)[2]
    r <- p-dim(gamma)[2] 
    maxdir <- length(object$result)   
   if(is.null(d)){
    h1 <- dim(object$zeta)[2] # h-1 for ire, sum(h-1) for pire
    H <- qr.Q(qr(gamma),complete=TRUE)[,-(1:(p-r)),drop=FALSE]
    n<-object$cases
    Gz <- object$Gz 
    zeta <- object$zeta
    m1 <- Gz %*% kronecker(diag(rep(1,h1)),H)
    m1 <- chol(t(m1) %*% m1)
    T_H <- n * sum (forwardsolve(t(m1),as.vector(t(H)%*%zeta))^2)
    df <- r*(h1)
    z <- data.frame(Test=T_H,df=df,p.value=pchisq(T_H,df,lower.tail=FALSE))
    z} 
   else {
    F0 <-if(maxdir >= d) object$result[[d]]$summary$Test else
          dr.iteration(object,object$Gz,d=d,...)$summary$Test
    F1 <- dr.joint.test(object,hypothesis,d=d,...)$summary$Test
    data.frame(Test=F1-F0,df=r*d,p.value=pchisq(F1-F0,r*d,lower.tail=FALSE))
    }}
 
# Unnumbered equation middle of second column, p. 415 of Cook and Ni (2004)   
dr.joint.test.ire<-function(object,hypothesis,d=NULL,...){
    if(is.null(d)) {dr.coordinate.test(object,hypothesis,...)} else {
     gamma <- if (class(hypothesis) == "formula")
         coord.hyp.basis(object, hypothesis)
         else as.matrix(hypothesis)
     gamma <- dr.R(object)%*%gamma  # Rotate to Q-coordinates:         
     dr.iteration(object,object$Gz,d=d,T=gamma)}}
    
### print/summary functions
print.ire <- function(x, width=50, ...) { 
    fout <- deparse(x$call,width.cutoff=width)
    for (f in fout) cat("\n",f)
    cat("\n")
    numdir <- length(x$result)
    tests <- x$indep.test
    for (d in 1:numdir) {
     tests <- rbind(tests,x$result[[d]]$summary)}   
    rownames(tests) <- paste(0:numdir,"D vs"," > ",0:numdir,"D",sep="")
    cat("Large-sample Marginal Dimension Tests:\n")
    print(tests)
    cat("\n")
    invisible(x)
}

"summary.ire" <- function (object, ...)
{   ans <- object[c("call")]
    result <- object$result
    numdir <- length(result)
    tests <- object$indep.test
    for (d in 1:numdir) {
       tests <- rbind(tests,result[[d]]$summary)}
    rownames(tests) <- paste(0:numdir,"D vs"," > ",0:numdir,"D",sep="")
    ans$method <- object$method
    ans$nslices <- object$slice.info$nslices
    ans$sizes <- object$slice.info$slice.sizes
    ans$weights <- dr.wts(object)
    ans$result <- object$result
    for (j in 1:length(ans$result)) {ans$result[[j]]$B <- dr.basis(object,j)}
    ans$n <- object$cases #NROW(object$model)
    ans$test <- tests
    class(ans) <- "summary.ire"
    ans
}

"print.summary.ire" <-
function (x, digits = max(3, getOption("digits") - 3), ...)
{
    cat("\nCall:\n")#S: ' ' instead of '\n'
    cat(paste(deparse(x$call), sep="\n", collapse = "\n"), "\n\n", sep="")
    cat("Method:\n")#S: ' ' instead of '\n'
    cat(x$method,"with",x$nslices, " slices, n =",x$n) 
    if(diff(range(x$weights)) > 0)cat(", using weights.\n") else cat(".\n")
    cat("\nSlice Sizes:\n")#S: ' ' instead of '\n'
    cat(x$sizes,"\n")
    cat("\nLarge-sample Marginal Dimension Tests:\n")
      print(as.matrix(x$test),digits=digits)
    cat("\n")       
    cat("\nSolutions for various dimensions:\n")
    for (d in 1:length(x$result)){
     cat("\n",paste("Dimension =",d),"\n")
     print(x$result[[d]]$B,digits=digits)}
     cat("\n")
     invisible(x)
}

dr.basis.ire <- function(object,numdir=length(object$result)) {
    fl <- function(z) apply(z,2,function(x){
        b <- x/sqrt(sum(x^2))
        if (b[1] < 0) -1*b else b}) 
    ans <- fl(backsolve(dr.R(object),object$result[[numdir]]$B))
    dimnames(ans) <- list(colnames(dr.x(object)),
                         paste("Dir",1:dim(ans)[2],sep=""))
    ans}

dr.direction.ire <- 
  function(object, which=1:length(object$result),x=dr.x(object)){
    d <- max(which)
    scale(x,center=TRUE,scale=FALSE) %*% dr.basis(object,d)
    }
   
#############################################################################
# partial ire --- see Wen and Cook (in press), Optimal sufficient dimension
# reduction in regressions with categorical predictors. Journal of Statistical
# planning and inference.
#############################################################################

dr.fit.pire <-function(object,numdir=4,nslices=NULL,slice.function=dr.slices,...){ 
    y <- dr.y(object)
    z <- dr.z(object) 
    p <- dim(z)[2]
    if(is.null(object$group)) object$group <- rep(1,dim(z)[1])
    group.names <- unique(as.factor(object$group))
    nw <- table(object$group)
    if (any(nw < p) ) stop("At least one group has too few cases")
    h <- if (!is.null(nslices)) nslices else NCOL(z)
    group.stats <- NULL
    for (j in 1:length(group.names)){
      name <- group.names[j] 
      group.sel <- object$group==name
      ans <- dr.compute(z[group.sel,],y[group.sel],method="ire",
            nslices=h,slice.function=slice.function,
            tests=FALSE,weights=object$weights[group.sel])
      object$zeta[[j]] <- ans$zeta
      group.stats[[j]] <- ans
      }
    object$group.stats <- group.stats
    numdir <- min(numdir,p-1)
    object$sir.raw.evectors <- dr.compute(z,y,nslices=h,
            slice.function=slice.function,
            weights=object$weights)$raw.evectors
    class(object) <- c("pire","ire","dr")
    object$indep.test <- dr.test(object,numdir=0,...)
    Gz <- Gzcomp(object,numdir)  # This is the same for all numdir > 0
    ans <- NULL
    for (d in 1:numdir){
       ans[[d]] <- dr.test(object,numdir=d,Gz,...)
       colnames(ans[[d]]$B) <- paste("Dir",1:d,sep="")     
    }
    object$Gz <- Gzcomp(object,d,span=ans[[numdir]]$B)
    aa<-c(object, list(result=ans,numdir=numdir))
    class(aa) <- class(object)
    return(aa)
}

dr.iteration.pire <- function(object,Gz,d=2,B=NULL,T=NULL,eps=1.e-6,itmax=200,
   verbose=FALSE){
   gsolve <- function(a1,a2){  #modelled after ginv in MASS
    Asvd <- svd(a1)
    Positive <- Asvd$d > max(sqrt(.Machine$double.eps) * Asvd$d[1], 0)
    if(all(Positive))
     Asvd$v %*% (1/Asvd$d * t(Asvd$u)) %*% a2
    else Asvd$v[, Positive, drop = FALSE] %*% ((1/Asvd$d[Positive]) * 
        t(Asvd$u[, Positive, drop = FALSE])) %*% a2}
   n <- object$cases
   zeta <- object$zeta
   n.groups <- length(zeta)
   p <- dim(zeta[[1]])[1]
   h1 <- 0
   h2 <- NULL
   for (j in 1:n.groups){
      h2[j] <- dim(zeta[[j]])[2]
      h1 <- h1 + h2[j]}
   if (d == 0){ 
     err <- 0
     for (j in 1:n.groups){
      err <- err + n*sum(forwardsolve(t(Gz[[j]]),as.vector(zeta[[j]]))^2)}
     data.frame(Test=err,df=(p-d)*(h1-d),
               p.value=pchisq(err,(p-d)*(h1-d),lower.tail=FALSE),iter=0)} else {
   T <- if(is.null(T)) diag(rep(1,p)) else T
   B <- if(is.null(B)) diag(rep(1,ncol(T)))[,1:d,drop=FALSE] else B
   fn <- function(B,C){ 
           ans <- 0
           for (j in 1:n.groups){ans <- ans +
            n * sum( forwardsolve(t(Gz[[j]]),as.vector(zeta[[j]])-
                as.vector(T%*%B%*%C[[j]]))^2) }
            ans}
   updateC <- function() {
        C <- NULL
        for (j in 1:n.groups){ C[[j]]<-
        matrix( qr.coef(qr(forwardsolve(t(Gz[[j]]),kronecker(diag(rep(1,h2[j])),T%*%B))),
                           forwardsolve(t(Gz[[j]]),as.vector(zeta[[j]]))), nrow=d)}
        C}
   updateB <- function() { 
      for (k in 1:d) { 
         PBk <- qr(B[,-k]) 
         a1 <- a2 <- 0
         for (j in 1:n.groups){
          alphak <- as.vector(zeta[[j]]-T%*%B[,-k,drop=FALSE]%*%C[[j]][-k,])
          m1 <-  forwardsolve(t(Gz[[j]]),
                    t(qr.resid(PBk,t(kronecker(C[[j]][k,],T)))))
          m2 <- forwardsolve(t(Gz[[j]]),alphak)
          a1 <- a1 + t(m1) %*% m1
          a2 <- a2 + t(m1) %*% m2}
         bk <- qr.resid(PBk, gsolve(a1,a2))
         bk[is.na(bk)] <- 0  # can use any OLS estimate; eg, set NA to 0
         bk <- qr.resid(PBk,bk)
         B[,k] <- bk/sqrt(sum(bk^2))}
         B}
     C <- updateC()  # starting values
     err <- fn(B,C)
     iter <- 0
     repeat{
       iter <- iter+1
       B <- updateB()
       C <- updateC()
       errold <- err
       err <- fn(B,C)
       if(verbose==TRUE) print(paste("Iter =",iter,"Fn =",err),quote=FALSE)
       if ( abs(err-errold)/errold < eps || iter > itmax ) break
       }  
     B <- T %*% B
     rownames(B) <- rownames(zeta)
     list(B=B,summary=data.frame(Test=err,df=(p-d)*(h1-d),
               p.value=pchisq(err,(p-d)*(h1-d),lower.tail=FALSE),iter=iter))
     }}
     
dr.coordinate.test.pire<-function(object,hypothesis,d=NULL,...){
    gamma <- if (class(hypothesis) == "formula")
        coord.hyp.basis(object, hypothesis)
        else as.matrix(hypothesis)
    gamma <- dr.R(object)%*%gamma  # Rotate to Q-coordinates:
    p <- object$qr$rank
    r <- p-dim(gamma)[2] 
    maxdir <- length(object$result) 
    n.groups <- length(object$group.stats)
    h1 <- 0
    h2 <- NULL
    zeta <- object$zeta
    for (j in 1:n.groups){
      h2[j] <- dim(zeta[[j]])[2]
      h1 <- h1 + h2[j]}  
   if(is.null(d)){
    H <- qr.Q(qr(gamma),complete=TRUE)[,-(1:(p-r)),drop=FALSE]
    n<-object$cases
    Gz <- object$Gz 
    T_H <- 0
    for (j in 1:n.groups){
      m1 <- Gz[[j]] %*% kronecker(diag(rep(1,h2[j])),H)
      m1 <- chol(t(m1) %*% m1)
     T_H <- T_H + n * sum (forwardsolve(t(m1),as.vector(t(H)%*%zeta[[j]]))^2)}
    df <- r*(h1)
    z <- data.frame(Test=T_H,df=df,p.value=pchisq(T_H,df,lower.tail=FALSE))
    z} 
   else {
    F0 <-if(maxdir >= d) object$result[[d]]$summary$Test else
          dr.iteration(object,object$Gz,d=d,...)$summary$Test
    F1 <- dr.joint.test(object,hypothesis,d=d,...)$summary$Test
    data.frame(Test=F1-F0,df=r*d,p.value=pchisq(F1-F0,r*d,lower.tail=FALSE))
    }}
    
Gzcomp.pire <- function(object,numdir,span=NULL){
    Gz <- NULL
    n.groups <- length(object$group.stats)
    pw <- sum(object$group.stats[[1]]$weights)/sum(object$weights)
    Gz[[1]] <- Gzcomp(object$group.stats[[1]],numdir=numdir,span=span)/sqrt(pw)
    if (n.groups > 1){
      for (j in 2:n.groups){
       pw <- sum(object$group.stats[[j]]$weights)/sum(object$weights)
       Gz[[j]] <- 
         Gzcomp(object$group.stats[[j]],numdir=numdir,span=span)/sqrt(pw)}}
    Gz
    }
    
"summary.pire" <- function (object, ...)
{   ans <- object[c("call")]
    result <- object$result
    numdir <- length(result)
    tests <- object$indep.test
    for (d in 1:numdir) {
       tests <- rbind(tests,result[[d]]$summary)}
    rownames(tests) <- paste(0:numdir,"D vs"," > ",0:numdir,"D",sep="")
    ans$method <- object$method
    ans$nslices <- ans.sizes <- NULL
    ans$n <-NULL
    for (stats in object$group.stats){
     ans$n <- c(ans$n,stats$n)
     ans$nslices <- c(ans$nslices,stats$slice.info$nslices)
     ans$sizes <- c(ans$sizes,stats$slice.info$slice.sizes)
     }
    ans$result <- object$result
    for (j in 1:length(ans$result)) {ans$result[[j]]$B <- dr.basis(object,j)}
    ans$weights <- dr.wts(object)
    ans$test <- tests
    class(ans) <- "summary.ire"
    ans
}
    
