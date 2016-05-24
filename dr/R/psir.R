#############################################################################
# partial sir --- see Chiaromonte, Cook and Li (2002).  The implemetation
# here follows Shao, Y. (2007) Topics on dimension reduction, unpublished PhD
# Dissertation, School of Statistics, University of Minnesota.
#############################################################################
dr.fit.psir <-function(object,numdir=4,nslices=2,pool=FALSE,
   slice.function=dr.slices,...){ 
    object <- psir(object,nslices,pool,slice.function)
    object$numdir <- numdir
    object$method <- "psir"
    class(object) <- c("psir", "sir","dr")
    return(object)
}

# The function psir() is the core function which does the actual fitting.
# The function takes a dr object as input.  The outputs are:
# 1. evectors and evalues: the estimated partial central subspace and
#    corresponding eigenvalues.
# 2. test(): the function to test marginal dimensional hypothesis.
# 3. coordinate.test(): the function to test coordinate hypothesis.
psir <- function(object,nslices,pool,slice.function) {
    Y <- dr.y(object)
    X <- dr.x(object)
    W <- dr.wts(object)
    n <- dim(X)[1]
    p <- dim(X)[2]
    group.names <- unique(as.factor(object$group))
    nG <- length(group.names)
    G <- numeric(n)
    for (j in 1:nG) G[object$group ==group.names[j]] <- j
    "%^%"<-function(A,n) { 
        if (dim(A)[2]==1) { A^n } else {
            eg<-eigen(A)
            (eg$vectors) %*% diag(abs(eg$values)^n) %*% t(eg$vectors) }}
    Sigma.pool <- matrix(0,p,p)
    Sigma <- array(0,c(p,p,nG))
    wt.cov <- function(x,w){ 
      xbarw <- apply(x,2,function(x) sum(w*x)/sum(w))
      xc <- t(apply(x,1,function(x) x-xbarw))
      (1/sum(w)) * t(xc) %*% apply(xc,2,function(x) w*x)
      }
    slice <- if (length(nslices)==nG) nslices else rep(nslices,nG)
    info <- NULL
    for (k in 1:nG){
      sel <- object$group == group.names[k]
      info[[k]]<- slice.function(Y[sel],slice[k])
      Sigma[,,k] <- wt.cov(X[sel,],W[sel])
      Sigma.pool <- Sigma.pool + sum(W[sel]) *Sigma[,,k]/ n
    }
    slice.info <- function() info
    # psir1 looks at just one level of group and returns the array of
    # slice means, and a function that when called computes a test stat.
    psir1 <- function(k) {
        group.sel <- object$group == group.names[k]
        Z <- X[group.sel,]
        y <- Y[group.sel]
        n <- length(y)
        Scale <- if(pool) Sigma.pool else Sigma[,,k]
        z <- (Z-matrix(1,n,1)%*%apply(Z,2,mean)) %*% (Scale %^% (-1/2))  
        zmeans <- array(0,c(p,info[[k]]$nslices))
        for (j in 1:info[[k]]$nslices) {
            slice.sel <- (info[[k]]$slice.indicator==j)
            if (sum(slice.sel)<1) mu <- rep(0,p)
            else mu <- apply(z[slice.sel,],2,mean)
            zmeans[,j] <- sqrt(sum(slice.sel)/n)*mu
        }
        # The function stat() computes the coorindate test statistic for the group,
        stat <- function(gamma) {
            r <- p - dim(gamma)[2]
            gamma <- Scale %^% (1/2) %*% gamma
            h <- info[[k]]$nslices
            H <- as.matrix(qr.Q(qr(gamma), complete = TRUE)[, (p - r + 1):p])
            st <- sum((t(H)%*%zmeans)^2)*n
            wts <- rep(1,min(p,h-1))
            wts[1:min(p,h-1)] <- 1-svd(zmeans)$d[1:min(p,h-1)]^2
            return(list(st=st,wts=wts))
        }
        return(list(zmeans=zmeans,stat=stat))
    }
# The following code computes the zmeans matrix.
    zmeans <- NULL
    for (k in 1:nG){zmeans <- cbind(zmeans,psir1(k)$zmeans)}
# The following code estimates the partial central mean subspace.
    D <- svd(zmeans,p)
    evectors <- 
     apply(Sigma.pool%^%(-1/2)%*%D$u,2,function(x)x/sqrt(sum(x^2)))
    
    dimnames(evectors) <- 
         list(attr(dr.x(object),"dimnames")[[2]],
              paste("Dir",1:dim(evectors)[2],sep=""))
# The function test() does the sequential marginal dimension testing 
# from 0 to d.  Function is public.
    test <- function(d) {
        h <- sum(slice)
        d <- min(d,h-nG)
        st <- df <- pv <- 0
        for (i in 0:(d-1)) {
            st[i+1] <- sum(D$d[(i+1):min(p,h-nG)]^2)*n
            df[i+1] <- (h-i-nG)*(p-i) 
            pv[i+1] <- 1 - pchisq(st[i+1], df[i+1])}
        z<-data.frame(cbind(st,df,pv))
        rr<-paste(0:(d-1),"D vs >= ",1:d,"D",sep="")
        dimnames(z)<-list(rr,c("Stat","df","p.value")) 
        return(z)
    }
# The function coordinate.test() tests the coordinate hypothesis.
# The test statistic is constructed by summing up the test statistics in each group.
# The weights are contructed by combining all the weights in each group.
# Function is public.
    coordinate.test <- function(H) {
        r <- p-dim(H)[2]
        st <- 0
        wts <- 0
        for (k in 1:nG) {
            tp <- psir1(k)$stat(H)
            st <- st+tp$st
            wts <- cbind(wts,tp$wts)
        }
        wts <- rep(wts,r)
        testr <- dr.pvalue(wts[wts>1e-5],st,a=object$chi2approx)
        df <- testr$df.adj
        pv <- testr$pval.adj
        return(data.frame(cbind(Test=st,P.value=pv)))
    }
    return(c(object, list(evectors=evectors,
       evalues=D$d^2,slice.info=slice.info,
       test=test,coordinate.test=coordinate.test)))
}

dr.test.psir <- function(object,numdir=object$numdir,...) 
  object$test(numdir)

dr.coordinate.test.psir <- function(object,hypothesis,d=NULL,...) {
    gamma <- if (class(hypothesis) == "formula")
        coord.hyp.basis(object, hypothesis)
        else as.matrix(hypothesis)
    object$coordinate.test(gamma)
}
    
summary.psir <- function(object,...) {
 ans <- summary.dr(object,...)
 ans$method <- "psir"
 gps<- sizes <- NULL
 for (g in 1:length(a1 <- object$slice.info())){
   gps<- c(gps,a1[[g]][[2]])
   sizes <- c(sizes,a1[[g]][[3]])}
 ans$nslices <- paste(gps,collapse=" ")
 ans$sizes <- sizes
 return(ans)
 }
 
summary.psir <- function (object, ...)
{   z <- object
    ans <- z[c("call")]
    nd <- min(z$numdir,length(which(abs(z$evalues)>1.e-15)))
    ans$evectors <- z$evectors[,1:nd]
    ans$method <- "psir"
    gps<- sizes <- NULL
    for (g in 1:length(a1 <- object$slice.info())){
      gps<- c(gps,a1[[g]][[2]])
      sizes <- c(sizes,a1[[g]][[3]])}
    ans$nslices <- paste(gps,collapse=" ")
    ans$sizes <- sizes
    ans$weights <- z$weights
    sw <- sqrt(ans$weights)
    y <- z$y #dr.y(z)
    ans$n <- z$cases #NROW(z$model)
    ols.fit <- qr.fitted(object$qr,sw*(y-mean(y)))
    angles <- cosangle(dr.direction(object),ols.fit)
    angles <- if (is.matrix(angles)) angles[,1:nd] else angles[1:nd]
    if (is.matrix(angles)) dimnames(angles)[[1]] <- z$y.name
    angle.names <- if (!is.matrix(angles)) "R^2(OLS|dr)" else
                        paste("R^2(",dimnames(angles)[[1]],"|dr)",sep="")
    ans$evalues <-rbind (z$evalues[1:nd],angles)
    dimnames(ans$evalues)<-
     list(c("Eigenvalues",angle.names),
          paste("Dir", 1:NCOL(ans$evalues), sep=""))
    ans$test <- dr.test(object,nd)
    class(ans) <- "summary.dr"
    ans
}
 
drop1.psir <- function(object,scope=NULL,...){
    drop1.dr(object,scope,...) }
