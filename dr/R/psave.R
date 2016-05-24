#############################################################################
# partial save --- Shao, Cook and Weisberg (submitted).  The implemetation
# here follows Shao, Y. (2007) Topics on dimension reduction, unpublished PhD
# Dissertation, School of Statistics, University of Minnesota.
#############################################################################

dr.fit.psave <-function(object,numdir=4,nslices=2,pool=FALSE,
   slice.function=dr.slices,...){
    object <- psave(object,nslices,pool,slice.function)
    object$result <- object$evectors[,1:numdir]
    object$numdir <- numdir
    object$method <- "psave"
    class(object) <- c("psave", "save", "dr")
    return(object)
}

# The function psave() is the core function which does the actual fitting and testing.
# dr.fit.psave() and dr.coordinate.test.save() are just wrappers.
# The function takes a data set as input.  The outputs are:
# 1. directions: the estimated partial central subspace.
# 2. test(): the function to test marginal dimensional hypothesis.
# 3. coordinate.test(): the function to test marginal coordinate hypothesis.
psave <- function(object,nslices,pool,slice.function) {
    Y <- dr.y(object)
    X <- dr.x(object)
    W <- dr.wts(object)
    n <- dim(X)[1]
    p <- dim(X)[2]
    group.names <- unique(as.factor(object$group))
    nG <- length(group.names)
    G <- numeric(n)
    for (j in 1:nG) G[object$group ==group.names[j]] <- j
    wt.cov <- function(x, w) {
        (if (sum(w) > 1) {(length(w) - 1)/sum(w)} else {0}) * 
                 var(sweep(x, 1, sqrt(w), "*"))}
    slice <- if (length(nslices)==nG) nslices else rep(nslices,nG)
    info <- NULL
    M <- matrix(0,p,p)
    Sigma.pool <- matrix(0,p,p)
    Sigma <- array(0,c(p,p,nG))
    "%^%"<-function(A,n) { 
        if (dim(A)[2]==1) { A^n } else {
            eg<-eigen(A)
            (eg$vectors) %*% diag(abs(eg$values)^n) %*% t(eg$vectors) }}
    info <- numslices <- NULL
    for (k in 1:nG){
      sel <- object$group == group.names[k]
      info[[k]]<- slice.function(Y[sel],slice[k])
      numslices <- c(numslices, info[[k]]$nslices)
      Sigma[,,k] <- wt.cov(X[sel,],W[sel])
      Sigma.pool <- Sigma.pool + sum(W[sel]) *Sigma[,,k]/ n
    }
    A <- array(0,c(sum(numslices),p,p))
    slice.info <- function() info
# One group at a time
    psave1 <- function(z,y,w,slices) {
        h <- slices$nslices
        n <- length(y)
        M <- matrix(0,p,p)
        A <- array(0,c(h,p,p))
        for (j in 1:h) {
            slice.ind <- (slices$slice.indicator==j)
            IminusC <- if (sum(slice.ind)<3)  diag(rep(1,p))
               else diag(rep(1,p))-wt.cov(z[slice.ind,],w[slice.ind]) 
            M <- M+sum(slice.ind)/n* IminusC%*%IminusC
            A[j,,] <- sqrt(sum(w[slice.ind])/sum(w))*IminusC
        }
        return(list(M=M,A=A))
    }
# Loop over groups to get M and A.
    for (k in 1:nG) {
        group.ind <- (G==k)
        Z <- X[group.ind,]
        start<-cumsum(c(0,numslices))[k]+1
        end<-cumsum(numslices)[k]
        Scale <- if(pool) Sigma.pool else Sigma[,,k]
        if (sum(group.ind)>=3) {
            Z <- (X[group.ind,]-matrix(1,sum(group.ind),1) %*%
                 apply(X[group.ind,],2,mean)) %*% (Scale %^% (-1/2))
            temp <- psave1(Z,Y[group.ind],W[group.ind],info[[k]])
            M <- M+temp$M*sum(group.ind)/n
            A[start:end,,]<- temp$A*sqrt(sum(W[group.ind])/n)
        }
    }
# The following code computes the eigenvectors of M, and 
# transform the eigenvectors into X scale.
    D <- eigen(M)
    or <- rev(order(abs(D$values)))
    evectors <- apply(Sigma.pool%^%(-1/2)%*%D$vectors, 2,
                  function(x) x/sqrt(sum(x^2)))
    dimnames(evectors) <- 
         list(attr(dr.x(object),"dimnames")[[2]],
              paste("Dir",1:dim(evectors)[2],sep=""))
    # The function test.psave() tests the marginal coordinate hypothesis.
    # The argument H is a matrix transormed from the actual hypothesis.
    # Function is private.
    test.psave <- function(H) {
        r <- dim(H)[2]
        st <- 0
        h <- sum(numslices)
        for (j in 1:h) {
            st <- st + sum((t(H) %*% A[j, , ] %*% H)^2) * n/2
        }
        df <- (h - nG) * r * (r + 1)/2
        pv <- 1 - pchisq(st, df)
        return(data.frame(Test=st,df=df,P.value=pv))
    }

    # The function coordinate.test() tests the marginal coordinate hypothesis.
    # The argument gamma is the actual hypothesis.
    # Function is public.
    coordinate.test <- function(gamma) {
        r <- p - dim(gamma)[2]
        gamma <- Sigma.pool %^% (1/2) %*% gamma
        H <- as.matrix(qr.Q(qr(gamma), complete = TRUE)[, (p - r + 1):p])
        return(test.psave(H))
    }
    # The function test() does the sequential marginal dimension testing from 0 to d.
    # Function is public.
    test <- function(d) {
        tp <- data.frame()
        for (i in 1:d) 
            tp <- rbind(tp,test.psave(D$vectors[,(i):p]))
        rr<-paste(0:(d-1),"D vs >= ",1:d,"D",sep="")
        dimnames(tp)<-list(rr,c("Stat","df","p-values"))
        return(tp)
    }
    return(c(object,
     list(evectors=evectors,evalues=D$values[or],slice.info=slice.info,
          test=test,coordinate.test=coordinate.test,M=M)))
}

dr.test.psave <- function(object,numdir=object$numdir,...) 
  object$test(numdir)

dr.coordinate.test.psave <- function(object,hypothesis,d=NULL,...) {
    gamma <- if (class(hypothesis) == "formula")
        coord.hyp.basis(object, hypothesis)
        else as.matrix(hypothesis)
    object$coordinate.test(gamma)
}

summary.psave <- function(object,...) {
 ans <- summary.psir(object,...)
 ans$method <- "psave"
 return(ans)
 }
 
drop1.psave <- function(object,scope=NULL,...){
    drop1.dr(object,scope,...) }
