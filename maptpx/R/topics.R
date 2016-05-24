##### Estimation for Topic Models ######

## intended main function; provides defaults and selects K via marginal lhd
topics <- function(counts, K, shape=NULL, initopics=NULL, tol=0.1, 
                   bf=FALSE, kill=2, ord=TRUE, verb=1, ...)
  ## tpxselect defaults: tmax=10000, wtol=10^(-4), qn=100, grp=NULL, admix=TRUE, nonzero=FALSE, dcut=-10
{
  X <- CheckCounts(counts)
  p <- ncol(X) 
  if(verb>0)
    cat(sprintf("\nEstimating on a %d document collection.\n", nrow(X)))

  ## check the prior parameters for theta
  if(prod(shape>0) != 1){ stop("use shape > 0\n") }
                
  ## check the list of candidate K values
  if(prod(K>1)!=1){ stop(cat("use K values > 1\n")) }
  K <- sort(K)
 
  ## initialize
  initopics <- tpxinit(X[1:min(ceiling(nrow(X)*.05),100),], initopics, K[1], shape, verb)
  
  ## either search for marginal MAP K and return bayes factors, or just fit
  tpx <- tpxSelect(X, K, bf, initopics, alpha=shape, tol, kill, verb, ...)
  K <- tpx$K
  
  ## clean up and out
  if(ord){ worder <- order(col_sums(tpx$omega), decreasing=TRUE) } # order by decreasing usage
  else{ worder <- 1:K }
  ## Main parameters
  theta=matrix(tpx$theta[,worder], ncol=K, dimnames=list(phrase=dimnames(X)[[2]], topic=paste(1:K)) )
  omega=matrix(tpx$omega[,worder], ncol=K, dimnames=list(document=NULL, topic=paste(1:K)) )
  if(nrow(omega)==nrow(X)){ dimnames(omega)[[1]] <- dimnames(X)[[1]] }
  
  ## topic object
  out <- list(K=K, theta=theta, omega=omega, BF=tpx$BF, D=tpx$D, X=X)
  class(out) <- "topics"
  invisible(out) }

 
## S3 method predict function
predict.topics <- function(object, newcounts, loglhd=FALSE, ...)
  ## tpxweights optional arguments and defauls are verb=FALSE, nef=TRUE, wtol=10^{-5}, tmax=1000
{
  if(is.vector(newcounts)){ newcounts <- matrix(newcounts, nrow=1) }
  if(class(newcounts)[1] == "TermDocumentMatrix"){ newcounts <- t(newcounts) }
  X <- as.simple_triplet_matrix(newcounts)
  
  if(!(class(object)%in%c("topics","matrix"))){ stop("object class must be `topics' or 'matrix'.") }

  if(class(object)=="topics"){
    theta <- object$theta
    if(nrow(theta) != ncol(X)){ stop("Dimension mismatch: nrow(theta) != ncol(X)") }
    if(nrow(object$X) != nrow(object$omega)) # simple mixture
      { Q <- matrix(tpxMixQ(X, omega=object$omega, theta=theta, ...)$lQ, ncol=ncol(theta))
        return( (1:ncol(theta))[apply(Q,1,which.max)] ) } 
  }
  else{ theta <- object }

  start <- tpxOmegaStart(X=X, theta=theta)
  
  ## re-order xvec in doc-blocks, and build indices
  doc <- c(0,cumsum(as.double(table(factor(X$i, levels=c(1:nrow(X)))))))
  xvo <- X$v[order(X$i)]
  wrd <- X$j[order(X$i)]-1

  W <- tpxweights(n=nrow(X), p=ncol(X), xvo=xvo, wrd=wrd, doc=doc, start=start, theta=theta, ...)

  if(loglhd){
    L <- sum( X$v*log(tpxQ(theta=theta, omega=W, doc=X$i, wrd=X$j)) )
    return(list(W=W, L=L)) }
  else{ return(W) }
}

## S3 method summary function
summary.topics <- function(object, nwrd=5, tpk=NULL, verb=TRUE, ...){

    
  K <- object$K
  if(is.null(tpk)){ tpk <- 1:K }
  else if(prod(tpk %in% 1:K)!=1){ stop("requested tpk's are not in 1:K") }
 
  if(nwrd>0){ if(verb){ cat(paste("\nTop", nwrd, "phrases by topic-over-null term lift (and usage %):\n\n")) }
              Q0 <- col_sums(object$X)/sum(object$X)
              topwords <- c()
              toplift <- c()
              usage <- round( (col_means(object$omega)*100), 1)
              for(k in tpk){
                odds <- (log(object$theta[,k]) - log(Q0))[Q0!=0]
                ko <- order(odds, decreasing=TRUE)
                topk <- dimnames(object$theta)[[1]][Q0!=0][ko[1:nwrd]]
                topwords <- cbind(topwords, topk)
                toplift <- cbind(toplift, exp(sort(odds, decreasing=TRUE)[1:nwrd]))
                if(verb){ cat(paste("[",k, "] '", sep=""))
                          cat(topk, sep="', '")
                          cat(paste("' (",usage[k],") \n",sep="")) }
              }
              topwords <- as.matrix(topwords)
              toplift <- as.matrix(toplift)
              dimnames(topwords)[[2]] <- dimnames(toplift)[[2]] <- paste("topic",tpk,sep="")
              dimnames(toplift)[[1]] <- dimnames(toplift)[[1]] <- 1:nwrd
            }
  else{ topwords <- toplift <- NULL }
              
  if(!is.null(object$BF) && !is.null(object$D) && verb)
    {
      cat("\nLog Bayes factor and estimated dispersion, by number of topics:\n\n")
      D <- rbind(logBF=object$BF,Disp=object$D[1,])
      if(verb>1){ D <- rbind(D,p.val=object$D[2,]) }
      print(round(D,2))
      cat(paste("\nSelected the K =",object$K,"topic model\n\n"))
    }
  else if(verb){
    n <- nrow(object$K)
    rho = pchisq(sum(object$r^2), df=length(object$r), lower.tail=FALSE)
    cat("\nDispersion = ", round(object$D$dispersion,2) ,"\n\n",sep="") }

  retobj <- data.frame(topic=rep(tpk, each=nwrd), phrase=c(topwords), lift=c(toplift))
  invisible(retobj)
}


## Colors for topic plotting
TOPICOLS <- matrix(nrow=6,
                   c(grey(c(.9,.8,.7,.55,.35,0)), #GREY
                     "#FFCCCC", "#FA8072", "#FF5333", "#EE0000", "#CC1100", "#800000", #RED
                     "#BDFCC9", "#98FB98", "#49E20E", "#009900", "#006400", "#004F00", #GREEN	
                     "#BBFFFF", "#67E6EC", "#00C5CD", "#0198E1", "#0147FA",  "#000080"))  #BLUE

              	

plot.topics <- function(x, type=c("weight","resid"), group=NULL, labels=NULL, 
                        col=NULL, xlab=NULL, ylab=NULL, main=NULL, tpk=NULL,
                        lgd.K=NULL, cex.lgdc = 1, cex.lgdt = 1, cex.rmar= 1, ...){

  if(type[1]=="resid"){

    if(is.null(col[1])){col <- 8}
    if(is.null(xlab)){ xlab="abs( adusted residuals )" }
    if(is.null(main)){ main="" }

    resids <- tpxResids(x$X, theta=x$theta, omega=x$omega, ...)$r

    hist(resids, col=col, border=grey(.9),
         xlab=xlab,
         main="", cex.lab=cex.lgdt, font.lab=3)
    
    return(invisible())
  }
  
  n <- nrow(x$omega)
  mmm <- n==nrow(x$X)
  if(n==1){
    if(is.null(main)){ main="" }
    if(is.null(xlab)){ xlab="topic" }
    if(is.null(ylab)){ ylab="weight" }
    if(is.null(col)){ col = 8 }
    return(plot(c(x$omega), type="h", lwd=10, col=col, main=main, ylab=ylab, xlab=xlab)) }

  if(is.null(tpk)){ tpk <- 1:x$K }
  if(is.null(lgd.K)){ lgd.K <- max(.1*length(tpk),.75) }
  
  if(is.null(group) || !mmm){
    ltop = .65*n
    lbot = .35*n
    
    if(is.null(col)){ col<-1 }
    tpcl <- c(0, TOPICOLS[1:6,col[1]%%5])
    W <- x$omega[,tpk]
    if(mmm){
      brks <- seq(0,1,length=8)
      tplg <- c("w=0","w=1") }
    else{
      brks <- seq(min(W),max(W),length=8)
      tplg <- c(round(min(W),3),round(max(W),3))
    }
  } else{ # bunch of extra commands to get shading for two groups
    group <- as.factor(group)
    ltop = .85*n
    lbot = .15*n
    
    if(length(group)!=n){ stop("Your group membership length doesn't match omega.") }
    if(nlevels(group)!=2){ stop("Sorry, grouping must be a two-level factor") }
    if(is.null(col) || length(col)<2){ col <- 1:2 }
    
    tpcl <- c(TOPICOLS[6:1,col[1]%%5],"#FFFFFF",TOPICOLS[,col[2]%%5])
    brks <- c(seq(-1, -0.1,length=7),seq(0.1,1,length=7))
    W <- x$omega[,tpk]*(c(-1,1)[as.logical(group)+1])
    if(is.null(labels)){labels <- c("F","T")}
    tplg=rep("w=1",2)
  }

  ## plot parameters
  xlg <- length(tpk)+lgd.K
  old.mar <- par()$mar
  par(xpd=TRUE, mar=c(5.1,4.1,2.1,5*cex.rmar))
  if(is.null(ylab)){
    if(nrow(x$X)==nrow(x$omega)){ ylab="Document" }
    else{ ylab="group" } }
  if(is.null(xlab)){ xlab="Topic" }
  if(is.null(main)){ main="Topic-Loading Weights" }

  ## finally plot
  image(y=1:n, x=1:length(tpk), t(W), ylab=ylab, xlab=xlab,
        main=main, col=tpcl, font.lab=3, xaxt="n", yaxt="n", breaks=brks, ...)
  axis(side=1, at=1:length(tpk), labels=tpk, tick=FALSE, line=-.5)

  if(!mmm){ axis(side=2, at=1:n) }
  else{axis(2)}
  
  points(rep(xlg,length(tpcl)), seq(lbot,ltop,length=length(tpcl)), col=tpcl, pch=15, cex=3*cex.lgdc)
  text(rep(xlg,2), y=c(lbot-.08*n, ltop+.08*n), tplg, cex=cex.lgdt)
  if(!is.null(labels)){ text(rep(xlg,2), y=c(lbot-.14*n, ltop+.14*n), labels, font=3, cex=cex.lgdt) }
  par(mar=old.mar)
}

## logit and expit, treating first element as null
logit <- function(prob){
  f <- function(p) .C("Rlogit", d=as.integer(d), 
    eta=double(d-1), prob=as.double(p), PACKAGE="maptpx")$eta
  prob[prob < 1e-10] <- 1e-10
  prob[1-prob < 1e-10] <- 1-1e-10
  if(is.matrix(prob) || is.data.frame(prob)){
    d <- ncol(prob)
    eta <- matrix(t(apply(prob,1,f)), ncol=d-1, dimnames=list(row=dimnames(prob)[[1]], col=dimnames(prob)[[2]][-1]))
  }
  else{
    d <- length(prob)
    eta <- f(prob)
  }
  return(eta) }

expit <- function(eta){
  f <- function(e) .C("Rexpit", d=as.integer(d), prob=double(d), 
    eta=as.double(e), PACKAGE="maptpx")$prob
  eta[eta==Inf] <- 1e10
  eta[eta==-Inf] <- -1e10
  if(is.matrix(eta) || is.data.frame(eta)){
    d <- ncol(eta)+1
    if(is.null(dimnames(eta)[[2]])) dimnames(eta)[[2]] <- paste("nef",1:ncol(eta),sep='')
    prob <- matrix(t(apply(eta,1,f)), ncol=d, dimnames=list(row=dimnames(eta)[[1]], col=c('nullcat',dimnames(eta)[[2]])))
  }
  else{
    d <- length(eta)+1
    prob <- f(eta)
  }
  return(prob) }


### topic weight variance matrix (in NEF parametrization)
topicVar <- function(counts, theta, omega){
  X <- CheckCounts(counts)
  if(nrow(omega) != nrow(X)) stop("omega does not match counts")
  if(ncol(omega) != ncol(theta)) stop("omega does not match theta")

  K <- ncol(omega)
  n <- nrow(X)
  p <- nrow(theta)
  
  q <- tpxQ(theta=theta, omega=omega, doc=X$i, wrd=X$j)
  H <- array(.C("RnegHW",
                  n = as.integer(n),
                  p = as.integer(p),
                  K = as.integer(K-1),
                  omeg = as.double(omega[,-1]),
                  thet = as.double(theta[,-1]),
                  doc = as.integer(X$i-1),
                  wrd = as.integer(X$j-1),
                  cnt = as.double(X$v),
                  q = as.double(q),
                  N = as.integer(length(q)),
                  H = double(n*(K-1)^2),
                  PACKAGE="maptpx")$H,
               dim=c(K-1,K-1,n))

  S <- array(apply(H, 3, function(h) tryCatch(solve(h), error=function(e) solve(h + diag(.00001,K-1))) ),
             dim=c(K-1,K-1,n), dimnames=list(topics=1:(K-1), topics=1:(K-1), doc=dimnames(X)[[1]]) )
  return(S)
}
