lplikint <-
function(jmptimes,jmpsizes=rep(1,length(jmptimes)),
                     Y=rep(1,length(jmptimes)),
                     K=function(x)3/4*(1-x^2)*(x<=1&x>= -1),
                     bw,adjust=1,nu=0,p=1,Tau=1,n=101,
                     tseq=seq(from=0,to=Tau,length=n),
                     tol=1e-5,maxit=100,us=10,gd=5){
  if(p<nu){p <- nu;warning("p<nu: p reset to nu!\n")}
  Ntau <- length(jmptimes)
  bw <- bw*adjust
  ##  W <- function(s,tt){1/bw*K((s-tt)/bw)}
  G <- function(x)x^(0:p) / gamma(0:p + 1)
  lplikint.fun <- function(tt){
    ##  if(tt==0.3)browser()
    G.mat <- sapply(jmptimes-tt,G)
    W.vec <- K((jmptimes-tt)/bw)/bw / Y
    valid <- !is.na(W.vec) & W.vec!=0
    mu.vec <- sapply(0:(p*2),
                     function(i)integrate(function(x)x^i*K(x),
                                          lower=-tt/bw,
                                          upper=(Tau-tt)/bw
                                          )$value
                     )
    d <- diag(bw^(0:p) / gamma(0:p + 1))
    score <- function(theta){
      denom <- t(G.mat)%*%theta
      G.mat[,valid,drop=FALSE] %*% (W.vec/denom)[valid] -
        d %*% mu.vec[0:p +1]
    }
    mHessian <- function(theta){
      denom <- t(G.mat)%*%theta
      G.mat[,valid,drop=FALSE] %*% diag(as.numeric((W.vec/denom^2))[valid]) %*%
        t(G.mat[,valid,drop=FALSE])
    }
    it <- 0;
    init <- Ntau/sum(c(Y,Y[Ntau])*diff(c(0,jmptimes,Tau)))*
      us^(1/gd * -gd:gd) ##us=upper search limit
    init.v <- apply(apply(cbind(init,matrix(0,nrow=gd*2+1,ncol=p)),1,score),
                    1,function(x)sum(abs(x)))

    theta1 <- c(init[which.min(init.v)],rep(0,p))
    theta0 <- rep(Inf,p+1)
    while(sqrt(sum((theta1-theta0)^2))>tol && it<maxit){
      theta0 <- theta1
      theta1 <- theta0+solve(mHessian(theta0),score(theta0))
      it <- it+1
    }
    if(it==maxit)
      warning(paste("max number of iterations reached at t=",tt,"!\n"))
    list(theta=theta1,vcov=solve(mHessian(theta1)),bw=bw*Tau)
  }
  raw.res <- lapply(tseq,lplikint.fun)
  list(x=tseq,
       y=sapply(seq(along=tseq),function(i)raw.res[[i]]$theta[nu+1]),
       se=sapply(seq(along=tseq),function(i)raw.res[[i]]$vcov[nu+1,nu+1]^0.5),
       bw=## sapply(seq(along=tseq),function(i)raw.res[[i]]$bw),
         bw,
       fun=lplikint.fun
       )         
}

