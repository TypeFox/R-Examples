test.equal <- function(obj) {
  if(!inherits(obj, "pendensity")) break("Not a object of class pendensity")
  if(is.null(obj$values$x)) stop("The input is an univariate pendensity-object. Nothing to test")
  help.env <- distr.func.help(obj)
  K <- obj$splines$K
  q <- obj$splines$q
  Z <- obj$splines$Z
  x.factor <- obj$values$covariate$x.factor
  l.Z <- length(x.factor[,1])
  ck <- obj$results$ck
  M <- floor(K/2)+1
  for(i in 1:l.Z) {
    assign(paste("ZZ.help",i,sep=""),kronecker(diag(1,K-1),x.factor[i,]),help.env)
    assign(paste("c",i,sep=""),c.temp <- c(ck[i,]),help.env)
    assign(paste("C",i,sep=""),(diag(c.temp)-tcrossprod(c.temp))[,-M],help.env)
  }
  
  knots.val <- obj$splines$knots.val
  mat1 <- matrix(0,K+q,K)
  for(i in 1:K) {
    vec <- matrix(0,K+q,1)
    for(j in 1:q) {
      vec[j+i,1] <- poly.part(i,j,knots.val,help.env,q,yi=NULL,poly=FALSE)
    }
    vec <- cumsum(vec)/sum(vec)
    mat1[,i] <- vec
  }
  mat1 <- mat1[1:K,]

  var.par <- obj$results$variance.par

  choose.temp <- choose(l.Z,2)

  help <- seq(1:l.Z)
  pvalues <- matrix(0,choose.temp,1)
  values <- c()
  l.Z.help <- l.Z
  start <- l.Z
  ind <- 1
  while(start>1) {
    for(j in start:2) {
      values[ind] <- paste(start," vs. ",j-1,sep="")
      CC.help <- get(paste("C",j,sep=""),help.env)%*%t(get(paste("ZZ.help",j,sep=""),help.env))-(get(paste("C",j-1,sep=""),help.env)%*%t(get(paste("ZZ.help",j-1,sep=""),help.env)))
      W.temp <- mat1%*%CC.help%*%var.par%*%t(CC.help)%*%t(mat1)
      svdW <- svd(W.temp)
      Tmax.b <- matrix(0,2000,1)  
      c.diff <- get(paste("c",j,sep=""),help.env)-get(paste("c",j-1,sep=""),help.env)
      Tmat <- matrix(0,K,1)
      for(i in 1:K) Tmat[i,1] <- sum(c.diff%*%mat1[,i])
      Tmax <- max(abs(Tmat))
      n.set <- 1000
      indi.help <- matrix(0,n.set,1)
      for(b in 1:n.set) {
        zb <- rnorm(K,0,1)
        Uw <- svdW$u
        dw <- diag(sqrt(svdW$d))
        Tb <- Uw%*%dw%*%zb
        Tmax.b <- max(abs(Tb))
        if(Tmax.b>=Tmax) indi.help[b,1] <- 1
      }
      pvalues[ind,1] <- (1/n.set)*sum(indi.help)
      l.Z.help <- l.Z.help-1
      ind <- ind+1
    }
    start <- start-1
    l.Z.help <- l.Z.help-1
  }
  rownames(pvalues) <- values
    
  print(pvalues)
}
