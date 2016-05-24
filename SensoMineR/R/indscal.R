indscal <- function(matrice,matrice.illu=NULL,maxit=200,coord=c(1,2),eps = 1/10^5){

transformdw<-function(X){# transforme les distances en produit scalaire pour un poids des individus constant =1/n
    taille<-dim(X)[[1]]
    d2ip<-(1/taille)*(rowSums((X^2),na.rm=TRUE))
    d2pp<-(1/taille)*sum(d2ip,na.rm=TRUE)
    X1<-X
    for (i in 1:(taille)) {
        for (j in( i:taille)) X1[i,j]<-(1/2)*(d2ip[i]+d2ip[j]-X[i,j]^2-d2pp)
    }
    for (i in 1:(taille-1)) {
        for (j in( i+1):(taille))  X1[j,i]<-   X1[i,j]
    }
    return(X1)
}

dist2full <- function(dis){
        n <- attr(dis, "Size")
        full <- matrix(0, n, n)
        full[lower.tri(full)] <- dis
        full + t(full)
}

k=max(coord)
stimuli=rownames(matrice)
subjects=colnames(matrice)[2*(1:(ncol(matrice)/2))]
d <- array(0, dim=c(nrow(matrice),nrow(matrice),ncol(matrice)/2))
for (i in 1:dim(d)[3]) d[,,i] <- dist2full(dist(matrice[,(2*i-1):(2*i)]))

      call <- match.call()
      call2 <- match.call(expand.dots =FALSE)
      n <- dim(d)[1]
      m <- dim(d)[3]
      b <- array(0, dim(d))
      ac <- rep(0, m)
      W <- matrix(0, m, k)
      diff <- 1
      r2 <- 0
      delta <- 1
      it <- numeric(1)
      record <- numeric(maxit)
      mode(d) <- "single"
      storage.mode(n) <- "integer"
      for(i in 1:m) b[, , i] <- transformdw(d[,,i])

      b <- array(apply(b, 3, function(x){(x)/sqrt(sum(x^2))}), dim = dim(b))
      F1 <- matrix(b, byrow =TRUE, nrow = m)
      F2 <- matrix(b, nrow = n)
      XL <- cmdscale(apply(d, 1:2, median), k = max(k, 2))  #
#
# remove comment from next line to perturb starting position
# for testing:
     XL <- XL * rnorm(n * k, 2, 0.3)     #
#
      if(k == 1) XL <- t(as.matrix(XL[, 1]))
      XR <- XL
      while((delta) > eps) {
            it <- it + 1
            if(it > maxit) stop("Process failed to converge")
            G <- apply(array(c(XL, XR), dim = c(n, k, 2)), 2, function(x){kronecker(x[, 1], x[, 2])})
            W <- F1 %*% G %*% solve(t(G) %*% G)
            bpred <- array(t(W %*% t(G)), dim = c(n, n, m))
            record[it] <- (1 - sum((b - bpred)^2)/m)
            delta <- abs(r2 - record[it])
            r2 <- record[it]
            G <- matrix(0, m * n, k)
            for(i in 1:k) G[, i] <- kronecker(W[, i], XR[, i])
            XL <- F2 %*% G %*% solve(t(G) %*% G)      #
            for(i in 1:k) G[, i] <- kronecker(W[, i], XL[, i])
            XR <- F2 %*% G %*% solve(t(G) %*% G)
      }
      XL <- XR
      G <- apply(array(c(XL, XR), dim = c(n, k, 2)), 2, function(x){kronecker(x[, 1], x[, 2])})
      W <- F1 %*% G %*% solve(t(G) %*% G) #
      W <- replace(W,W<0,0)
      bpred <- array(t(W %*% t(G)), dim = c(n, n, m))
      r2 <- (1 - sum((b - bpred)^2)/m)
      r2ind <- 1 - apply(b - bpred, 3, function(x){sum(x^2)})     #
# final normalisation:
      X <- scale(XR)/sqrt((n - 1))
      W <- W/matrix(rep(((X/XR)[1,  ]^2),m),byrow=TRUE,nrow=m)    #
      rarr <- array(0, dim = c(n, k, m))
      z <- array(0, dim = c(n, k, m))
      pd <- matrix(0, m, (n * (n - 1))/2)
      stress.ind <- numeric(m)
      for(i in 1:m) {rarr[,  , i] <- X * sqrt(matrix(W[i,  ], n, k, byrow =TRUE))}
      names(r2ind) <- subjects
      dimnames(X) <- list(stimuli, paste("Dim", 1:k))
      dimnames(W) <- list(subjects, paste("Dim", 1:k))
dev.new()
  plot( W[,coord], xlim=c(0,1),ylim=c(0,1),type="n",xlab=paste("Dim",coord[1]),ylab=paste("Dim",coord[2]), main = "Weight representation")
  points(W[,coord[1]],W[,coord[2]],cex=1.2,pch=20)
  text( W[,coord[1]], W[,coord[2]], labels=subjects, cex = 0.8, pos = 4, offset = 0.2)

dev.new()
  X <- as.matrix(X[,])
  aa=cor(matrice,X[,coord])

  plot(0, 0, xlab = paste("Dim",coord[1]), ylab = paste("Dim",coord[2]), xlim = c(-1,1), ylim = c(-1,1), col = "white", asp=1, main="Correlation circle")
  x.cercle <- seq(-1, 1, by = 0.01)
  y.cercle <- sqrt(1 - x.cercle^2)
  lines(x.cercle, y = y.cercle)
  lines(x.cercle, y = -y.cercle)
  abline(v=0,lty=2)
  abline(h=0,lty=2)
  for (v in 1:nrow(aa)){
    arrows(0, 0, aa[v, 1], aa[v, 2], length = 0.1, angle = 15, code = 2)
    text(aa[v, 1], y = aa[v, 2], labels = rownames(aa)[v], pos = 1, offset=0.1)
  }

  if (length(matrice.illu)!=0){
    bb=cor(matrice.illu,X[,coord])
    for (v in 1:nrow(bb)){
      arrows(0, 0, bb[v, 1], bb[v, 2], length = 0.1, angle = 15, code = 2, col = "blue")
      text(bb[v, 1], y = bb[v, 2], labels = rownames(bb)[v], pos = 1, offset=0.1, col="blue")
    }
  }

dev.new()
  plot(X[,coord], main = "Stimuli map", xlab = paste("Dim",coord[1]), ylab = paste("Dim",coord[2]), asp=1, cex=0.8, pch=20)
  text( X[,coord[1]], X[,coord[2]], labels = stimuli, cex = 0.8, pos = 4, offset = 0.2)
  abline(v=0,lty=2)
  abline(h=0,lty=2)

     dimnames(rarr) <- list(stimuli, paste("Dim", 1:k), subjects)
      out <- list(W = as.matrix(W), points = as.matrix(X), subvar = r2ind, r2 = r2,
            dfr = (k * (m + n - 2))/((m * n * (n - 1))/2))
      class(out) <- "indscal"
      out
}
