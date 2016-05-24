tau.within <- function(x){
  
  
  if (!inherits(x, "netmeta"))
    stop("Argument 'x' must be an object of class \"netmeta\"")

  
  nmak <- nma.krahn(x)
  
  
  t <- diag(rep(1,length(as.character(nmak$design$design))))
  rownames(t) <- paste(as.character(nmak$design$comparison),
                       as.character(nmak$design$design),sep="_")
  colnames(t) <- paste(as.character(nmak$design$comparison),
                       as.character(nmak$design$design),sep="_")
  X_incon_studies <- t[paste(as.character(nmak$studies$comparison),
                             as.character(nmak$studies$design),sep="_"),]
  
  
  DiagonalBlocks <- function(V){
    if (length(V)==1) B<-1
    else{
      m <- NULL
      b <- NULL
      d <- NULL
      B <- NULL
      for (i in 1:nrow(V)){
        h <- which(abs(V[i,])>0)
        b[i] <- h[1]
        d[i] <- max(h)-min(h)+1
      }
      ##
      m <- d[as.numeric(names(table(b)))]
      b <- b[as.numeric(names(table(b)))]
      ##
      for (i in 1:length(m))
        B[i] <- list(b[i]:(b[i]+m[i]-1))
    }
    B
  }

  
  V0 <- function(m)
    diag(m)/2 + matrix(.5, m, m)
  
  
  W <- function(V, tausq){
    W <- diag(nrow(as.matrix(V)))
    b <- DiagonalBlocks(V)
    ##
    for (i in 1:length(b))
      W[b[[i]], b[[i]]] <- solve(as.matrix(V)[b[[i]], b[[i]]]+tausq*V0(length(b[[i]])))
    ##
    W
  }
  
  
  tau.DL <- function(X, y, V){
    C <- solve(t(X)%*%(W0<-W(V,0))%*%X)
    theta <- C%*%t(X)%*%W0%*%y
    r <- y-X%*%theta
    Q <- t(r)%*%W0%*%r
    V.0 <- diag(nrow(V))
    b <- DiagonalBlocks(V)
    ##
    for (i in 1:length(b))
      V.0[b[[i]],b[[i]]]<-V0(length(b[[i]]))
    ##
    trace1 <- sum(diag(W0%*%V.0))
    trace2 <- sum(diag(C%*%t(X)%*%W0%*%V.0%*%W0%*%X))
    ##
    if (trace1==trace2)
      tausq <- 0
    else
      tausq <- max(0, (Q-nrow(V)+ncol(X))/(trace1-trace2))
    ##
    sqrt(tausq)
  }
  
  
  tau.within <- tau.DL(X_incon_studies, nmak$studies$TE, nmak$V.studies)
  ##
  tau.within
}
