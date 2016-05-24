PCArot<-function(obj,dim,itermax=100,graph=TRUE)
{
  cl <- match.call()
  if (!inherits(obj, "PCAmix")) 
    stop("use only with \"PCAmix\" objects")
  if (dim <= 1) stop("'dim' must be an interger greater or equal to 2")
  if (dim > ncol(obj$A)) stop("'dim' must be an integer smaller or equal to 'ndim' specified in 'PCAmix'")
  A <- obj$A[,1:dim] 
  V<- obj$V[,1:dim] 
  V<-sqrt(obj$M)*V
  scores.stand <- obj$scores.stand[,1:dim] 
  indexj <- obj$rec$indexj
  p1 <- obj$rec$p1
  p <- obj$rec$p
  p2 <- p-p1
  Z <- obj$rec$Z
  G <- obj$rec$G
  m <- ncol(Z)-p1
  X.quali <- obj$rec$X.quali
  n <- obj$rec$n
  
  if (dim==2)
  {
    res <- sol.2dim(A,indexj,p,p1)
    theta <- res$theta
    iter <- 1
    T <- res$T
    A.rot <- A%*%T
    scores.stand.rot <- scores.stand%*%T
  } else {
    matcombn <- combn(1:dim,2) 
    nbpaires <- ncol(matcombn) 
    A.rot <- A
    scores.stand.rot <- scores.stand
    cptthetazero <- 0
    iter <- 0
    while (cptthetazero<nbpaires){
      cptthetazero<-0
      iter <- iter+1
      for (j in 1:nbpaires){
        indicescol <- matcombn[,j]
        res <- sol.2dim(A.rot[,indicescol],indexj,p,p1)
        theta <- res$theta
        T <- res$T
        A.rot[,indicescol] <- A.rot[,indicescol]%*%T
        scores.stand.rot[,indicescol] <- scores.stand.rot[,indicescol] %*%T
        if (round(theta,digits=3)==0)
          cptthetazero <- cptthetazero+1  
      }
      if (iter>itermax)	stop("Stop: maximum number of iterations reached.")	
    }
    theta <- NULL
    T <- t(scores.stand)%*%scores.stand.rot/n
  }
  A1 <- NULL
  A2 <- NULL
  A2coord <- NULL
  C <- NULL
  ps <- NULL
  if(p1>0){
    A1 <- A.rot[1:p1,] 
    colnames(A1) <- paste("dim", 1:dim, sep = "",".rot")
  }
  if (p1!=p)	{
    if(p1>0) A2 <- A.rot[-(1:p1),] else A2 <- A.rot
    ns <- apply(G,2,sum)
    ps <- ns/n
    A2coord <- sweep(A2,MARGIN=1,STATS=sqrt(ps),FUN="/") 
    colnames(A2coord) <- paste("dim", 1:dim, sep = "",".rot")
    C <- matrix(NA,(p-p1),dim)
    rownames(C) <- colnames(X.quali)
    colnames(C) <- paste("dim", 1:dim, sep = "",".rot")
    for (j in 1:(p-p1)) {
      C[j,] <- apply(A2[which(indexj==(j+p1))-p1,],2,FUN=function(x) {sum(x^2)}) }
  }	
  sqload.rot <- rbind(A1^2,C)
  var.rot <- apply(sqload.rot,2,sum)
  eig <- matrix(0,dim,2)
  colnames(eig) <- c("Variance","Proportion")
  rownames(eig) <- names(var.rot)
  eig[,1] <- var.rot
  eig[,2] <- eig[,1]/sum(obj$eig[,1])*100
  scores.rot <- sweep(scores.stand.rot,2,STATS=sqrt(var.rot),FUN="*")
  V.rot <- V%*%diag(1/sqrt(obj$eig[1:dim,1]))%*%T%*%diag(sqrt(eig[,1]))
  colnames(sqload.rot) <- names(var.rot) <- colnames(scores.rot) <- colnames(scores.stand.rot)<- colnames(V.rot)<-paste("dim", 1:dim, sep = "",".rot")
  rownames(V.rot) <- colnames(Z)
  coef <- structure(vector(mode = "list", length = dim), names = paste("dim", 1:dim, sep = "",".rot"))
  gc <- obj$rec$g #gravity center
  sdev <- obj$rec$s #standard deviations
  for (g in 1:dim)  {
    beta <- matrix(NA,p1+m+1,1)
    beta[1,1] <- -sum(V.rot[,g]*gc/sdev)
    beta[2:(p1+m+1),1] <- V.rot[,g]/sdev
    rownames(beta) <- c("const",colnames(Z))
    coef[[g]] <- beta
  }
  V.rot<-(obj$M)^(-1/2)*V.rot
  
  res.ind<-list(coord=scores.rot)
  res.quanti<-list(coord=A1)
  res.levels<-list(coord=A2coord)
  res.quali<-list(coord=C)
  res <- list(call = cl,theta=theta,T=T,eig=eig,ind=res.ind,quanti=res.quanti,levels=res.levels,quali=res.quali,
              coef=coef,sqload=sqload.rot, scores.stand=scores.stand.rot,scores=scores.rot,
              categ.coord=A2coord,quanti.cor=A1,quali.eta2=C,ndim=dim,rec=obj$rec,iter=iter,V=V.rot)
  class(res) <- "PCAmix"
  if (graph) {
    plot.PCAmix(res,main="Scores after rotation")
    if (p1!=p) {
      dev.new()
      plot.PCAmix(res,choice="levels",main="Categories after rotation")
    }
    if (p1!=0) {
      dev.new()
      plot.PCAmix(res,choice="cor",main="Correlation circle after rotation")
    }
    dev.new()
    plot.PCAmix(res,choice="sqload",main="Squared loadings after rotation")
  }
  return(res)	
}
