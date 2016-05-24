PCAmix<- function (X.quanti=NULL,X.quali=NULL,ndim=5,rename.level=FALSE,weight.col=NULL,weight.row=NULL,graph=TRUE)
{
  cl <- match.call()
  rec <- recod(X.quanti, X.quali,rename.level)
  n <- rec$n
  p <- rec$p
  p1 <- rec$p1
  p2 <- p - p1
  X <- rec$X
  G <- rec$G
  W <- rec$W
  m <- ncol(W) - p1
  q <- qr(W)$rank
  indexj <- rec$indexj
  if (!is.null(X.quali)) 
  {
    ns <- apply(G, 2, sum)
    ps <- ns/nrow(G)
  }
  else {
    ns <- NULL
    ps <- NULL
  }
  M1 <- c(rep(1, p1))
  M2 <- ps
  M2.inv <- 1/M2
  N <- rep(1/n, n)
  if (!is.null(weight.col) == TRUE) 
  {
    weight.col.quant <- weight.col[1:p1]
    weight.col.qual <- weight.col[(p1 + 1):(p1 + m)]
    M1 <- M1 * weight.col.quant
    M2.inv <- M2.inv * weight.col.qual
  }
  if (!is.null(weight.row) == TRUE) 
  {
    N <- (N * weight.row)/sum(weight.row)
  }
  Met.global <- c(M1, M2.inv)
  names(Met.global) <- colnames(W)
  e <- svd.triplet(W, N, Met.global)
  eig <- matrix(0, q, 3)
  colnames(eig) <- c("Eigenvalue", "Proportion", "Cumulative")
  rownames(eig) <- paste("dim", 1:q, sep = " ")
  eig[, 1] <- e$vs[1:q]^2
  eig[, 2] <- 100 * eig[, 1]/sum(eig[, 1], na.rm = T)
  eig[1, 3] <- eig[1, 2]
  if (q > 1) 
  {
    for (j in 2:q) eig[j, 3] <- eig[j, 2] + eig[j - 1, 3]
  }
  if (ndim <= 1) 
    stop("\"ndim\" must be an interger greater or equal to 2")
  ndim <- min(ndim, q)
  U <- data.frame(e$U[, 1:ndim])
  rownames(U) <- rownames(W)
  colnames(U) <- rownames(eig)[1:dim(U)[2]]
  d <- e$vs[1:ndim]
  V.total.dim <- data.frame(e$V)
  U.total.dim <- data.frame(e$U)
  d.total.dim <- e$vs
  if (q != 1) {
    F <- as.matrix(U) %*% diag(d)
    F.total.dim <- as.matrix(U.total.dim) %*% diag(d.total.dim)
    colnames(F) <- rownames(eig)[1:dim(F)[2]]
  }
  else {
    F <- data.frame(U * d)
    F.total.dim <- data.frame(U.total.dim * d.total.dim)
    colnames(F) <- rownames(eig)[1:dim(F)[2]]
  }
  V <- data.frame(e$V[, 1:ndim])
  A <- NULL
  A1 <- NULL
  A2 <- NULL
  C <- NULL
  if (p1 > 0 & p2 == 0) {
    V1 <- as.matrix(V[1:p1, ])
    V1.total.dim <- as.matrix(V.total.dim[1:p1, ])
    if (p1 > 1) 
    {
      A1 <- V1 %*% diag(d)
      A1.total.dim <- V1.total.dim %*% diag(d.total.dim)
    }
    else {
      A1 <- data.frame(V1 * d)
      A1.total.dim <- data.frame(V1.total.dim * d.total.dim)
    }
    colnames(A1) <- paste("dim", 1:ndim, sep = "")
    rownames(A1) <- colnames(W)[1:p1]
    contrib.quanti <- sweep(A1^2, 1, STATS = M1, FUN = "*")
    contrib.quanti.pct <- sweep(contrib.quanti, 2, STATS = d^2, 
                                FUN = "/")
    colnames(contrib.quanti) <-colnames(contrib.quanti.pct) <- paste("dim", 1:ndim, sep = "")
    contrib.moda <- NULL
    contrib.moda.pct <- NULL
    A <- A1
    rownames(A) <- colnames(W)
    cos2.quanti <- sweep(A1^2, 1, STATS = apply(A1.total.dim, 
                                                1, function(v) {
                                                  return(sum(v^2))
                                                }), FUN = "/")
    contrib.moda <- NULL
    cos2.moda <- NULL
    rownames(A) <- colnames(W)
    rownames(contrib.quanti) <- rownames(contrib.quanti.pct) <-rownames(cos2.quanti) <- colnames(W)
    colnames(contrib.quanti) <- colnames(contrib.quanti.pct) <- colnames(cos2.quanti) <- colnames(A)
  }
  if (p1 == 0 & p2 > 0) {
    V2 <- as.matrix(V[(p1 + 1):(p1 + m), ])
    V2.total.dim <- as.matrix(V.total.dim[(p1 + 1):(p1 + 
                                                      m), ])
    if (p2 > 1) {
      A2 <- diag(M2.inv) %*% V2 %*% diag(d)
      A2.total.dim <- diag(M2.inv) %*% V2.total.dim %*% 
        diag(d.total.dim)
    }
    else {
      A2 <- data.frame(diag(M2.inv) %*% V2 * d)
      A2.total.dim <- data.frame(diag(M2.inv) %*% V2.total.dim * 
                                   d.total.dim)
    }
    if (!is.null(weight.col)) {
      A2 <- sweep(A2, 1, STATS = weight.col.qual, FUN = "/")
      A2.total.dim <- sweep(A2.total.dim, 1, STATS = weight.col.qual, 
                            FUN = "/")
    }
    colnames(A2) <- paste("dim", 1:ndim, sep = "")
    rownames(A2) <- colnames(W)[(p1 + 1):(p1 + m)]
    contrib.moda <- sweep(A2^2, 1, STATS = ps, FUN = "*")
    contrib.moda.pct <- sweep(contrib.moda, 2, STATS = d^2, FUN = "/")
    if (!is.null(weight.col)) {
      contrib.moda <- sweep(contrib.moda, 1, STATS = weight.col.qual, 
                            FUN = "*")
      contrib.moda.pct <- sweep(contrib.moda, 2, STATS = d^2, FUN = "/")
    }
    colnames(contrib.moda) <- colnames(contrib.moda.pct) <-paste("dim", 1:ndim, sep = "")
    C <- matrix(NA, p2, length(d))
    rownames(C) <- colnames(X.quali)
    colnames(C) <- paste("dim", 1:ndim, sep = "")
    for (j in 1:(p - p1)) {
      C[j, ] <- apply(as.matrix(contrib.moda.pct[which(indexj == 
                                                         j), ]), 2, FUN = sum)
    }
    C <- sweep(C, 2, STATS = d^2, FUN = "*")
    contrib.quanti <- NULL
    contrib.quanti.pct <- NULL
    A <- diag(sqrt(ps)) %*% as.matrix(A2)
    cos2.moda <- sweep(A2^2, 1, STATS = apply(A2.total.dim, 
                                              1, function(v) {
                                                return(sum(v^2))
                                              }), FUN = "/")
    contrib.quanti <- NULL
    cos2.quanti <- NULL
    rownames(A) <- colnames(W)
    rownames(contrib.moda) <-rownames(contrib.moda.pct) <- colnames(W)
    colnames(contrib.moda) <-colnames(contrib.moda.pct) <- colnames(A)
  }
  if (p1 > 0 & p2 > 0) {
    V1 <- as.matrix(V[1:p1, ])
    V1.total.dim <- as.matrix(V.total.dim[1:p1, ])
    V2 <- as.matrix(V[(p1 + 1):(p1 + m), ])
    V2.total.dim <- as.matrix(V.total.dim[(p1 + 1):(p1 + 
                                                      m), ])
    A1 <- (V1) %*% diag(d)
    A1.total.dim <- (V1.total.dim) %*% diag(d.total.dim)
    colnames(A1) <- paste("dim", 1:ndim, sep = "")
    rownames(A1) <- colnames(W)[1:p1]
    A2 <- diag(M2.inv) %*% V2 %*% diag(d)
    A2.total.dim <- diag(M2.inv) %*% V2.total.dim %*% diag(d.total.dim)
    if (!is.null(weight.col)) {
      A2 <- sweep(A2, 1, STATS = weight.col.qual, FUN = "/")
      A2.total.dim <- sweep(A2.total.dim, 1, STATS = weight.col.qual, 
                            FUN = "/")
    }
    colnames(A2) <- paste("dim", 1:ndim, sep = "")
    rownames(A2) <- colnames(W)[(p1 + 1):(p1 + m)]
    contrib.quanti <- sweep(A1^2, 1, STATS = M1, FUN = "*")
    contrib.quanti.pct <- sweep(contrib.quanti, 2, STATS = d^2, 
                                FUN = "/")
    colnames(contrib.quanti) <- colnames(contrib.quanti.pct) <- paste("dim", 1:ndim, sep = "")
    contrib.moda <- sweep(A2^2, 1, STATS = ps, FUN = "*")
    contrib.moda.pct <- sweep(contrib.moda, 2, STATS = d^2, FUN = "/")
    if (!is.null(weight.col)) {
      contrib.moda <- sweep(contrib.moda, 1, STATS = weight.col.qual, 
                            FUN = "*")
      contrib.moda.pct <- sweep(contrib.moda, 2, STATS = d^2, FUN = "/")
    }
    colnames(contrib.moda) <- colnames(contrib.moda.pct) <- paste("dim", 1:ndim, sep = "")
    C <- matrix(NA, p2, length(d))
    rownames(C) <- colnames(X.quali)
    colnames(C) <- paste("dim", 1:ndim, sep = "")
    for (j in 1:(p - p1)) {
      C[j, ] <- apply(contrib.moda.pct[which(indexj == (j + 
                                                          p1)) - p1, ], 2, FUN = sum)
    }
    C <- sweep(C, 2, STATS = d^2, FUN = "*")
    cos2.moda <- sweep(A2^2, 1, STATS = apply(A2.total.dim, 
                                              1, function(v) {
                                                return(sum(v^2))
                                              }), FUN = "/")
    A <- rbind(A1, diag(sqrt(ps)) %*% A2)
    cos2.quanti <- sweep(A1^2, 1, STATS = apply(A1.total.dim, 
                                                1, function(v) {
                                                  return(sum(v^2))
                                                }), FUN = "/")
    rownames(A) <- colnames(W)
    rownames(contrib.quanti) <- rownames(contrib.quanti.pct) <- rownames(cos2.quanti) <- colnames(W)[1:p1]
    rownames(contrib.moda) <- rownames(contrib.moda.pct) <- colnames(W)[(p1 + 1):(p1 + 
                                                                                    m)]
    colnames(contrib.moda) <-colnames(contrib.moda.pct) <- colnames(contrib.quanti)<- colnames(contrib.quanti.pct) <- colnames(cos2.quanti) <- colnames(A)
  }
  if (p2 > 0) {
    contrib.quali <- matrix(NA, p2, length(d))
    for (j in 1:(p - p1)) {
      contrib.quali[j, ] <- apply(as.matrix(contrib.moda[which(indexj == 
                                                                 (j + p1)) - p1, ]), 2, FUN = sum)
    }
    contrib.quali.pct<-sweep(contrib.quali, 2, STATS = d^2, FUN = "/")
    colnames(contrib.quali) <-colnames(contrib.quali.pct) <- colnames(C)
    rownames(contrib.quali) <-rownames(contrib.quali.pct) <- rownames(C)
  }
  else {
    contrib.quali <- NULL
    contrib.quali.pct <- NULL
  }
  contrib.ind<-(1/n) * F^2
  contrib.ind.pct <- sweep(contrib.ind, 2, STATS = d^2, FUN = "/")
  
  cos2.ind <- sweep(F^2, 1, STATS = apply(F.total.dim, 1, function(v) {
    return(sum(v^2))
  }), FUN = "/")
  rownames(contrib.ind) <- rownames(contrib.ind.pct) <- rownames(cos2.ind) <- rownames(W)
  colnames(contrib.ind) <- colnames(contrib.ind.pct) <-colnames(cos2.ind) <- colnames(A)
  result.ind <- list(coord = F, contrib = contrib.ind, contrib.pct = 100 * contrib.ind.pct,
                     cos2 = cos2.ind)
  result.quanti <- list(coord = A1, contrib= contrib.quanti, contrib.pct = 100 * contrib.quanti.pct, 
                        cos2 = cos2.quanti)
  result.levels <- list(coord = A2, contrib=contrib.moda, contrib.pct = 100 * contrib.moda.pct, 
                        cos2 = cos2.moda)
  result.quali<-list(contrib = contrib.quali, contrib.pct=contrib.quali.pct)
  
  rownames(V) <- rownames(A)
  colnames(V) <- colnames(A)
  sqload <- rbind(A1^2, C)
  names(d) <- colnames(U) <- colnames(F) <- colnames(sqload) <- colnames(cos2.ind) <- colnames(contrib.ind) <- paste("dim", 
                                                                                                                     1:ndim, sep = "")
  
  coef <- structure(vector(mode = "list", length = ndim), names = paste("dim", 
                                                                        1:ndim, sep = ""))
  ind.quant<-1:p1
  ind.qual<-(p1+1):(p1+m)
  
  gc<-rec$g
  gc.quant<-gc[ind.quant]
  gc.qual<-gc[ind.qual]
  
  sdev<-rec$s
  sdev.quant<-sdev[ind.quant]
  sdev.qual<-sdev[ind.qual]
  
  V.quant<-V[ind.quant, ,drop=FALSE]
  V.qual<-V[ind.qual, ,drop=FALSE]
  
  M<-Met.global
  M.quant<-M[ind.quant]
  M.qual<-M[ind.qual]
  
  for (g in 1:ndim) {
    beta <- matrix(NA, p1 + m + 1, 1)
    
    if (p1>0 & p2>0){
      beta[1, 1] <- -(sum(V.quant[, g]* M.quant * gc.quant/sdev.quant) + sum(V.qual[, g]* M.qual * gc.qual))
      beta[2:(p1+1),]<-(V.quant[, g]*M.quant)/sdev.quant #coeffs beta pour les quantis
      beta[(p1+2):(p1+m+1),1]<-(V.qual[, g]*M.qual) #coeffs beta pour les qualis
    }    
    if(p1>0 & p2==0){
      beta[1, 1] <- -sum(V.quant[, g]* M.quant * gc.quant/sdev.quant)    
      beta[2:(p1+1),]<-(V.quant[, g]*M.quant)/sdev.quant #coeffs beta pour les quantis
    }   
    if(p1==0 & p2>0){
      beta[1, 1] <- -sum(V.qual[, g]* M.qual * gc.qual)
      beta[(p1+2):(p1+m+1),1]<-(V.qual[, g]*M.qual) #coeffs beta pour les qualis
    }
    rownames(beta) <- c("const", colnames(W))
    coef[[g]] <- beta
  }
  
  Z <- rec$Z
  res <- list(call = cl, eig = eig, scores.stand = as.matrix(U), 
              scores = F, V = as.matrix(V), sqload = sqload, A = as.matrix(A), 
              categ.coord = A2, quanti.cor = A1, quali.eta2 = C, rec = rec, 
              ndim = ndim, W = W, ind = result.ind, quanti = result.quanti, 
              levels = result.levels,quali=result.quali, coef = coef, Z = Z, M = Met.global)
  class(res) <- "PCAmix"
  if (graph==TRUE) {
    plot.PCAmix(res)
    if (p1 != p) {
      dev.new()
      plot.PCAmix(res, choice = "levels")
    }
    if (p1 != 0) {
      dev.new()
      plot.PCAmix(res, choice = "cor")
    }
    dev.new()
    plot.PCAmix(res, choice = "sqload")
  }
  return(res)
}