ITP1bspline <-
function(data,mu=0,order=2,nknots=dim(data)[2],B=10000){
  fisher_cf_L <- function(L){ #fisher on rows of the matrix L
    return(-2*rowSums(log(L)))
  }
  fisher_cf <- function(lambda){ #fisher on vector lambda
    return(-2*sum(log(lambda)))
  }
  pval.correct <- function(pval.matrix){
    matrice_pval_2_2x <- cbind(pval.matrix,pval.matrix)
    p <- dim(pval.matrix)[2]
    matrice_pval_2_2x <- matrice_pval_2_2x[,(2*p):1]
    corrected.pval <- numeric(p)
    for(var in 1:p){
      pval_var <- matrice_pval_2_2x[p,var]
      inizio <- var
      fine <- var #inizio fisso, fine aumenta salendo nelle righe
      for(riga in (p-1):1){
        fine <- fine + 1
        pval_cono <- matrice_pval_2_2x[riga,inizio:fine]
        pval_var <- max(pval_var,pval_cono)
      }
      corrected.pval[var] <- pval_var
    }
    corrected.pval <- corrected.pval[p:1]
    return(corrected.pval)
  }
  data <- as.matrix(data)
  
  n <- dim(data)[1]
  J <- dim(data)[2]
  labels <- rep(1,n)
  data <- data - matrix(data=mu,nrow=n,ncol=J,byrow=TRUE)
  
  print('First step: basis expansion')
  #splines coefficients:
  bspl.basis <- create.bspline.basis(c(1,J),norder=order,breaks=seq(1,J,length.out=nknots))
  ascissa <- seq(1,J,1)
  
  data.fd <- Data2fd(t(data),ascissa,bspl.basis)
  coeff <- t(data.fd$coef)
  
  p <- dim(coeff)[2]
  
  #functional data
  npt <- 1000
  ascissa.2 <- seq(1,J,length.out=npt)
  bspl.eval.smooth <- eval.basis(ascissa.2,bspl.basis)
  data.eval <- t(bspl.eval.smooth %*% t(coeff) )
  data.eval <- data.eval + matrix(data=mu,nrow=n,ncol=npt)
  
  #univariate permutations
  print('Second step: joint univariate tests')
  T0 <- abs(colMeans(coeff))  #sample mean
  T_coeff <- matrix(ncol=p,nrow=B)
  for (perm in 1:B){
    signs <- rbinom(n,1,0.5)*2 - 1
    coeff_perm <- coeff*signs
    T_coeff[perm,] <- abs(colMeans(coeff_perm))
  }
  pval <- numeric(p)
  for(i in 1:p){
    pval[i] <- sum(T_coeff[,i]>=T0[i])/B
  }
  
  #combination
  print('Third step: interval-wise combination and correction')
  q <- numeric(B)
  L <- matrix(nrow=B,ncol=p)
  for(j in 1:p){
    ordine <- sort.int(T_coeff[,j],index.return=T)$ix
    q[ordine] <- (B:1)/(B)
    L[,j] <- q
  }
  
  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=p,ncol=p)
  matrice_pval_asymm[p,] <- pval[1:p]
  pval_2x <- c(pval,pval)
  L_2x <- cbind(L,L)
  for(i in (p-1):1){
    for(j in 1:p){
      inf <- j
      sup <- (p-i)+j
      T0_temp <- fisher_cf(pval_2x[inf:sup])
      T_temp <- fisher_cf_L(L_2x[,inf:sup])
      pval_temp <- sum(T_temp>=T0_temp)/B
      matrice_pval_asymm[i,j] <- pval_temp
    }
    print(paste('creating the p-value matrix: end of row ',as.character(p-i+1),' out of ',as.character(p),sep=''))
  }
  
  #symmetric combination matrix
  matrice_pval_symm <- matrix(nrow=p,ncol=4*p)
  for(i in 0:(p-1)){
    for(j in 1:(2*p)){
      matrice_pval_symm[p-i,j+i+p] <- matrice_pval_asymm[p-i,(j+1)%/%2]
      if(j+i>2*p-i){
        matrice_pval_symm[p-i,j+i-p] <- matrice_pval_asymm[p-i,(j+1)%/%2]
      }
    }
  }
  
  corrected.pval <- pval.correct(matrice_pval_asymm)
  print('Interval Testing Procedure completed')
  ITP.result <- list(basis='B-spline',test='1pop',mu=mu,coeff=coeff,pval=pval,pval.matrix=matrice_pval_asymm,corrected.pval=corrected.pval,labels=labels,data.eval=data.eval,heatmap.matrix=matrice_pval_symm)
  class(ITP.result) = 'ITP1'
  return(ITP.result)
}
