ITP1fourier <-
function(data,mu=0,maxfrequency=floor(dim(data)[2]/2),B=10000){
  fisher_cf_L <- function(L){ #fisher on rows of the matrix L
    return(-2*rowSums(log(L)))
  }
  fisher_cf <- function(lambda){ #fisher on vector lambda
    return(-2*sum(log(lambda)))
  }
  calcola_hotelling <- function(mean0,x){
    x_mean <- colMeans(x)
    x_cov <- cov(x)
    x_invcov <- solve(x_cov)
    x_T2 <- n * (x_mean-mean0) %*% x_invcov %*% (x_mean-mean0) 
    return(x_T2)
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
  data <- data - matrix(data=mu,nrow=n,ncol=J)
  labels <- rep(1,n)
  
  print('First step: basis expansion')
  #fourier coefficients:
  ak_hat <- NULL
  bk_hat <- NULL
  
  for(unit in 1:n){ 
    #indice <- 1
    data_temp <- data[unit,]
    Period <- length(data_temp)
    abscissa <- 0:(Period-1)
    trasformata <- fft(data_temp)/length(abscissa)
    ak_hat <- rbind(ak_hat,2*Re(trasformata)[1:(maxfrequency+1)])
    bk_hat <- rbind(bk_hat,-2*Im(trasformata)[2:(maxfrequency+1)]) 
  }
  coeff <- cbind(ak_hat,bk_hat)
  p <- dim(coeff)[2]
  
  a0 <- coeff[,1]
  ak <- coeff[,2:((p+1)/2)]
  bk <- coeff[,((p+1)/2+1):p]
  dim <- (p+1)/2
  
  #functional data
  K <- p
  if(K %% 2 ==0){
    K <- K+1
  }
  npt <- 1000
  ascissa.smooth <- seq(0, Period, length.out=npt)
  basis <- matrix(0,nrow=npt,ncol=K)
  basis[,1] <- 1/2
  for(i in seq(2,(K-1),2)){
    basis[,i] <- sin(2*pi*(i/2)*ascissa.smooth/Period)
  }
  for(i in seq(3,(K),2)){
    basis[,i] <- cos(2*pi*((i-1)/2)*ascissa.smooth/Period)
  }
  basis.ord <- cbind(basis[,seq(1,K,2)],basis[,seq(2,K-1,2)])
  data.eval <- coeff %*% t(basis.ord)
  data.eval <- data.eval + matrix(data=mu,nrow=n,ncol=npt)
  
  #univariate permutations
  print('Second step: joint univariate tests')
  T0 <- numeric(dim)
  for(freq in 2:dim){
    T0[freq] <- calcola_hotelling(c(0,0),cbind(ak[,freq-1],bk[,freq-1]))
  }
  T0[1] <- abs(sum(coeff[,1]))
  
  T_hotelling <- matrix(nrow=B,ncol=dim)
  
  for (perm in 1:B){
    signs <- rbinom(n,1,0.5)*2-1
    coeff_perm <- coeff*signs
    ak_perm <- coeff_perm[,2:((p+1)/2)]
    bk_perm <- coeff_perm[,((p+1)/2+1):p]
    T_hotelling_temp <- numeric(dim)
    for(freq in 2:dim){
      T_hotelling_temp[freq] <- calcola_hotelling(c(0,0),cbind(ak_perm[,freq-1],bk_perm[,freq-1]))
    }
    T_hotelling_temp[1] <- abs(sum(coeff_perm[,1]))
    T_hotelling[perm,] <- T_hotelling_temp
  }
  
  pval <- numeric(dim)
  for(i in 1:dim){
    pval[i] <- sum(T_hotelling[,i]>=T0[i])/B
  }
  
  #combination
  print('Third step: interval-wise combination and correction')
  q <- numeric(B)
  L <- matrix(nrow=B,ncol=dim)
  for(j in 1:dim){
    ordine <- sort.int(T_hotelling[,j],index.return=T)$ix
    q[ordine] <- (B:1)/(B)
    L[,j] <- q
  }
  
  #asymmetric combination matrix:
  matrice_pval_asymm <- matrix(nrow=dim,ncol=dim)
  matrice_pval_asymm[dim,] <- pval[1:(dim)]
  pval_2x <- c(pval,pval)
  L_2x <- cbind(L,L)
  for(i in (dim-1):1){
    for(j in 1:dim){
      inf <- j
      sup <- (dim-i)+j
      T0_temp <- fisher_cf(pval_2x[inf:sup])
      T_temp <- fisher_cf_L(L_2x[,inf:sup])
      pval_temp <- sum(T_temp>=T0_temp)/B
      matrice_pval_asymm[i,j] <- pval_temp
    }
    print(paste('creating the p-value matrix: end of row ',as.character(dim-i+1),' out of ',as.character(dim),sep=''))
  }
  
  #symmetric combination matrix
  matrice_pval_symm <- matrix(nrow=dim,ncol=4*dim)
  for(i in 0:(dim-1)){
    for(j in 1:(2*dim)){
      matrice_pval_symm[dim-i,j+i+dim] <- matrice_pval_asymm[dim-i,(j+1)%/%2]
      if(j+i>2*dim-i){
        matrice_pval_symm[dim-i,j+i-dim] <- matrice_pval_asymm[dim-i,(j+1)%/%2]
      }
    }
  }
  
  corrected.pval <- pval.correct(matrice_pval_asymm)
  print('Interval Testing Procedure completed')
  ITP.result <- list(basis='Fourier',test='1pop',mu=mu,coeff=coeff,pval=pval,pval.matrix=matrice_pval_asymm,corrected.pval=corrected.pval,labels=labels,data.eval=data.eval,heatmap.matrix=matrice_pval_symm)
  class(ITP.result) = 'ITP1'
  return(ITP.result)
}
