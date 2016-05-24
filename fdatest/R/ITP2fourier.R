ITP2fourier <-
function(data1,data2,mu=0,maxfrequency=floor(dim(data1)[2]/2),B=10000,paired=FALSE){
  fisher_cf_L <- function(L){ #fisher on rows of the matrix L
    return(-2*rowSums(log(L)))
  }
  fisher_cf <- function(lambda){ #fisher on vector lambda
    return(-2*sum(log(lambda)))
  }
  calcola_hotelling_2pop <- function(pop1,pop2,delta.0){ #calcola T^2 pooled per due pop
    mean1 <- colMeans(pop1)
    mean2 <- colMeans(pop2)
    cov1 <- cov(pop1)
    cov2 <- cov(pop2)
    n1 <- dim(pop1)[1]
    n2 <- dim(pop2)[1]
    p <- dim(pop1)[2]
    Sp      <- ((n1-1)*cov1 + (n2-1)*cov2)/(n1+n2-2)
    Spinv   <- solve(Sp)
    T2 <- n1*n2/(n1+n2) * (mean1-mean2-delta.0) %*% Spinv %*% (mean1-mean2-delta.0)
    return(T2)
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
  
  data1 <- as.matrix(data1)
  data2 <- as.matrix(data2)
  
  n1 <- dim(data1)[1]
  n2 <- dim(data2)[1]
  J <- dim(data1)[2]
  n <- n1+n2
  etichetta_ord <- c(rep(1,n1),rep(2,n2))
  data2 <- data2 - matrix(data=mu,nrow=n2,ncol=J)
  
  print('First step: basis expansion')
  #splines coefficients:
  eval <- rbind(data1,data2)
  ak_hat <- NULL
  bk_hat <- NULL
  Period <- J
  
  for(unit in 1:n){ 
    #indice <- 1
    data_temp <- eval[unit,]
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
  data.eval[(n1+1):n,] <- data.eval[(n1+1):n,] + matrix(data=mu,nrow=n2,ncol=npt)
  
  print('Second step: joint univariate tests')
  #univariate permutations
  T0 <- numeric(dim)
  for(freq in 2:dim){
    T0[freq] <- calcola_hotelling_2pop(cbind(ak[1:n1,freq-1],bk[1:n1,freq-1]),cbind(ak[(n1+1):n,freq-1],bk[(n1+1):n,freq-1]),c(0,0))
  }
  T0[1] <- abs(mean(coeff[1:n1,1]) - mean(coeff[(n1+1):n,1]))
  
  T_hotelling <- matrix(nrow=B,ncol=dim)
  for (perm in 1:B){
    if(paired==TRUE){
      if.perm <- rbinom(n1,1,0.5) 
      coeff_perm <- coeff
      for(couple in 1:n1){
        if(if.perm[couple]==1){
          coeff_perm[c(couple,n1+couple),] <- coeff[c(n1+couple,couple),]
        }
      }
    }else if(paired==FALSE){
      permutazioni <- sample(n)
      coeff_perm <- coeff[permutazioni,]
    }
    ak_perm <- coeff_perm[,2:((p+1)/2)]
    bk_perm <- coeff_perm[,((p+1)/2+1):p]
    T_hotelling_temp <- numeric(dim)
    for(freq in 2:dim){
      T_hotelling_temp[freq] <- calcola_hotelling_2pop(cbind(ak_perm[1:n1,freq-1],bk_perm[1:n1,freq-1]),cbind(ak_perm[(n1+1):n,freq-1],bk_perm[(n1+1):n,freq-1]),c(0,0))
    }
    T_hotelling_temp[1] <- abs(mean(coeff_perm[1:n1,1]) - mean(coeff_perm[(n1+1):n,1]))
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
  ITP.result <- list(basis='Fourier',test='2pop',mu=mu,paired=as.character(paired),coeff=coeff,pval=pval,pval.matrix=matrice_pval_asymm,corrected.pval=corrected.pval,labels=etichetta_ord,data.eval=data.eval,heatmap.matrix=matrice_pval_symm)
  class(ITP.result) = 'ITP2'
  return(ITP.result)
}
