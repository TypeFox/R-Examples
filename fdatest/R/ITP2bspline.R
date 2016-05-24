ITP2bspline <-
function(data1,data2,mu=0,order=2,nknots=dim(data1)[2],B=10000,paired=FALSE){
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
  bspl.basis <- create.bspline.basis(c(1,J),norder=order,breaks=seq(1,J,length.out=nknots))
  ascissa <- seq(1,J,1)
  
  data.fd <- Data2fd(t(eval),ascissa,bspl.basis)
  coeff <- t(data.fd$coef)
  p <- dim(coeff)[2]
  
  #functional data
  npt <- 1000
  ascissa.2 <- seq(1,J,length.out=npt)
  bspl.eval.smooth <- eval.basis(ascissa.2,bspl.basis)
  data.eval <- t(bspl.eval.smooth %*% t(coeff))
  data.eval[(n1+1):n,] <- data.eval[(n1+1):n,] + matrix(data=mu,nrow=n2,ncol=npt)
  
  print('Second step: joint univariate tests')
  #univariate permutations
  T0 <- abs(colMeans(coeff[1:n1,]) - colMeans(coeff[(n1+1):n,])) #sample mean difference
  T_coeff <- matrix(ncol=p,nrow=B)
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
    T_coeff[perm,] <- abs(colMeans(coeff_perm[1:n1,]) - colMeans(coeff_perm[(n1+1):n,]))
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
  ITPresult <- list(basis='B-spline',test='2pop',mu=mu,paired=as.character(paired),coeff=coeff,pval=pval,pval.matrix=matrice_pval_asymm,corrected.pval=corrected.pval,labels=etichetta_ord,data.eval=data.eval,heatmap.matrix=matrice_pval_symm)
  class(ITPresult) = 'ITP2'
  return(ITPresult)
}
