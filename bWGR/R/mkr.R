mkr = function(Y,K=NULL,eK=NULL,it=500,bu=200,th=3,
               df=5,R2=0.5,EigT=0.05,verb=FALSE){
  
  if(is.null(K)&is.null(eK)) stop('Either K or eK have to be provided')
  if(is.null(eK)) eK = eigen(K)
  
  # Inputs
  U = eK$vectors
  D = eK$values
  if(!is.null(EigT)) {U=U[,1:sum(D>EigT)];D=D[D>EigT]}
  k = ncol(Y)
  n = nrow(Y)
  p = ncol(U)  
  sk_prior = R2*apply(Y,2,var,na.rm=T)*(df+2)
  s = diag(sk_prior*df,k)
  N = crossprod(!is.na(Y))
  mN = mean(N)
  B = matrix(0,p,k)
  E = apply(Y,2,function(x)x-mean(x,na.rm = T))
  E[is.na(E)]=0
  xx = rep(1,p)
  mu = colMeans(Y,na.rm = T)
  v = var(Y,na.rm = T)
  VA = 0.1*diag(diag(var(Y,na.rm = T)))
  VE = 0.1*diag(diag(var(Y,na.rm = T)))
  
  # Store
  mc = seq(bu,it,th)
  lmc = length(mc)
  mcMu = rep(0,k)
  mcB = B
  mcVA = mcVE = matrix(0,k,k)
  
  # Indicators
  y = z = list()
  for(i in 1:k){
    z[[i]] = which(!is.na(Y[,i]))
    y[[i]] = Y[z[[i]],i]
  }
  
  # MCMC
  if(verb) pb = txtProgressBar(style=3)
  
  for(j in 1:it){
    
    mEA = diag(solve(VA)%*%VE)
    
    # Update regression coefficients
    for(i in 1:k){
      
      L = mEA[i]/D
      up = KMUP(X=U[z[[i]],],b=B[,i],xx=xx,
                E=E[z[[i]],i],L=L,p=p,Ve=VE[i,i],pi=0)
      B[,i] = up[[1]]    
      err = up[[3]]
      mu[i] = rnorm(1,mu[i]+mean(err),VE[i,i]/N[i,i])
      E[z[[i]],i] = y[[i]]-U[z[[i]],]%*%B[,i]-mu[i] 
    }
    
    # Update variance components
    SSa = crossprod(B/D,B)
    VA = solve(rWishart(1,p+df,solve(SSa+s))[,,1])
    SSe = (crossprod(E)/N)*mN
    VE = solve(rWishart(1,mN+2,solve(SSe))[,,1])
    
    if(j%in%mc){
      mcMu = mcMu + mu
      mcB = mcB + B
      mcVA = mcVA + VA
      mcVE = mcVE + VE
    }
    
    if(verb) setTxtProgressBar(pb, j/it)
  }
  
  if(verb) close(pb)
  
  # Posterior means
  mcMu = mcMu/lmc
  mcB = mcB/lmc
  mcVA = mcVA/lmc
  mcVE = mcVE/lmc
  A = U%*%mcB
  
  final = list('BV'=A,'VA'=mcVA,'VE'=mcVE)
  return(final)
}


Hmat = function (ped,gen=NULL,Diss=FALSE,inbred=FALSE){
  
  # ped - Pedigree. Numeric matrix n by 3.
  # gen - Numeric matrix where rownames regards the ped ID
  # Diss - Logical. If true, it ignores inbreeding.
  # inbred - Logical. Use genomic for FIS or set main diagonal as 2.
  
  n = nrow(ped)
  A = diag(0, n)
  
  if(is.null(gen)){
    
    # WITHOUT GENOTYPES
    for (i in 1:n) {
      for (j in 1:n) {
        if (i > j) {
          A[i, j] = A[j, i]
        } else {
          d = ped[j, 2]
          s = ped[j, 3]
          if (d == 0) Aid = 0 else Aid = A[i, d]
          if (s == 0) Ais = 0 else Ais = A[i, s]
          if (d == 0 | s == 0) Asd = 0 else Asd = A[d, s]
          if (i == j) Aij = 1 + 0.5 * Asd else Aij = 0.5 * (Aid + Ais)
          A[i, j] = Aij
        }}}
    
  }else{
    
    # WITH GENOTYPES 
    G = as.numeric(rownames(gen))
    if(any(is.na(gen))){
      ND = function(g1,g2){
        X = abs(g1-g2)
        L = 2*sum(!is.na(X))
        X = sum(X,na.rm=TRUE)/L
        return(X)} 
    }else{
      ND = function(g1,g2) sum(abs(g1-g2))/(2*length(g1))
    }
    Inbr = function(g) 2-mean(g==1)
    
    # LOOP     
    for (i in 1:n) {
      for (j in 1:n) {
        
        #######################
        if (i > j) {
          A[i, j] = A[j, i]
        } else {
          
          d = ped[j, 2]
          s = ped[j, 3]
          
          if(j%in%G){
            
            if(d!=0&d%in%G){ Aid=2*ND(gen[paste(j),],gen[paste(d),]) }else{if(d==0) Aid=0 else Aid=A[i,d]}
            if(s!=0&s%in%G){ Ais=2*ND(gen[paste(j),],gen[paste(s),]) }else{if(s==0) Ais=0 else Ais=A[i,s]}
            Asd = Inbr(gen[paste(j),])
            
          }else{
            if (d == 0) Aid = 0 else Aid = A[i, d]
            if (s == 0) Ais = 0 else Ais = A[i, s]
            if (d == 0 | s == 0) Asd = 0 else Asd = A[d, s]
          }
          
          if (i == j) Aij = ifelse(Diss,1,ifelse(inbred,ifelse(i%in%G,Asd,2),1+0.5*Asd))
          else Aij = 0.5*(Aid+Ais)
          A[i, j] = round(Aij,6)
        }
        #######################
        
      }}
  }
  return(A)
}
