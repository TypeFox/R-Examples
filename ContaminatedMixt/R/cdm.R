.ncovpar <- function (modelname = NULL, p = NULL, G = NULL) {
  if (is.null(p)) 
    stop("p is null")
  if (is.null(G)) 
    stop("G is null")
  if (is.null(modelname)) 
    stop("modelname is null")
  if (modelname == "EII") 
    npar = 1
  else if (modelname == "VII") 
    npar = G
  else if (modelname == "EEI") 
    npar = p
  else if (modelname == "VEI") 
    npar = p + G - 1
  else if (modelname == "EVI") 
    npar = p * G - G + 1
  else if (modelname == "VVI") 
    npar = p * G
  else if (modelname == "EEE") 
    npar = p * (p + 1)/2
  else if (modelname == "EEV") 
    npar = G * p * (p + 1)/2 - (G - 1) * p
  else if (modelname == "VEV") 
    npar = G * p * (p + 1)/2 - (G - 1) * (p - 1)
  else if (modelname == "VVV") 
    npar = G * p * (p + 1)/2
  else if (modelname == "EVE") 
    npar = p * (p + 1)/2 + (G - 1) * (p - 1)
  else if (modelname == "VVE") 
    npar = p * (p + 1)/2 + (G - 1) * p
  else if (modelname == "VEE") 
    npar = p * (p + 1)/2 + (G - 1)
  else if (modelname == "EVV") 
    npar = G * p * (p + 1)/2 - (G - 1)
  else stop("modelname is not correctly defined")
  return(npar)
}

model.type <- function(modelname=NULL, Sk=NULL, ng=NULL, D=NULL, mtol=1e-10, mmax=10) {
  if (is.null(modelname)) stop("modelname is null")
  
  if (modelname == "EII") val = msEII(Sk=Sk, ng=ng)
  else if (modelname == "VII") val = msVII(Sk=Sk, ng=ng)
  else if (modelname == "EEI") val = msEEI(Sk=Sk, ng=ng)
  else if (modelname == "VEI") val = msVEI(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
  else if (modelname == "EVI") val = msEVI(Sk=Sk, ng=ng)
  else if (modelname == "VVI") val = msVVI(Sk=Sk, ng=ng)
  else if (modelname == "EEE") val = msEEE(Sk=Sk, ng=ng)
  else if (modelname == "EEV") val = msEEV(Sk=Sk, ng=ng)
  else if (modelname == "VEV") val = msVEV(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
  else if (modelname == "VVV") val = msVVV(Sk=Sk, ng=ng)
  else if (modelname == "EVE") val = msEVE(Sk=Sk, ng=ng, D0=D, eplison= mtol, max.iter= mmax)
  else if (modelname == "VVE") val = msVVE(Sk=Sk, ng=ng, D0=D, eplison= mtol, max.iter= mmax)
  else if (modelname == "VEE") val = msVEE(Sk=Sk, ng=ng, eplison= mtol, max.iter= mmax)
  else if (modelname == "EVV") val = msEVV(Sk=Sk, ng=ng)
  else stop("modelname or covtype is not correctly defined")
  
  if (!is.list(val)) val = list(sigma=val)
  return(val)		
}

msEEE <- function(Sk=NULL, ng=NULL) {
  # Sk is an array of with dim (p x p x G)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  W = sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)
  
  val = array(0, c(d,d,G))
  for (g in 1:G) val[,,g] = W
  
  logdetW = log(det(W))
  invW    = solve(W)
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = W
    val$invSigma[,,g] = invW
    val$logdet[g]     = logdetW
  }
  return(val)
}

sumSk.wt <- function(Sk=NULL, wt=NULL, d=NULL, G=NULL) {
  # Sum Sk over the groups with weights wt.
  W = matrix(0, nrow=d, ncol=d )
  for (g in 1:G) W = W + Sk[,,g]* wt[g]
  return(W)
}

msEEV <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  EWk = Sk
  A  = matrix(0,d,d)
  for (g in 1:G) {
    Wk = Sk[,,g]*ng[g]
    EWk[,,g] = eigen(Wk)$vectors
    A = A + t(EWk[,,g]) %*% Wk %*% EWk[,,g]
  }
  
  lam = prod(diag(A))^(1/d)
  A   = A/lam
  lam = lam/sum(ng)
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    =   lam *( EWk[,,g] %*% A %*% t(EWk[,,g]) )
    val$invSigma[,,g] = 1/lam *( EWk[,,g] %*% diag(1/diag(A),d) %*% t(EWk[,,g]) )
    val$logdet[g]     =  d*log(lam) 
  }
  return(val)
}

getA <- function(Ok=NULL, lam=NULL, G=NULL, d=NULL) {
  A  = matrix(0,d,d)
  for (g in 1:G) A = A + Ok[,,g]/lam[g]
  A= diag(A)
  A= diag( A/prod(A)^(1/d) )
  return( A )	
}

getOk <- function(Sk=NULL, ng=NULL, G=NULL) {
  Ok  = Sk
  for (g in 1:G) {
    Wk  = Sk[,,g]*ng[g]
    EWk = eigen(Wk)$vectors
    Ok[,,g] = t(EWk) %*% Wk %*% EWk
    Ok[,,g] = diag(diag(Ok[,,g]))
  }
  return(Ok)	
}

getEkOk <- function(Sk=NULL, ng=NULL, G=NULL) {
  Ok = Sk
  EWk = Sk
  for (g in 1:G) {
    Wk  = Sk[,,g]*ng[g]
    EWk[,,g] = eigen(Wk)$vectors
    Ok[,,g] = t(EWk[,,g]) %*% Wk %*% EWk[,,g]
  }
  return(list(Ok=Ok,EWk=EWk))	
}

msVEV <- function(Sk=NULL, ng=NULL, eplison=1e-14, max.iter= 100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  temp = getEkOk(Sk=Sk, ng=ng, G=G)
  Ok  = temp$Ok
  EWk = temp$EWk
  lam = apply(Ok,3,function(z) { sum(diag(z)) } )/(ng*d)
  A   = getA(Ok=Ok, lam=lam, d=d, G=G) 
  lam = apply(Ok,3,function(z, invA) { sum(diag(z*invA)) }, invA=diag(1/diag(A)) )/(ng*d)
  
  conv = c( d*sum(ng*(1+log(lam))), Inf  )
  count = 1 
  while ( diff(conv)/conv[1] > eplison & count < max.iter) {
    A   = getA(Ok=Ok, lam=lam, d=d, G=G) 
    lam = apply(Ok,3,function(z, invA) { sum(diag(z*invA)) }, invA=diag(1/diag(A)) )/(ng*d)
    conv = c(d*sum(ng*(1+log(lam))), conv[1] )
    count = count +1
  }
  
  val = array(0, c(d,d,G))
  for (g in 1:G) val[,,g] = lam[g] * ( EWk[,,g] %*% A %*% t(EWk[,,g]) )
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    =    lam[g] * ( EWk[,,g] %*% A %*% t(EWk[,,g]) )
    val$invSigma[,,g] =  1/lam[g] * ( EWk[,,g] %*% diag(1/diag(A),d) %*% t(EWk[,,g]) )
    val$logdet[g]     =  d*log(lam[g])
  }
  return(val)
}

msVVV <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    =  Sk[,,g]
    val$invSigma[,,g] =  solve(Sk[,,g])
    val$logdet[g]     =  log(det(Sk[,,g]) )
  }
  return(val)
}

msEEI <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  W = sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)
  B = diag(diag(W))
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G) )
  for (g in 1:G) { 
    val$sigma[,,g]    = B
    val$invSigma[,,g] = diag( 1/diag(B),d )
    val$logdet[g] =  sum(log( diag(B) ))
  }
  return(val)
}

msVEI <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  lam = apply(Sk,3,function(z) { sum(diag(z)) } )/d
  #	W = sumSk.wt(Sk=Sk, wt=ng*lam, d=d, G=G)
  W = sumSk.wt(Sk=Sk, wt=ng/lam, d=d, G=G)
  W = diag(W)
  B = diag(W/prod(W)^(1/d))
  lam = apply(Sk,3,function(z, invB) { sum(diag(z*invB)) }, invB=diag(1/diag(B) ) )/d
  
  conv = c( d*sum(ng*(1+log(lam))), Inf  )
  count = 1 
  while ( abs(diff(conv)) > eplison & count < max.iter) {
    #	while ( count < max.iter) {
    #		W = sumSk.wt(Sk=Sk, wt=ng*lam, d=d, G=G)
    W = sumSk.wt(Sk=Sk, wt=ng/lam, d=d, G=G)
    
    W = diag(W)
    B = diag(W/prod(W)^(1/d))
    lam = apply(Sk,3,function(z, invB) { sum(diag(z*invB)) }, invB=diag(1/diag(B) ) )/d
    
    conv = c(d*sum(ng*(1+log(lam))), conv[1] )
    count = count +1
  }
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam[g]*B
    val$invSigma[,,g] = diag( 1/diag(B) * 1/lam[g],d )
    val$logdet[g]     = d*log(lam[g])
  }
  return(val)
}

msEVI <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  Bk = matrix(0,nrow=d, ncol=G)
  for (g in 1:G) Bk[,g] = diag(Sk[,,g]*ng[g])
  lam = apply(Bk, 2, prod)^(1/d) 
  Bk  = sweep(Bk, 2, 1/lam, FUN="*")	
  lam = sum(lam)/sum(ng)
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam * diag(Bk[,g],d)
    val$invSigma[,,g] = 1/lam * diag(1/Bk[,g],d)
    val$logdet[g]     = d*log(lam)
  }
  return(val)
}

msVVI <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  #	Bk = matrix(0,nrow=d, ncol=G)
  #	for (g in 1:G) Bk[,g] = diag(Sk[,,g]*ng[g])
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = diag(diag(Sk[,,g]),d)
    val$invSigma[,,g] = diag(1/diag(Sk[,,g]),d)
    val$logdet[g]     = sum( log(diag(Sk[,,g])) )
  }
  return(val)
}

msEII <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  W = sumSk.wt(Sk=Sk, wt=ng, d=d, G=G)/sum(ng)
  lam = sum(diag(W))/(sum(ng)*d)
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G) )
  for (g in 1:G) { 
    val$sigma[,,g] = diag(rep(lam,d), d)
    val$invSigma[,,g] = diag(rep(1/lam,d),d)
    val$logdet[g] = d*log(lam) 
  }
  return(val)
}

msVII <- function(Sk=NULL, ng=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    sumdiagSkg = sum(diag(Sk[,,g]))
    val$sigma[,,g] = diag(rep(sumdiagSkg/d,d))
    val$invSigma[,,g] = diag(rep( d/sumdiagSkg, d))
    val$logdet[g] =  d* log( sumdiagSkg ) - d*log(d)
  }
  return(val)
}

msVVE <- function(Sk=NULL, ng=NULL) {
  # Sk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  #	lam = numeric(G)
  #	for (g in 1:G) 	lam[g] = sum(diag(Sk[,,g]))/d
  
  val = array(0, c(d,d,G))
  for (g in 1:G) val[,,g] = diag(rep(sum(diag(Sk[,,g]))/d,d))
  return(val)
}

msEVE <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 10, D0=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1]; G = dim(Sk)[3];

  Wk = Sk
  W  = matrix(0,d,d)
  for (g in 1:G) {
    Wk[,,g] = Sk[,,g]*ng[g]
    W  = W + Wk[,,g]
  }
  #	D  = diag(rep(1,d))
  D = D0
  if (is.null(D)) D  = diag(rep(1,d))
  #	D  = t(eigen(W)$vectors)
  Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% (D) ) }, D=D )
  Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
  #print(c(0, 1, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
  D = newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d)
  #print(c(0, 2, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
  
  conv = c( testval(Wk=Wk,Ak=Ak,D=D,G=G), Inf  )
  count = 1  
  while ( diff(conv)/abs(conv[1]) > eplison & count < max.iter) {
    D = newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d) 
    Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% (D) ) }, D=D )
    Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
    
    #print(c(count, 0, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
    conv = c(testval(Wk=Wk,Ak=Ak,D=D,G=G), conv[1] )
    count = count +1
  }
  lam =  0
  for (g in 1:G) lam = lam  + sum(diag( D %*% diag(1/Ak[,g])%*% t(D) %*% Wk[,,g] ))
  lam = lam/(sum(ng)*d)
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G), D=D  )
  for (g in 1:G) { 
    val$sigma[,,g]    = D %*% diag(lam*Ak[,g])%*% t(D)
    val$invSigma[,,g] = D %*% diag(1/lam*1/Ak[,g])%*% t(D)
    val$logdet[g]     = d*log( lam ) 
  }
  return(val)
  
}

msVVE <- function(Sk=NULL, ng=NULL, eplison=1e-20, max.iter= 10, D0=NULL) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1]; G = dim(Sk)[3];
  Wk = Sk
  W  = matrix(0,d,d)
  for (g in 1:G) {
    Wk[,,g] = Sk[,,g]*ng[g]
    W  = W + Wk[,,g]
  }
  
  D = D0
  if (is.null(D)) D  = diag(rep(1,d))
  
  Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% D ) }, D=D )
  Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
  
  #print(c(0, 1, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
  D = newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d)
  #print(c(0, 2, testval(Wk=Wk,Ak=Ak,D=D,G=G)))		
  
  conv = c( testval(Wk=Wk,Ak=Ak,D=D,G=G), Inf  )
  count = 1  
  while ( diff(conv)/abs(conv[1]) > eplison & count < max.iter) {
    Ak = apply(Wk,3, function(z,D) { diag( t(D) %*% z %*% D ) }, D=D )
    Ak = apply(Ak,2,function(z) { z/prod(z)^(1/length(z))})
    D = newD(D=D, Wk=Wk, Ak=Ak, G=G, d=d) 
    
    conv = c(testval(Wk=Wk,Ak=Ak,D=D,G=G), conv[1] )
    count = count +1
  }
  #print(count)
  lam = numeric(G) 
  for (g in 1:G) lam[g] =sum(diag( D %*% diag(1/Ak[,g])%*% t(D) %*% Sk[,,g] ))/d
  #print(apply(Ak,2,prod))
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = (D %*% diag(lam[g]*Ak[,g])%*% t(D))
    val$invSigma[,,g] = (D %*% diag(1/lam[g]*1/Ak[,g])%*% t(D))
    val$logdet[g]     =  d*log(lam[g])
  }
  return(val)
}

msVEE <- function(Sk=NULL, ng=NULL, eplison=1e-14, max.iter=100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  Wk = Sk
  W  = matrix(0,d,d)
  for (g in 1:G) {
    Wk[,,g] = Sk[,,g]*ng[g]
    W  = W + Wk[,,g]
  }
  
  C = W/det(W)^(1/d)
  invC = solve(C)
  lam = apply(Sk,3,function(z, invC) { sum(diag(z %*% invC)) }, invC=invC)/d
  
  val1 = sum(apply(Wk,3,function(z, invC) { sum(diag(z*invC)) }, invC= invC)/lam) +d*sum(ng*lam)
  conv = c(val1, Inf  )
  count = 1 
  while ( diff(conv)/conv[1]  > eplison & count < max.iter) {
    #for (i in 1:max.iter) {
    C = sumSk.wt(Sk=Wk, wt=1/lam, d=d,G=G)
    C = C/det(C)^(1/d)		
    invC = solve(C)
    
    lam = apply(Sk,3,function(z, invC) { sum(diag(z %*% invC)) }, invC=invC)/d
    val1 = sum(apply(Wk,3,function(z, invC) { sum(diag(z*invC)) }, invC= invC )/lam) +d*sum(ng*lam)
    conv = c(val1, conv[1] )
    count = count +1
  }
  #print(c(count,det(C), lam)	)
  
  invC = solve(C)
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam[g]*C
    val$invSigma[,,g] = 1/lam[g]*invC
    val$logdet[g]     = d*log(lam[g])
  }
  return(val)
}

msEVV <- function(Sk=NULL, ng=NULL, eplison=1e-12, max.iter= 100) {
  # Wk is a list of length G with matrix (p x p)
  # ng is a vector of length G of weights of Wk
  d = dim(Sk)[1];
  G = dim(Sk)[3];
  
  Wk = Sk
  Ck = Sk
  lam = numeric(G)
  for (g in 1:G) {
    Wk[,,g] = Sk[,,g]*ng[g]
    lam[g]    = det(Wk[,,g])^(1/d)
  }
  Ck = sweep(Wk, 3, 1/lam, FUN="*")	
  lam = sum(lam)
  
  
  val = list(sigma=array(0, c(d,d,G)), invSigma=array(0, c(d,d,G)), logdet=numeric(G)  )
  for (g in 1:G) { 
    val$sigma[,,g]    = lam * Ck[,,g]
    val$invSigma[,,g] = 1/lam * solve(Ck[,,g])
    val$logdet[g]     =  d* log(lam ) 
  }
  return(val)
}

newD3.MM <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
  z = matrix(0,d,d)
  lambda =0 
  for (g in 1:G) {
    lambdak = max(eigen(Wk[,,g])$values)
    z = z + diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]  -lambdak *(diag(1/Ak[,g])%*% t(D) )
  } 
  z1 = svd(z)
  Xk1 = (z1$v) %*% t(z1$u) 
  return( Xk1 )
}

newD4.MM <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
  z = matrix(0,d,d)
  lambda =0 
  for (g in 1:G) {
    lambdak = max(1/Ak[,g])
    #z = z + diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]  -lambdak *(diag(1/Ak[,g])%*% t(D) )
    z = z + Wk[,,g] %*% (D) %*% diag(1/Ak[,g])  -lambdak *(Wk[,,g]%*% (D) )
  } 
  z1 = svd(z)
  #	Xk1 = (z1$v) %*% t(z1$u) 
  #	return( t(Xk1) )
  # OR 
  Xk1 = (z1$v) %*% t(z1$u) 
  return( t(Xk1) )
  
}

newD <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL, tmax = 100) {
  D6 =D
  D6 = newD3.MM(D=D6, d=d, G=G, Wk=Wk, Ak=Ak, tmax = 100)
  D6 = newD4.MM(D=D6, d=d, G=G, Wk=Wk, Ak=Ak, tmax = 100)
  return(D6)
}

testval <- function(Wk=NULL, Ak=NULL, D=NULL, G=NULL) {
  z = numeric(G)
  #	for (g in 1:G) z[g] = sum(diag( D %*%  diag(1/Ak[,g]) %*% t(D) %*% Wk[,,g]))
  for (g in 1:G) z[g] = sum(diag( t(D) %*% Wk[,,g] %*% D %*%  diag(1/Ak[,g])  ))
  return(sum(z))
}

testgrad.D <- function(D=NULL, d=NULL, G=NULL, Wk=NULL, Ak=NULL) {
  z = matrix(0,d,d)
  for (g in 1:G) z = z + Wk[,,g] %*% D %*% diag(1/Ak[,g]) 
  return(2*z)
}