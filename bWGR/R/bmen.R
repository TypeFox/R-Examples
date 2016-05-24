# Metropolis Elastic Net
BagMEN = function(y,X,it=500,bi=100,bag=0.5,alpha=0.5,wpe=1){
  # Function to update beta
  upB = function(e,mu,X,b,l,a,xx,p,E2,X2,bag,pi,wpe,O){
    xx = xx*bag
    a_new = rbeta(1,10*pi,10*(1-pi))
    b1 = b2 = rep(NA,p)
    e1 = e2 = e
    L1_1 = (1-a)*l/(2*xx)
    L2_1 = (a)*(xx/(xx*l))
    L1_2 = (1-a_new)*l/(2*xx)
    L2_2 = (1-a_new)*(xx/(xx*l))
    # Update regression coefficients
    for(j in O){
      # Old alpha
      xy1 = (crossprod(e1,X[,j])+b[j]*xx[j])/xx[j]
      s1 = sign(xy1)
      beta = abs(xy1)-L1_1[j]
      ben = s1*beta*L2_1[j]
      b1[j] = ifelse(beta>0,ben,0)
      e1 = e1 - X[,j]*(b1[j]-b[j])
      # New alpha
      xy2 = (crossprod(e2,X[,j])+b[j]*xx[j])/xx[j]
      s2 = sign(xy2)
      beta = abs(xy2)-L1_2[j]
      ben = s2*beta*L2_2[j]
      b2[j] = ifelse(beta>0,ben,0)
      e2 = e2 - X[,j]*(b2[j]-b[j])
    }
    # Loss function
    SSPE_1 = sum(as.vector(tcrossprod(b1,X2)-E2)^2)+0.00001
    SSPE_2 = sum(as.vector(tcrossprod(b2,X2)-E2)^2)+0.00001
    LOSS1 = wpe*SSPE_1+crossprod(e1)+l*(0.5*crossprod(b1)*(1-a)+sum(abs(b1))*a)
    LOSS2 = wpe*SSPE_2+crossprod(e2)+l*(0.5*crossprod(b2)*(1-a_new)+sum(abs(b2))*a_new)
    LR = LOSS2/(LOSS1+LOSS2)
    # METROPOLIS
    if(LR<runif(1)){
      P=list('b'=b2,'a'=a_new,'e'=e2,'oob'=SSPE_2)
    }else{
      P=list('b'=b1,'a'=a,'e'=e1,'oob'=SSPE_2)}
    return(P)
  }
  # Missing
  if(anyNA(y)){
    mis = which(is.na(y))
    Y = y
    XX = X[mis,]
    y = y[-mis]
    X = X[-mis,]
    MISS = TRUE
  }else{
    MISS = FALSE
  }
  # Data
  O = order(cor(X,y)^2,decreasing = TRUE)
  n = nrow(X)
  p = ncol(X)
  bn = round(n*bag)
  xx = apply(X,2,function(x)crossprod(x))
  MC = it-bi
  # Parameters
  mu = mean(y)
  e = y-mu
  b = rep(0,p)
  a = 1
  l = 1
  # Store posterior
  B = rep(0,p)
  MU = A = L = SSPE = 0
  # Loop
  pb = txtProgressBar(style = 3)
  for(i in 1:it){
    
    if(bag==1){
      a2 = rbeta(1,10*alpha,10*(1-alpha))
      e1=e2=e
      UP1 = GSEN(X=X,b=b,O=O,xx=xx,e=e1,L=l,A=a,p=p)
      UP2 = GSEN(X=X,b=b,O=O,xx=xx,e=e1,L=l,A=a2,p=p)
      b1 = UP1$b
      e1 = UP1$e
      b2 = UP2$b
      e2 = UP2$e
      LOSS1 = crossprod(e1)+l*(0.5*crossprod(b1)*(1-a)+sum(abs(b1))*a)
      LOSS2 = crossprod(e2)+l*(0.5*crossprod(b2)*(1-a2)+sum(abs(b2))*a2)
      LR = LOSS2/(LOSS1+LOSS2)
      ru = runif(1)
      if(LR<ru){
        b=b2;a=a2;e=e2
      }else{
        b=b1;a=a;e=e1}
      
    }else{
      # Bagging
      s = sort(sample(1:n,n-bn,replace=FALSE))
      # UPDATE
      UP = upB(e[-s],mu,X[-s,],b,l,a,xx*bag,p,e[s],
               X[s,],bag,alpha,wpe,O=O)
      b = UP[[1]]
      a = UP[[2]]
      e = UP[[3]]
    }
    mu = mu + mean(e)
    S_prior = runif(1)
    df_prior = rpois(1,4)
    ve = crossprod(e+S_prior)/(bn+df_prior)
    vb = crossprod(b+S_prior)/(p+df_prior)
    l = ve/vb
    e = as.vector(y-(mu+tcrossprod(b,X)))
    # STORE
    if(it>bi){
      B = B+b
      MU = MU+mu
      A = A+a
      L = L+l
      if(bag<1) SSPE = SSPE+UP$oob/(n-bn)
    }
    setTxtProgressBar(pb, i/it)
  }
  close(pb)
  # Posterior
  Bhat = B/MC
  MUhat = MU/MC
  Ahat = A/MC
  Lhat = L/MC
  MSPEout = mean(SSPE)/MC
  # Prediction
  if(MISS){
    Yhat = Y
    Yhat[-mis] = as.vector(mu+tcrossprod(Bhat,X))
    Yhat[mis] = as.vector(mu+tcrossprod(Bhat,XX))
  }else{
    Yhat = as.vector(mu+tcrossprod(Bhat,X))
  }
  # OUTPUT
  LIST = list('hat'=Yhat,'coef'=Bhat,'b0'=MUhat,
              'alp'=Ahat,'lmb'=Lhat,'MSPEoob'=MSPEout)
  
}
