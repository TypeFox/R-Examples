# Genomic Ensamble Machine Learning Algorithm
# Ensembles Bag+Ridge+SSVS+RF+KNN+RKHS
# Algorithm steps
# 1. Eigendecompose a kinship (Gaussian kernel)
# 2. Sample X observations (Bootstrapping)
# 3. Sample U Eigenvectors with Pr(Eigenvalue) for clustering X and unobserved via KNN in K clusters
# 4. Sample V Eigenvector of X with Pr(R2(Eig,Y))
# 5. Sample S marker of X with Pr(R2(snp,Y))
# 6. Fit the model y = mu + V + K:S, rigding V+K:S
# 7. Predict unobserved values
# 8. Repeat 2-7 several times
# 9. Aggregate empirical distribution on predicted Ys

GELA = function(y,x,It=500,CL=4,PC=15,Par=20,Bg=0.75,rdg=1e-6){
  mi = which(is.na(y))
  nmi = length(mi)
  cat('Get Gauss. Kernel and its Eigen.\n')
  D2 = as.matrix(dist(x)^2)
  GK = exp(-D2/mean(D2))
  eig = eigen(GK,T)
  Y = y[-mi]
  X = x[-mi,]
  xx = x[mi,]
  E = eig$vectors[-mi,]
  ee = eig$vectors[mi,]
  R2x = cor(Y,X)^2
  R2x = as.vector(R2x/sum(R2x))
  R2e = cor(Y,E)^2
  R2e = as.vector(R2e/sum(R2e))
  R2v = as.vector(eig$values/sum(eig$values))
  Bag = round(Bg*length(Y))
  yhat = matrix(NA,nmi,It)
  Xb = matrix(NA,Bag,Par)
  Eb = matrix(NA,Bag,PC)
  Pb = matrix(NA,Bag+nmi,PC)
  n = length(Y)
  nW = (1+PC+CL*Par)
  Wp = matrix(NA,Bag,nW)
  Wv = matrix(NA,nmi,nW)
  pb = txtProgressBar(style = 3)
  for(i in 1:It){
    o = sample(1:n,Bag)
    t = sample(1:ncol(X),Par,prob=R2x)
    g = sample(1:ncol(E),PC,prob=R2e)
    p = sample(1:ncol(E),PC,prob=R2v)
    Xb[1:Bag,1:Par] = X[o,t]
    Eb[1:Bag,1:PC] = E[o,g]
    Pb[1:Bag,1:PC] = E[o,p]
    Pb[(Bag+1):(Bag+nmi),1:PC] = ee[,p]
    K = factor(kmeans(Pb,CL)$cluster)
    K1 = K[1:Bag]
    K2 = K[(Bag+1):(Bag+nmi)]
    Wp[1:Bag,1:nW] = model.matrix(~Eb+K1:Xb)
    Wv[1:nmi,1:nW] = model.matrix(~ee[,g]+K2:xx[,t])
    WW = crossprod(Wp)
    diag(WW) = diag(WW)+c(0,rep(rdg,nW-1))
    Wy = crossprod(Wp,Y[o])
    fit = solve(WW,Wy)
    yhat[,i] = Wv %*% fit
    setTxtProgressBar(pb, i/It)}
  close(pb)
  Hat = rowMeans(yhat)
  yFit = y
  yFit[mi] = Hat
  return(yFit)
}

wGELA = function(y,x,It=500,CL=4,PC=15,Par=20,Bg=0.75,rdg=1e-6){
  mi = which(is.na(y))
  nmi = length(mi)
  cat('Get Gauss. Kernel and its Eigen.\n')
  D2 = as.matrix(dist(x)^2)
  GK = exp(-D2/mean(D2))
  eig = eigen(GK,T)
  Y = y[-mi]
  X = x[-mi,]
  xx = x[mi,]
  E = eig$vectors[-mi,]
  ee = eig$vectors[mi,]
  R2x = cor(Y,X)^2
  R2x = as.vector(R2x/sum(R2x))
  R2e = cor(Y,E)^2
  R2e = as.vector(R2e/sum(R2e))
  R2v = as.vector(eig$values/sum(eig$values))
  Bag = round(Bg*length(Y))
  yhat = matrix(NA,nmi,It)
  Xb = matrix(NA,Bag,Par)
  Eb = matrix(NA,Bag,PC)
  Pb = matrix(NA,Bag+nmi,PC)
  n = length(Y)
  nW = (1+PC+CL*Par)
  Wp = matrix(NA,Bag,nW)
  Wv = matrix(NA,nmi,nW)
  WR2 = rep(NA,It)
  pb = txtProgressBar(style = 3)
  for(i in 1:It){
    o = sample(1:n,Bag)
    t = sample(1:ncol(X),Par,prob=R2x)
    g = sample(1:ncol(E),PC,prob=R2e)
    p = sample(1:ncol(E),PC,prob=R2v)
    Xb[1:Bag,1:Par] = X[o,t]
    Eb[1:Bag,1:PC] = E[o,g]
    Pb[1:Bag,1:PC] = E[o,p]
    Pb[(Bag+1):(Bag+nmi),1:PC] = ee[,p]
    K = factor(kmeans(Pb,CL)$cluster)
    K1 = K[1:Bag]
    K2 = K[(Bag+1):(Bag+nmi)]
    Wp[1:Bag,1:nW] = model.matrix(~Eb+K1:Xb)
    Wv[1:nmi,1:nW] = model.matrix(~ee[,g]+K2:xx[,t])
    WW = crossprod(Wp)
    diag(WW) = diag(WW)+c(0,rep(rdg,nW-1))
    Wy = crossprod(Wp,Y[o])
    fit = solve(WW,Wy)
    WR2[i] = cor(Wp%*%fit,Y[o])^2
    yhat[,i] = Wv %*% fit
    setTxtProgressBar(pb, i/It)}
  close(pb)
  Hat = apply(yhat,1,weighted.mean,w=WR2)
  yFit = y
  yFit[mi] = Hat
  return(yFit)
}
