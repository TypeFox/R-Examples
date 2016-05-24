
#########################
#   Hypothesis Testing
#########################

#p.s: the current code is only for a single environmental exposure

GA.prelim<-function(Y,time,X=NULL,corstr="exchangeable"){
  ##Preliminary
  Y<-as.matrix(Y);N<-nrow(Y)
  cluster.id<-unique(time[,1]);m<-length(cluster.id)
  #sort the data by cluster.id
  order.index<-order(match(time[,1],cluster.id))
  if (sum(time[,1]!=time[order.index,1])>0){
    msg<-sprintf("time was not sorted, now data is grouped by subject id")
    warning(msg,call.=F)
    Y<-as.matrix(Y[order.index,])
    if(length(X)!=0){X<-as.matrix(X[order.index])}
    time<-time[order.index,]
  }

  if(length(X)!=0){
    Z.A<-svd(as.matrix(X))$u
    nullgee<-geem(Y~Z.A,family=gaussian,corstr=corstr,id=time[,1])
  }else{
    nullgee<-geem(Y~1,family=gaussian,corstr=corstr,id=time[,1])
  }

  #prepare the intermediate results
  result.prelim<-list(Y=Y,time=time,X=X,cluster.id=cluster.id,m=m,corstr=corstr,Z.A=Z.A,nullgee=nullgee)
  return(result.prelim)
}

#G is an m*q matrix, each row coresponds to one subject.
GA.test<-function(result.prelim,G,Gsub.id=NULL,bootstrap=NULL,MinP.adjust=NULL,impute.method='fixed'){
  ## Load preliminary data
  Y<-result.prelim$Y;corstr<-result.prelim$corstr
  m<-result.prelim$m;X<-result.prelim$X
  time<-result.prelim$time;cluster.id<-result.prelim$cluster.id
  Z.A<-result.prelim$Z.A #Z.A0=svd(cbind(rep(1,N),X, M.E))$u
  nullgee<-result.prelim$nullgee

  ## Deal with the genotype
  SNP.list<-colnames(G)
  # match the phenotype and genotype subject ID
  if(length(Gsub.id)==0){Gsub.id<-cluster.id}
  G<-as.matrix(G[match(cluster.id,Gsub.id),])

  # missing genotype imputation
  G[G==9]<-NA
  N_MISS<-sum(is.na(G))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(G)/ncol(G))
    warning(msg,call.=F)
    G<-Impute(G,impute.method)
  }
  # centered genotype
  G<-as.matrix(G);center.G<-t(t(G)-colMeans(G))
  MAF<-colMeans(as.matrix(G[,colMeans(center.G^2)!=0]))/2;MAF[MAF>0.5]<-1-MAF[MAF>0.5]
  G<-as.matrix(center.G[,colMeans(center.G^2)!=0]);n.G<-ncol(G)
  SNP.name<-SNP.list[colMeans(center.G^2)!=0]
  if(ncol(G) == 0){
    msg<-sprintf("G does not include any heterozygous variant")
    stop(msg)
  }

  # generate long form genotype
  Z<-as.matrix(G[match(time[,1],cluster.id),])

  if(length(bootstrap)==0){
    # calculate score
    GA.fit<-GA.Score(nullgee,Y,time,Z.A,Z,corstr) #Z.A does not include the intercept
    # calculate score statistics and their covariance matrix
    score<-GA.fit$score;M<-GA.fit$M
    # calculate p-value
    TestStat<-t(score)%*%score
    EigenM<-eigen(M,symmetric=TRUE,only.values=TRUE)$values
    p.value<-davies(TestStat,EigenM)$Qq

    # implement a single SNP based analysis by GEE score tests
    p.single<-cbind(MAF,as.vector(pchisq(as.vector(score)^2/diag(M),df=1,lower.tail=F)))
    colnames(p.single)<-c('MAF','p.value')
    rownames(p.single)<-SNP.name
  }else{
    # calculate score
    GA.fit<-GA.bootstrap(nullgee,Y,time,Z.A,Z,corstr,bootstrap) #Z.A does not include the intercept
    score<-GA.fit$score;perturb.score<-GA.fit$perturb.score
    Q<-t(score)%*%score
    # calculate bootstrap based mean, variance and kurtosis
    perturb.Q<-apply(perturb.score^2,2,sum)
    perturb.mean<-mean(perturb.Q)
    perturb.variance<-var(perturb.Q)
    perturb.kurtosis<-mean((perturb.Q-perturb.mean)^4)/perturb.variance^2-3
    # calculate p.value
    if (perturb.kurtosis>0){df<-12/perturb.kurtosis}else{
      df<-100000
    }
    p.value<-1-pchisq((Q-perturb.mean)*sqrt(2*df)/sqrt(perturb.variance)+df,df)
    # implement a single SNP based analysis by GEE score tests
    perturb.variance.single<-apply(perturb.score,1,var)
    p.single<-cbind(MAF,as.vector(pchisq(as.vector(score)^2/perturb.variance.single,df=1,lower.tail=F)))
    colnames(p.single)<-c('MAF','p.value')
    rownames(p.single)<-SNP.name
  }
  if(length(MinP.adjust)!=0){
    eigen.G<-eigen(cor(G),only.values=TRUE)$values
    me.G<-effect.n(eigen.G,MinP.adjust)
    p.MinP<-min(min(p.single[,2])*me.G,1)
  }else{p.MinP<-NA}

  return(list(n.marker=n.G,p.single=p.single,p.MinP=p.MinP,p.value=p.value))
}

#################################################
#   Calculate Score Statistics and their Variance
#################################################

GA.bootstrap<-function(nullgee,Y,time,Z.A,Z.I,corstr,bootstrap) #Z.A does not include the intercept
{
  rho<-as.numeric(nullgee$alpha);N<-nrow(time);Z.A<-cbind(rep(1,N),Z.A) #add the intercept
  mu<-colSums(nullgee$beta*t(Z.A));Y.res<-Y-mu
  #rho<-nullgee$geese$alpha;mu<-nullgee$fitted.values;Y.res<-Y-mu
  cluster.id<-unique(time[,1]);m<-length(cluster.id);dimI<-ncol(Z.I)
  ##Generalized score test
  score.array<-matrix(0,dimI,m)
  n.total<-1;n.rep<-as.numeric(table(time[,1]))
  for (i in 1:m)
  {
    ni<-n.rep[i]
    index<-n.total:(n.total+ni-1)
    n.total<-n.total + ni
    #working covariance matrix
    if (corstr=="exchangeable"){Vi<-diag(1-rho,ni)+ matrix(rho, ni, ni);Vi.inv<-solve(Vi)}
    if (corstr=="ar1"){Vi<-rho^abs(outer(time[index,2], time[index,2],"-"));Vi.inv<-solve(Vi)}
    if (corstr=="independence"){Vi<-diag(1,ni);Vi.inv<-Vi}
    #if (link=="logit"){Ci<-mu[index]*(1-mu[index]);Vi.inv<-outer(sqrt(Ci),sqrt(Ci)^-1,"*")*Vi.inv}
    #score matrices
    Z.Ii<-Z.I[index,];if (ni==1) {Z.Ii<-t(Z.Ii)}
    score.array[,i]<-t(Z.Ii)%*%(Vi.inv%*%Y.res[index])/sqrt(m)
  }
  score<-as.matrix(apply(score.array,1,sum))
  # null distribution
  perturb.coef<-matrix(rbinom(m*bootstrap,1,0.5)*2-1,m,bootstrap)
  perturb.score<-score.array%*%perturb.coef

  return(list(score=score,perturb.score=perturb.score))#Q: test statistic; M: covariance matrix
}

GA.Score<-function(nullgee,Y,time,Z.A,Z.I,corstr) #Z.A does not include the intercept
{
  rho<-as.numeric(nullgee$alpha);N<-nrow(time);Z.A<-cbind(rep(1,N),Z.A) #add the intercept
  mu<-colSums(nullgee$beta*t(Z.A));Y.res<-Y-mu
  #rho<-nullgee$geese$alpha;mu<-nullgee$fitted.values;Y.res<-Y-mu
  cluster.id<-unique(time[,1]);m<-length(cluster.id)
  ##Generalized score test
  dimY<-nrow(Y.res);dimA<-ncol(Z.A);dimI<-ncol(Z.I);score<-matrix(0,dimI,1)
  An<-matrix(0,dimA+dimI,dimA+dimI);Bn<-matrix(0,dimA+dimI,dimA+dimI)
  n.total<-1;n.rep<-as.numeric(table(time[,1]))
  for (i in 1:m)
  {
    ni<-n.rep[i]
    index<-n.total:(n.total+ni-1)
    n.total<-n.total + ni
    #working covariance matrix
    if (corstr=="exchangeable"){Vi<-diag(1-rho,ni)+ matrix(rho, ni, ni);Vi.inv<-solve(Vi)}
    if (corstr=="ar1"){Vi<-rho^abs(outer(time[index,2], time[index,2],"-"));Vi.inv<-solve(Vi)}
    if (corstr=="independence"){Vi<-diag(1,ni);Vi.inv<-Vi}
    #if (link=="logit"){Ci<-mu[index]*(1-mu[index]);Vi.inv<-outer(sqrt(Ci),sqrt(Ci)^-1,"*")*Vi.inv}
    #score matrices
    Z.Ai<-Z.A[index,];Z.Ii<-Z.I[index,]
    if (ni==1) { Z.Ai<-t(Z.Ai); Z.Ii<-t(Z.Ii)}
    Zi<-cbind(Z.Ii,Z.Ai)
    #An<-An+t(Zi)%*%Vi.inv%*%Zi
    An[(dimI+1):(dimI+dimA),(dimI+1):(dimI+dimA)]<-An[(dimI+1):(dimI+dimA),(dimI+1):(dimI+dimA)]+t(Z.Ai)%*%Vi.inv%*%Z.Ai
    An[1:dimI, (dimI+1):(dimI+dimA)]<-An[1:dimI, (dimI+1):(dimI+dimA)]+t(Z.Ii)%*%(Vi.inv%*%Z.Ai)
    Bntemp<-t(Zi)%*%(Vi.inv%*%Y.res[index])
    Bn<-Bn+Bntemp%*%t(Bntemp)
    score<-score+Bntemp[1:dimI,]#t(Z.Ii)%*%Vi.inv%*%Y.res[index]
  }
  An<-An/m; Bn<-Bn/m
  #Calculate p-value
  A00<-An[(dimI+1):(dimI+dimA),(dimI+1):(dimI+dimA)];A00.inv<-solve(A00);A10<-An[1:dimI, (dimI+1):(dimI+dimA)]
  B11<-Bn[1:dimI,1:dimI];B10<-Bn[1:dimI,(dimI+1):(dimI+dimA)];B00<-Bn[(dimI+1):(dimI+dimA),(dimI+1):(dimI+dimA)]
  P<-A10%*%A00.inv;M<-B11-t(B10%*%t(P))-B10%*%t(P)+P%*%B00%*%t(P)
  #Revise the covariance matrix using m-q
  M<-M*m/(m-dimA)
  #Revise the score using standard representation
  score<-score/sqrt(m)
  return(list(score=score,M=M))#Q: test statistic; M: covariance matrix
}

# calculate number of independent tests
effect.n<-function(x,MinP.adjust){
  temp<-0;sum.EV<-sum(x) #summation of eigen values
  for (i in 1:length(x)){
    temp<-temp+x[i];if (temp>sum.EV*MinP.adjust){break}
  }
  return(i)
}

#Data management
Minor.allele<-function(x){if(mean(x)>1){x<-2-x};return(x)}
Standardize<-function(x){return((x-mean(x))/sd(x))}
Center<-function(x){return(x-mean(x))}
Variation<-function(x){x<-apply(x,2,Minor.allele);x<-x[,which(apply(x,2,sd)>0)];return(x)}
Common<-function(x){x<-apply(x,2,Minor.allele);freq<-apply(x,2,mean)/2;x<-x[,which(freq>0.05)];return(x)}
##Matrix calculation
#Get sqrt.root of a matrix
Get.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}
#Get inverse of a matrix; numerically robust
Get.inverse<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix inverse!")}
  a.inverse <- a.eig$vectors[,ID1] %*% diag(a.eig$values[ID1]^-1) %*% t(a.eig$vectors[,ID1])
  return(a.inverse)
}
#Get -1/2 of a matrix; numerically robust
Get.inv.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){stop("Error to obtain matrix square!")}
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])^-2) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}
#Time similarity
Check.same<-function(x,y){return(as.numeric(x==y))}
Check.near<-function(x,y){return(as.numeric(abs(x-y)==1))}
D.matrix<-function(X){
  n<-length(X[,1]);temp<-matrix(0,n,n)
  #checking time
  temp<-outer(X[,1],X[,1],Check.same)*outer(X[,2],X[,2],Check.near)
  diag(temp)<-0
  return(temp)
}
# Simple Imputation (from SKAT package)
# Z : an m x p genotype matrix with m samples and p SNPs

Impute<-function(Z, impute.method){
  p<-dim(Z)[2]
  if(impute.method =="random"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-rbinom(length(IDX),2,maf1)
      }
    }
  } else if(impute.method =="fixed"){
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-2 * maf1
      }
    }
  } else if(impute.method =="bestguess") {
    for(i in 1:p){
      IDX<-which(is.na(Z[,i]))
      if(length(IDX) > 0){
        maf1<-mean(Z[-IDX,i])/2
        Z[IDX,i]<-round(2 * maf1)
      }
    }
  } else {
    stop("Error: Imputation method shoud be \"fixed\", \"random\" or \"bestguess\" ")
  }
  return(as.matrix(Z))
}



