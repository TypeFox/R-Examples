
#########################
#   Null model fitting
#########################

null.LGRF<-function(Y,time,X=NULL) #Y: n*1 matrix, time:n*2 matrix, X:Environment n*p matrix
{
  ##Preliminary
  Y<-as.matrix(Y);n<-nrow(Y);index.cluster<-unique(time[,1])
  if(length(X)==0){is.naX<-0;nrowX<-nrow(Y)}else{
    is.naX<-is.na(X);nrowX<-nrow(X)
  }
  N_MISS<-sum(is.naX+sum(is.na(Y))+sum(is.na(time)))
  if(N_MISS>0){
    msg<-sprintf("X, Y or time has missing values")
    stop(msg)
  }
  if((nrow(Y)-nrow(time))^2+(nrow(Y)-nrowX)^2>0){
    msg<-sprintf("Number of observations in Y, X and time does not match")
    stop(msg)
  }
  
  X1<-cbind(rep(1,n),X);p<-ncol(X1)-1
  #Iteration to get beta and eta
  step.error<-Inf;beta<-rep(0,p+1);eta<-0;iter<-0 #iteration time
  while (step.error>0.00001)
  {
    A<-0;B<-0;x.field<-NULL
    for (i in index.cluster)
    {
      index<-which(time[,1]==i);ni<-length(index)
      #time similarity
      S.subject<-matrix(1,ni,ni);diag(S.subject)<-0
      V<-diag(1,length(index))-eta*S.subject
      X.temp<-X1[index,];Y.temp<-Y[index,]
      if (ni==1) {X.temp<-t(X.temp)}
      A<-A+t(X.temp)%*%V%*%X.temp;B<-B+t(X.temp)%*%V%*%Y.temp
    }
    #calculating residual
    beta.new<-solve(A)%*%B;Y.res<-Y-X1%*%beta.new
    if(length(index.cluster)==nrow(time)) #check the existence of repeated measurements
    {
      eta.new<-0
    }
    if(length(index.cluster)!=nrow(time))
    {
      #Generate RF psuedo variables
      for (i in index.cluster)
      {
        index<-which(time[,1]==i);ni<-length(index)
        #time similarity
        S.subject<-matrix(1,ni,ni);diag(S.subject)<-0
        x.field<-rbind(x.field,S.subject%*%Y.res[index,])
      }
      eta.new<-lm(Y.res~0+x.field)$coefficients
    }
    #calculate step error
    step.error<-sum(c(beta.new-beta,eta.new-eta)^2);iter<-iter+1
    beta<-beta.new;eta<-eta.new
  }
  return(list(Y.res=Y.res,time=time,X1=X1,eta=eta))
}

########################
#   Region test
########################

# note that genotype Z is not in long form.
test.LGRF<-function(Z,result.null,Gsub.id=NULL,interGXT=FALSE,similarity='GR',impute.method="fixed") #Z: n*t matrix, result.null: result of null.longGRF
{

  #interGXT=FALSE; similarity='GR'
  # recall data
  Y.res<-result.null$Y.res;time<-result.null$time
  X1<-result.null$X1;eta<-result.null$eta
  # match the phenotype and genotype subject ID
  if(length(Gsub.id)==0){Gsub.id<-unique(time[,1])}
  Z<-as.matrix(Z[match(unique(time[,1]),Gsub.id),])
  
  # missing genotype imputation
  N_MISS<-sum(is.na(Z))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(Z)/ncol(Z))
    warning(msg,call.=F)
    Z<-Impute(Z,impute.method) 
  }
  #if((length(unique(time[,1]))-nrow(Z))^2>0){
  #  msg<-sprintf("Number of subjects in Y, X, time and Z does not match")
  #  stop(msg)
  #}
  # generate long form genotype
  ID_G<-unique(time[,1]);Z<-as.matrix(Z[match(time[,1],ID_G),])
  
  # tests
  if (interGXT==FALSE)
  {
    # generate IBS pseudo variable
    if (similarity=='IBS'){G<-IBS_pseudo(Z)}
    # centered genotype
    if (similarity=='GR'){G<-apply(Z,2,center)}
    # preliminary
    dimY<-nrow(Y.res);m<-length(unique(time[,1]));dimX<-ncol(X1);dimG<-ncol(G) #number of SNP/pseudo variable
    index.cluster<-unique(time[,1])
    # generalized score test
    dS<-matrix(0,2*dimG+dimX,dimX);
    S1<-matrix(0,2*dimG,1);S0<-matrix(0,dimX,1);D<-matrix(0,2*dimG+dimX,2*dimG+dimX)
    for (i in index.cluster)
    {
      index<-which(time[,1]==i);ni<-length(index)
      #time similarity
      S.subject<-matrix(1,ni,ni);diag(S.subject)<-0
      Vi<-diag(1,ni)-eta*S.subject
      Gi<-G[index,];Xi<-X1[index,];Yi<-as.matrix(Y.res[index,])
      if (ni==1) {Gi<-t(Gi);Xi<-t(Xi)}
      G.temp<-cbind(Vi%*%Gi,Gi,Xi)
      # calculate cluster based statistics
      Si<-t(G.temp)%*%Yi;dSi<--t(G.temp)%*%Xi;Di<-Si%*%t(Si)
      # calculate overall statistics
      S1<-S1+Si[1:(2*dimG),];dS<-dS+dSi;D<-D+Di
    }
    D<-D/m;S1<-S1/sqrt(m);dS<-dS/sqrt(m)
    # test statistic
    Q<-t(S1[(dimG+1):(2*dimG),])%*%S1[1:dimG,]
    # null distribution
    P<-dS[1:(2*dimG),]%*%solve(dS[(2*dimG+1):(2*dimG+dimX),])
    D11<-D[1:(2*dimG),1:(2*dimG)];D10<-D[1:(2*dimG),(2*dimG+1):(2*dimG+dimX)];D00<-D[(2*dimG+1):(2*dimG+dimX),(2*dimG+1):(2*dimG+dimX)]
    M<-D11-P%*%t(D10)-D10%*%t(P)+P%*%D00%*%t(P)
    M.rotate<-cbind(M[,(dimG+1):(2*dimG)],M[,1:dimG])
    eigen.M<-Re(eigen(M.rotate,symmetric=FALSE,only.values=TRUE)$values)
    # calculate p.value
    p.value<-davies(2*Q,eigen.M)$Qq
  }
  else
  {
    # generate IBS pseudo variable
    if (similarity=='IBS'){G<-IBS_pseudo(Z)}
    # centered genotype
    if (similarity=='GR'){G<-apply(Z,2,center)}
    # interaction terms
    Gt<-apply(Z*time[,2],2,center)
    # preliminary
    dimY<-nrow(Y.res);m<-length(unique(time[,1]));dimX<-ncol(X1);dimG<-(ncol(G)+ncol(Gt))/2 #number of SNP/pseudo variable. PS: for convinience, we use this for the revision without changing other part
    index.cluster<-unique(time[,1])
    # calculate weights
    K.g<-t(G)%*%G;K.gt<-t(Gt)%*%Gt
    s2.g<-2*sum(diag(K.g%*%K.g));s2.gt<-2*sum(diag(K.gt%*%K.gt))
    alpha.g<-sqrt(s2.gt/(s2.g+s2.gt));alpha.gt<-sqrt(s2.g/(s2.g+s2.gt))
    # generalized score test
    dS<-matrix(0,4*dimG+dimX,dimX);
    S1<-matrix(0,4*dimG,1);S0<-matrix(0,dimX,1);D<-matrix(0,4*dimG+dimX,4*dimG+dimX)
    for (i in index.cluster)
    {
      index<-which(time[,1]==i);ni<-length(index)
      #time similarity
      S.subject<-matrix(1,ni,ni);diag(S.subject)<-0
      Vi<-diag(1,ni)-eta*S.subject
      Gi<-G[index,];Gti<-Gt[index,];Xi<-X1[index,];Yi<-as.matrix(Y.res[index,])
      if (ni==1) {Gi<-t(Gi);Gti<-t(Gti);Xi<-t(Xi)}
      G.temp<-cbind(sqrt(alpha.g)*Vi%*%Gi,sqrt(alpha.gt)*Vi%*%Gti,sqrt(alpha.g)*Gi,sqrt(alpha.gt)*Gti,Xi)
      #calculate cluster based statistics
      Si<-t(G.temp)%*%Yi;dSi<--t(G.temp)%*%Xi;Di<-Si%*%t(Si)
      #calculate overall statistics
      S1<-S1+Si[1:(4*dimG),];dS<-dS+dSi;D<-D+Di
    }
    D<-D/m;S1<-S1/sqrt(m);dS<-dS/sqrt(m)
    #test statistic
    Q<-t(S1[(2*dimG+1):(4*dimG),])%*%S1[1:(2*dimG),]
    #null distribution
    P<-dS[1:(4*dimG),]%*%solve(dS[(4*dimG+1):(4*dimG+dimX),])
    D11<-D[1:(4*dimG),1:(4*dimG)];D10<-D[1:(4*dimG),(4*dimG+1):(4*dimG+dimX)];D00<-D[(4*dimG+1):(4*dimG+dimX),(4*dimG+1):(4*dimG+dimX)]
    M<-D11-P%*%t(D10)-D10%*%t(P)+P%*%D00%*%t(P)
    M.rotate<-cbind(M[,(2*dimG+1):(4*dimG)],M[,1:(2*dimG)])
    eigen.M<-Re(eigen(M.rotate,symmetric=FALSE,only.values=TRUE)$values)
    #calculate p.value
    p.value<-davies(2*Q,eigen.M)$Qq
  }
  return(list(p.value=p.value,n.marker=dimG))
}

# generate IBS pseudo variable
IBS_pseudo<-function(x)
{
  pseudo1<-(0<=x & x<0.5)*sqrt(2)+(0.5<=x & x<1.5)*sqrt(2)/2
  pseudo2<-(0.5<=x & x<1.5)*sqrt(2)/2+(1.5<=x & x<2.5)*sqrt(2)
  pseudo3<-(0.5<=x & x<1.5)
  pseudo<-cbind(pseudo1,pseudo2,pseudo3)
  return(pseudo)
}

# minP test for comparision
test.MinP<-function(Z,result.null,Gsub.id=NULL,corstr="exchangeable",MinP.adjust=0.95,impute.method="fixed"){
  Y.res<-result.null$Y.res;time<-result.null$time
  X1<-result.null$X1;eta<-result.null$eta;dimG<-ncol(Z)
  # match the phenotype and genotype subject ID
  if(length(Gsub.id)==0){Gsub.id<-unique(time[,1])}
  Z<-as.matrix(Z[match(unique(time[,1]),Gsub.id),])
  # missing genotype imputation
  N_MISS<-sum(is.na(Z))
  if(N_MISS>0){
    msg<-sprintf("The missing genotype rate is %f. Imputation is applied.", N_MISS/nrow(Z)/ncol(Z))
    warning(msg,call.=F)
    Z<-Impute(Z,impute.method) 
  }
  # check if the matrices match
  if((length(unique(time[,1]))-nrow(Z))^2>0){
    msg<-sprintf("Number of subjects in Y, X, time and Z does not match")
    stop(msg)
  }
  # check the allele frequency
  MAF<-apply(Z,2,mean)/2
  if(min(MAF)<0.01){
    msg<-sprintf("There are rare variants with MAF<0.01. MinP by GEE might be unstable", N_MISS/nrow(Z)/ncol(Z))
    warning(msg,call.=F)
  }
  # generate long form genotype
  ID_G<-unique(time[,1]);Z<-Z[match(time[,1],ID_G),]
  # MinP tests by GEE
  temp.pvalue<-rep(0,ncol(Z))
  for (k in 1:ncol(Z)){
    if(ncol(X1)!=1){temp.fit<-geeglm(Y.res~X1[,-1]+Z[,k],id=time[,1],corstr=corstr)}else{
      temp.fit<-geeglm(Y.res~Z[,k],id=time[,1],corstr=corstr)
    }
    temp.pvalue[k]<-summary(temp.fit)$coefficient[nrow(summary(temp.fit)$coefficient),4]
  }
  # calculate number of independent tests
  eigen.G<-eigen(cor(Z),only.values=TRUE)$values
  me.G<-effect.n(eigen.G,MinP.adjust)
  # calculate adjusted p.value
  p.value<-min(1,min(temp.pvalue)*me.G)
  return(list(p.value=p.value,n.marker=dimG))
}

# calculate number of independent tests
effect.n<-function(x,MinP.adjust){
  temp<-0;sum.EV<-sum(x) #summation of eigen values
  for (i in 1:length(x)){
    temp<-temp+x[i];if (temp>sum.EV*MinP.adjust){break}
  }
  return(i)
}


# Some functions used below is from the SKAT package

Get.sqrt<-function(A){
  a.eig <- eigen(A,symmetric=TRUE)
  ID1<-which(a.eig$values > 0)
  if(length(ID1)== 0){
    stop("Error to obtain matrix square!")
  }
  a.sqrt <- a.eig$vectors[,ID1] %*% diag(sqrt(a.eig$values[ID1])) %*% t(a.eig$vectors[,ID1])
  return(a.sqrt)
}

#Time similarity
Check.same<-function(x,y)
{return(as.numeric(x==y))}

Check.near<-function(x,y)
{return(as.numeric(abs(x-y)==1))}

#Centering
center<-function(x)
{return(x-mean(x))}

# Simple Imputation (from SKAT package)
# Z : an n x p genotype matrix with n samples and p SNPs
# Missing : a missing genotype value. Default is 9

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



