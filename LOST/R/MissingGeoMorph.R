MissingGeoMorph <-
function(x, nlandmarks, method="BPCA"){
  
  ##########################################
  ########### internal functions############
  ##########################################
  
  format.array<-function(dataset,nlandmarks){
    dataset<-as.matrix(dataset)
    nspecimen<-nrow(dataset)/nlandmarks
    start<-seq(from=1,to=nrow(dataset),by=nlandmarks)
    sparray<-array(dataset,dim<-c(nlandmarks,2,nspecimen))
    for (k in 1:nspecimen){
      x<-start[k]
      y<-x+nlandmarks-1
      single<-dataset[x:y,]
      sparray[,,k]<-single
    }
    return(sparray)
  }
  
  
  format.matrix<-function(x,nlandmarks){
    N=nlandmarks*2
    nspecimen<-nrow(x)/nlandmarks
    datavector<-unmatrix(x,byrow=TRUE)
    new.matrix<-matrix(datavector,nspecimen,N,byrow=TRUE)
    return(new.matrix)
  }
  
  
  
  
  
  unformat<-function(x,nlandmarks){
    nsp<-length(x)/(nlandmarks*2)
    rows<-nsp*nlandmarks
    if(nsp==1){matrix<-matrix(x,rows,2,byrow=TRUE); return(matrix)}
    else{
      vector<-unmatrix(x,byrow=TRUE)
      matrix<-matrix(vector,rows,2,byrow=TRUE)
      return(matrix)
    }
  }
  
  
  
  
  
  
  
  ##### tps2d is from Claude, J. (2008). Thin-Plate Splines. Morphometrics with R pp. 181-188. Springer, New York.  
  tps2d<-function(M,matr,matt){
    p<-dim(matr)[1];q<-dim(M)[1];n1<-p+3
    P<-matrix(NA,p,p)
    for (i in 1:p){
      for (j in 1:p){
        r2<-sum((matr[i,]-matr[j,])^2)
        P[i,j]<-r2*log(r2)
      }
    }
    P[which(is.na(P))]<-0
    Q<-cbind(1,matr)
    L<-rbind(cbind(P,Q),cbind(t(Q),matrix(0,3,3)))
    m2<-rbind(matt,matrix(0,3,2))
    coefx<-ginv(L)%*%m2[,1]
    coefy<-ginv(L)%*%m2[,2]
    fx<-function(matr,M,coef1){
      Xn<-numeric(q)
      for (i in 1:q){
        Z<-apply((matr-matrix(M[i,],p,2,byrow=T))^2,1,sum)
        U<-sum(coef1[1:p]*(Z*log(Z+1)))
        U[which(is.na(U))]<-0
        Xn[i]<-coef1[p+1]+coef1[p+2]*M[i,1]+coef1[p+3]*M[i,2]+U
      }
      Xn
    }
    matg<-matrix(NA,q,2)
    matg[,1]<-fx(matr,M,coefx)
    matg[,2]<-fx(matr,M,coefy)
    return(matg)
    
  }
  
  
  missing.tps<-function(X,nlandmarks){
    
    
    
    
    char.com<-matrix(ncol=ncol(X),nrow=nrow(X))
    xvalues<-X[,1]
    nspecimen<-nrow(X)/nlandmarks
    xmat<-matrix(xvalues,ncol=nlandmarks,nrow=nspecimen,byrow=TRUE)
    missguide<-ifelse((is.na(xmat))==TRUE,1,0)
    whichpoints<-missguide*matrix(1:nlandmarks,ncol=nlandmarks,nrow=nspecimen,byrow=TRUE)
    eachmissing<-apply(missguide,1,sum)
    incompletes<-setdiff(ifelse(eachmissing==0,0,1)*1:nspecimen,0)
    
    completes<-complete.specimens(X,nlandmarks)
    GPA.com<-procGPA(format.array(completes,nlandmarks),pcaoutput=FALSE,distances=FALSE)
    aligned.com<-GPA.com$rotated
    mean.aligned<-GPA.com$mshape
    
    aligned.mat<-aligned.com[,,1]
    for(z in 2:dim(aligned.com)[3]){
      aligned.mat<-rbind(aligned.mat,aligned.com[,,z])
    }
    
    for (i in 1:length(incompletes)){
      current<-incompletes[i]  
      cur.spec<-X[((current-1)*nlandmarks+1):(current*nlandmarks),]
      cur.miss<-setdiff(whichpoints[current,],0)
      cur.spec.miss<-cur.spec[-cur.miss,]
      mean.miss<-mean.aligned[-cur.miss,]
      tps.aligned<-tps2d(mean.aligned,mean.miss,cur.spec.miss)
      cur.spec2<-cur.spec
      cur.spec2[is.na(cur.spec2)==TRUE]<-0
      missing<-ifelse(is.na(cur.spec)==TRUE,1,0)
      opposite<-ifelse(is.na(cur.spec)==TRUE,0,1)
      new.spec<-cur.spec2*opposite+tps.aligned*missing    
      starts<-(current-1)*nlandmarks+1
      stops<-current*nlandmarks
      char.com[starts:stops,]<-as.matrix(new.spec)
    }
    
    completed<-setdiff(1:nspecimen,incompletes)
    for (f in 1:length(completed)){
      com.n<-completed[f]
      com.spec<-aligned.com[,,f]
      starts<-(com.n-1)*nlandmarks+1
      stops<-com.n*nlandmarks
      char.com[starts:stops,]<-com.spec
    }
    
    return(char.com)
  }
  
  
  est.reg1<-function(crocmiss,col_indep){
    ##creating a fill matrix
    cols<-ncol(crocmiss);rows<-nrow(crocmiss);estimated_matrix<-matrix(ncol=cols,nrow=rows)
    if (col_indep==1){deps<-2:cols} else 
      if (col_indep==cols){deps<-1:(cols-1)} else {
        deps<-c(1:(col_indep-1),(col_indep+1):cols)} 
    ndeps<-length(deps)
    indep<-crocmiss[,col_indep]
    rcoefs<-numeric()
    for (i in 1:ndeps){
      vari<-deps[i]
      lm_fit<-lm(crocmiss[,vari]~indep); lm_sum<-summary.lm(lm_fit)
      rcoefs[i]<-sqrt(lm_sum$r.squared)
    }
    ranks<-rank(rcoefs)
    newindep<-indep
    ## run lms and fill in missing values for independent variable
    for (m in 1:length(deps)){
      a<-m-1
      whichvari<-ifelse(ranks==(length(rcoefs)-a),1,0); strongest<-sum(deps*whichvari)
      indeps_lm<-lm(newindep~crocmiss[,strongest])
      indeps_coef<-indeps_lm$coefficients
      logestimate_indep<-indeps_coef[1]+indeps_coef[2]*(crocmiss[,strongest])
      estimate_indep<-logestimate_indep
      missings<-ifelse(is.na(newindep),1,0)
      nonmissings<-ifelse(is.na(newindep),0,1)
      fillnonmissing<-ifelse(is.na(newindep),0,newindep)
      fillmissing<-ifelse(is.na(newindep),estimate_indep,0)
      newindep<-fillmissing+fillnonmissing
    }
    return(newindep)
  }
  
  best.reg<-function(x){
    estimated.matrix<-matrix(ncol=ncol(x),nrow=nrow(x))
    for (i in 1:ncol(x)){
      estimated.matrix[,i]<-est.reg1(x,i)
    }
    return(estimated.matrix)
  }
  
  complete.specimens<-function(dataset,nlandmarks){
    base<-c(1,1)
    included<-base
    excluded<-base
    nspecimen<-nrow(dataset)/nlandmarks
    start<-seq(from=1,to=nrow(dataset),by=nlandmarks)
    nsp<-1:nspecimen
    for (k in 1:nspecimen){
      x<-start[k]
      y<-x+nlandmarks-1
      single<-dataset[x:y,]
      reduced<-na.omit(single)
      rows<-nrow(reduced)
      if (rows==nlandmarks){included<-rbind(included,single)}
      else {excluded<-rbind(excluded,single)}
    }
    end<-nrow(included)
    included<-included[2:end,]
    return (included)
  }
  
########### MissingGeoMorph ##############
####### applying internal functions #####
  
  
  aligned<-align.missing(x,nlandmarks)
  new.matrix<-format.matrix(aligned,nlandmarks)
  
  
  if(method=="BPCA"){
    aligned<-align.missing(x,nlandmarks)
    new.matrix<-format.matrix(aligned,nlandmarks)
    estimator<-pca(new.matrix,method="bpca")
    estimated.values<-completeObs(estimator)
    results<-unformat(estimated.values,nlandmarks)
  } else if (method=="mean"){
    aligned<-align.missing(x,nlandmarks)
    new.matrix<-format.matrix(aligned,nlandmarks)
    estimated.values<-impute(new.matrix,what="mean")
    results<-unformat(estimated.values,nlandmarks)
  } else if (method=="reg"){
    aligned<-align.missing(x,nlandmarks)
    new.matrix<-format.matrix(aligned,nlandmarks)
    estimated.values<-best.reg(new.matrix)
    results<-unformat(estimated.values,nlandmarks)
  } else if (method=="TPS"){
    results<-missing.tps(aligned,nlandmarks)
  }
  return(results)
  
  
  
  
}
