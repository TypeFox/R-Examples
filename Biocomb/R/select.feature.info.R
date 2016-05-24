compute.aucs <- function(dattable){

  labs <- dattable[,ncol(dattable)]
  aucvals <- rep(0,ncol(dattable)-1)

  val=levels(labs)
  #pos <- rep(val[2],ncol(dattable)-1)
  pos<-factor(rep(val[2],ncol(dattable)-1),levels=val)
  if(length(val)==2)
  {
    for (i in 1:(ncol(dattable)-1)){
      pred <- prediction(dattable[,i],labs)
      aucv <- performance(pred,"tpr", "fpr",measure="auc")
      aucval <- attr(aucv,"y.values")[[1]]
      if (aucval<0.5){
        aucval <- 1-aucval
        pos[i] <- val[1]  ####Positive Correlation ist richtig, AUC-Wert >=0.5, sonst AUC-Wert <0.5
      }
      aucvals[i] <- aucval
    }
    auctab<-data.frame(names(dattable)[1:(ncol(dattable)-1)],aucvals,pos)
    names(auctab)<-c("Biomarker","AUC","Positive class")
  }
  else
  {
    for (i in 1:(ncol(dattable)-1)){
      aucval <- multiclass.roc(labs,dattable[,i])$auc
      aucval2 <- multiclass.roc(labs,dattable[,i])$auc
      if (aucval<aucval2){
        aucval <- aucval2
      }
      aucvals[i] <- aucval
    }
    auctab<-data.frame(names(dattable)[1:(ncol(dattable)-1)],aucvals)
    names(auctab)<-c("Biomarker","AUC")
  }

  return(auctab)
}

chi2.algorithm<- function(matrix,attrs.nominal,threshold)
{
  dd=dim(matrix)

  if(length(attrs.nominal)>0)
  {
  for(i in 1:length(attrs.nominal))
  {
    matrix[,attrs.nominal[i]]=as.factor(matrix[,attrs.nominal[i]])
  }
  }
  #for inconsistency for nominal
  vrem.nominal=matrix[,attrs.nominal,drop=FALSE]
  if(length(attrs.nominal)>0)
  {
  for(i in 1:length(attrs.nominal))
  {
    vrem.nominal[,i]=as.numeric(vrem.nominal[,i])
  }
  }
  #-------
  data=matrix[,-c(attrs.nominal,dd[2]),drop=FALSE]
  data.start=data
  class=matrix[,dd[2]]
  class=as.character(class)

  d1=dim(data)

  label=unique(class)
  mat.int=matrix(0,2,length(label))
  colnames(mat.int)=label

  int.list=fun1_chi(data,class)
  int.list.start=int.list

  #Phase 1
  sig.value=0.6
  df=length(label)-1
  chi.value=qchisq(1-sig.value, df=df)

  chi.stat=fun2_chi(int.list,mat.int)
  chi.stat.start=chi.stat
  len_chi=sapply(chi.stat, function(z) length(z))

  incons=0

  step=0.1
  delta=6
  shag=1
  calc=0
  while(incons<=threshold)
  {
    sig.value0=sig.value

    if(shag==delta)
    {
      step=step*0.1
      delta=delta+9
    }

    sig.value=sig.value-step
    shag=shag+1

    if(sig.value<0.000000000011)
    {
      #browser()
    }
    chi.value=qchisq(1-sig.value, df=df)

    check=sapply(chi.stat,function(z) length(z))
    if(all(check==0))
    {
      break
    }

    out3=fun3_chi(chi.stat,int.list,data, chi.value, mat.int)
    data=out3$data
    chi.stat=out3$chi_stat
    int.list=out3$int_list

    incons=check_incons(data, vrem.nominal,class)
    calc=calc+1
  }

  #Phase 2
  data=data.start

  sig.attr=rep(sig.value0,d1[2])
  chi.value=qchisq(1-sig.value0, df=df)
  chi.attr=rep(chi.value,d1[2])

  int.list=int.list.start
  chi.stat=chi.stat.start

  data=fun4_chi(chi.stat,int.list,data,vrem.nominal,chi.attr,sig.attr,class,mat.int,threshold,df,step,delta,shag)

  rr=sapply(1:d1[2],function(z) length(unique(data[,z])))
  data.out=data[,which(rr>1),drop=FALSE]
  data.out=cbind(data.out,matrix[,attrs.nominal,drop=FALSE])

  return(list(data.out=data.out,subset=colnames(data.out)))
}


select.forward.Corr<- function(matrix,disc.method,attrs.nominal)
{
  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3

  dd=dim(m3)
  if(dd[2]>1)
  {
  subset <- forward_path(0:(ncol(m3)-2),m3)
  subset=subset+1
  subset<-names(m3)[subset]
  }
  else
  {
    subset <- NULL
  }
  return(subset)
}

select.forward.wrapper<- function(dattable)
{
  evaluator <- function(subset) {
    #k-fold cross validation
    results = sapply(1:k, function(i) {
      test.idx <- testind[,i]
      train.idx <- !test.idx
      test <- dattable[test.idx, , drop=FALSE]
      train <- dattable[train.idx, , drop=FALSE]
      tree <- rpart(as.simple.formula(subset, names(dattable)[ncol(dattable)]), train,method = "class")
      error.rate = sum(test[,ncol(dattable)] != predict(tree, test, type="c")) / nrow(test)
      return(1 - error.rate)
    })
    return(mean(results))
  }

  k <- 5
  splits <- runif(nrow(dattable))
  testind <- sapply(1:k, function(i) {(splits >= (i - 1) / k) & (splits < i / k)})
  subset <- forward.search(names(dattable)[-ncol(dattable)], evaluator)
}


CalcGain<-function(m1,m2,symm)
{
  dd=length(m1)
  fq1=table(m1)
  fq1=fq1/dd[1]

  entropyF1=-sapply(fq1, function(z) if(z==0) 0 else z*log(z))
  entropyF1=sum(entropyF1)

  fq2=table(m2)
  fq2=fq2/dd[1]

  entropyF2=-sapply(fq2, function(z) if(z==0) 0 else z*log(z))
  entropyF2=sum(entropyF2)

  fq=table(m1,m2)

  entropyF12=0

  for(i in 1:length(fq2))
  {
    fq0=fq[,i]/sum(fq[,i])
    vrem=-sapply(fq0,function(z) if(z==0) 0 else z*log(z))
    entropyF12=entropyF12+(fq2[i])*sum(vrem)
  }

  entropy=entropyF1-entropyF12
  if(symm)
  {
    if((entropyF1+entropyF2)==0)
    {
      entropy=0
    }
    else
    {
      entropy=2*entropy/(entropyF1+entropyF2)
    }
  }
  return(entropy)
}


ProcessData1<-function(matrix,disc.method,attrs.nominal)
{
  dd=dim(matrix)
  matrix=data.frame(matrix)

  matrix[,dd[2]]=as.factor(matrix[,dd[2]])
  #data=matrix

  if(disc.method=="MDL")
  {
    m3 <- Discretize(as.formula(paste(names(matrix)[dd[2]],"~.")), data = matrix)
    #m3<-mdlp(matrix)$Disc.data
  }

  if(disc.method=="equal frequency")
  {
    m3=matrix
    for(i in 1:(dd[2]-1))
    {
      if(!(i%in%attrs.nominal))
      {
        m3[,i] <- discretize(matrix[,i], method="frequency",categories=3)
      }
    }
  }

  if(disc.method=="equal interval width")
  {
    m3=matrix
    for(i in 1:(dd[2]-1))
    {
      if(!(i%in%attrs.nominal))
      {
        m3[,i] <- discretize(matrix[,i], categories=3)
      }
    }
  }
  #-------------------

  #extract the features with one interval

  sel.one=lapply(m3, function(z) (length(levels(z))==1)&&(levels(z)=="'All'"))

  sel.one=which(unlist(sel.one)==TRUE)

  #selected features
  sel.feature=1:dd[2]
  if(length(sel.one)>0)
  {
    sel.feature=sel.feature[-sel.one]

    matrix=matrix[,-sel.one,drop=FALSE]
    m3=m3[,-sel.one,drop=FALSE]
  }
  return (list(m3=m3,sel.feature=sel.feature))
}

ProcessData<-function(matrix,disc.method,attrs.nominal,flag=FALSE)
{
  dd=dim(matrix)
  matrix=data.frame(matrix)

  matrix[,dd[2]]=as.factor(matrix[,dd[2]])
  #data=matrix

  if(disc.method=="MDL")
  {
    m3 <- Discretize(as.formula(paste(names(matrix)[dd[2]],"~.")), data = matrix)
    #m3<-mdlp(matrix)$Disc.data
  }

  if(disc.method=="equal frequency")
  {
    m3=matrix
    for(i in 1:(dd[2]-1))
    {
      if(!(i%in%attrs.nominal))
      {
        m3[,i] <- discretize(matrix[,i], method="frequency",categories=3)
      }
    }
  }

  if(disc.method=="equal interval width")
  {
    m3=matrix
    for(i in 1:(dd[2]-1))
    {
      if(!(i%in%attrs.nominal))
      {
        m3[,i] <- discretize(matrix[,i], categories=3)
      }
    }
  }
  #-------------------
  sel.feature=1:dd[2]
  if(flag)
  {
  #extract the features with one interval

   sel.one=lapply(m3, function(z) (length(levels(z))==1)&&(levels(z)=="'All'"))

   sel.one=which(unlist(sel.one)==TRUE)

  #selected features

  if(length(sel.one)>0)
  {
    sel.feature=sel.feature[-sel.one]

    matrix=matrix[,-sel.one,drop=FALSE]
    m3=m3[,-sel.one,drop=FALSE]
  }
  }
  return (list(m3=m3,sel.feature=sel.feature))
}

select.cfs<-function(matrix)
{
val <- cfs(as.formula(paste(names(matrix)[ncol(matrix)]," ~ .")), matrix)
val<-sapply(val, function(z) which(names(matrix)==z))
info.val <- data.frame(names(matrix)[val],val)
names(info.val) <- c("Biomarker","Index")
return(info.val)
}


select.relief<-function(matrix)
{
val <- relief(as.formula(paste(names(matrix)[ncol(matrix)]," ~ .")), matrix,neighbours.count = 5, sample.size = 10)
val <- sort(val[[1]],decreasing=T,index.return=TRUE)
info.val <- data.frame(names(matrix)[val$ix[1:(ncol(matrix)-1)]],val$x,val$ix)
names(info.val) <- c("Biomarker","Weights","NumberFeature")
return(info.val)
}

select.inf.chi2<-function(matrix,disc.method,attrs.nominal)
{
  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  #algorithm
  dd=dim(m3)
  if(dd[2]>1)
  {
  #stat=sapply(1:(dd[2]-1), function(z) chisq.test(table(m3[,z],m3[,dd[2]]))$statistic)
  #names(stat)=colnames(m3[,1:(dd[2]-1)])
  #to compare
  weights <- chi.squared(as.formula(paste(names(m3)[dd[2]],"~.")), m3)

  #what features are selected
  res=sort(weights$attr_importance,decreasing = TRUE,index.return=TRUE)
  val=res$ix
  weights.sort=res$x
  num.feature=sel.feature[val] #val - sorting
  info=data.frame(names(m3)[val],weights.sort,num.feature)
  }
  else
  {
    info=data.frame(character(),numeric(),numeric())
  }
  names(info) <- c("Biomarker","ChiSquare","NumberFeature")
  return(info)
}

select.inf.symm<-function(matrix,disc.method,attrs.nominal)
{
  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  #algorithm
  dd=dim(m3)

  if(dd[2]>1)
  {
  #SU1=information.gain(names(matrix)[dd[2]]~., matrix) #package "FSelect"

  #entropy of feature
  entropy=c()
  class=m3[,dd[2]]

  for(j in 1:(dd[2]-1))
  {
    feature=m3[,j]

    #Function
    out=CalcGain(feature,class,TRUE)
    entropy=c(entropy,out)
    #--------
  }
  #what features are selected
  res=sort(entropy,decreasing = TRUE,index.return=TRUE)
  val=res$ix
  entropy.sort=res$x
  num.feature=sel.feature[val] #val - sorting
  info=data.frame(names(m3)[val],entropy.sort,num.feature)
  }
  else
  {
    info=data.frame(character(),numeric(),numeric())
  }
  names(info) <- c("Biomarker","SymmetricalUncertainty","NumberFeature")
  return(info)
}

select.inf.gain<-function(matrix,disc.method,attrs.nominal)
{
  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  #algorithm
  dd=dim(m3)
  if(dd[2]>1)
  {
  #SU1=information.gain(names(matrix)[dd[2]]~., matrix) #package "FSelect"

  #entropy of feature
  entropy=c()
  class=m3[,dd[2]]

  for(j in 1:(dd[2]-1))
  {
    feature=m3[,j]

    #Function
    out=CalcGain(feature,class,FALSE)
    entropy=c(entropy,out)
    #--------
  }
  #what features are selected
  res=sort(entropy,decreasing = TRUE,index.return=TRUE)
  val=res$ix
  entropy.sort=res$x
  num.feature=sel.feature[val] #val - sorting
  info=data.frame(names(m3)[val],entropy.sort,num.feature)
  }
  else
  {
    info=data.frame(character(),numeric(),numeric())
  }
  names(info) <- c("Biomarker","Information.Gain","NumberFeature")
  return(info)
}

select.fast.filter<-function(matrix,disc.method,threshold,attrs.nominal)
{

  #second package "RWeka"

  out=ProcessData(matrix,disc.method,attrs.nominal,FALSE)
  m3=out$m3
  sel.feature=out$sel.feature
  #algorithm
  dd=dim(m3)
  if(dd[2]>1)
  {
  #SU1=information.gain(names(matrix)[dd[2]]~., matrix) #package "FSelect"

  #entropy of feature
  entropy=c()
  class=m3[,dd[2]]

  for(j in 1:(dd[2]-1))
  {
  feature=m3[,j]

  #Function
  out=CalcGain(feature,class,FALSE)
  entropy=c(entropy,out)
  #--------
  }

  ind=sapply(entropy,function(z) z>=threshold)

  entropy=entropy[ind]
  m3=m3[,ind,drop=FALSE]

  index.F1=1

  res=sort(entropy,decreasing = TRUE,index.return=TRUE)
  val=res$ix
  entropy.sort=res$x

  while(index.F1<=length(val))
  {
    Fp=m3[,val[index.F1]]
    index.F2=index.F1+1
    while(index.F2<=length(val))
    {
    Fq=m3[,val[index.F2]]
    SUpq=CalcGain(Fp,Fq,FALSE)
    if(SUpq>=entropy.sort[index.F2])
    {
      val=val[-index.F2]
      entropy.sort=entropy.sort[-index.F2]
      index.F2=index.F2-1
    }
    index.F2=index.F2+1
    }
    index.F1=index.F1+1
  }

  #what features are selected, ind-features with SU(p,c)>threshold
  num.feature=sel.feature[ind]
  num.feature=num.feature[val] #val - sorting
  info=data.frame(names(m3)[val],entropy.sort,num.feature)
  }
  else
  {
    info=data.frame(character(),numeric(),numeric())
  }
  names(info) <- c("Biomarker","Information.Gain","NumberFeature")
  return(info)
}
