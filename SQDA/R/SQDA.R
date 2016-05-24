#' @importFrom limma lmFit eBayes makeContrasts contrasts.fit topTable
#' @importFrom PDSCE pdsoft
#' @importFrom mvtnorm dmvnorm
NULL
#' simulated example data
#' @title exampeldata
#' @name exampledata
#' @docType data
#' @keywords data
NULL

#' simulated training data
#' @title training data
#' @name train.data
#' @docType data
#' @keywords data
NULL

#' simulated testing data
#' @title testing data
#' @name test.data
#' @docType data
#' @keywords data
NULL


#' generate cross-validation ids
#' @param data data matrix with column names being the class labels and row names being the genes.
#' @param cv the cross-validation folds
#' @return cross-validation ids that can be used to split data into training data and testing data.
cross<-function(data=NULL,cv=5){
  if(is.null(data)){stop("input your data matrix!")}
  group<-unique(colnames(data))
  gid=1:ncol(data)
  for(g in group){
    id<-gid[colnames(data)==g]
    times=as.integer(length(id)/cv)+1
    gid[id]<-rep(1:cv,times=times)[1:length(id)]
  }
  return(gid)
}
#' Gene sorter
#' @param data data matrix with column names being the class labels and row names being the genes.
#' @return topTable data structure from limma.
sortgene<-function(data=NULL){
  if(is.null(data)){stop("input your data matrix!")}
  f<-as.factor(colnames(data))
  design<-model.matrix(~0+f)
  colnames(design)<-sub("^f","",x=colnames(design))
  name=unique(colnames(data))
  x=unlist(sapply(1:(length(name)-1),function(x){paste(name[x],"-",name[(x+1):length(name)],sep="")} ))
  contrast.matrix <- makeContrasts(contrasts=x,levels=design)
  fit <- lmFit(data, design)
  fit2 <- contrasts.fit(fit, contrast.matrix)
  fit2 <- eBayes(fit2)
  if(length(name)>2){
    table<-as.data.frame(cbind(F=fit2$F,p.value=fit2$F.p.value))
    table<-table[sort(table$F,decreasing=T,index.return=T)$ix,]
  }else{table<-topTable(fit2, coef=1, adjust.method="fdr",number=nrow(data))}
  return(table)
}
#' Cross-validation function
#' @param data data matrix with column names being the class labels and row names being the genes.
#' @param cv cross-validation folds.
#' @param lam a sequence of \code{lambda}'s. 
#' @return returns a list object with following item. 
#' \item{cv.error}{cross-validation errors for each \code{lambda}}
sGdaCV2<-function(data=NULL,cv=5,lam=0){
  gid<-cross(data,cv)
  pred<-NULL
  rate<-NULL
  name<-NULL
  temp=NULL
  for(j in 1:cv){
    test.data<-data[,gid==j]
    train.data<-data[,gid!=j]
    name<-c(name,colnames(test.data))
    temp<-c(temp,sGda(data=train.data,data.new=test.data,lam=lam)$pred)
  }
  pred<-temp
  rate<-sum(pred!=name)/length(name)
  result<-list(cv.error=rate)
  return(result)
}


#' Prediction function
#' @param data data matrix with column names being the class labels and row names being the genes.
#' @param data.new the new data needs to be predicted.
#' @param lam optimal \code{lambda} from cross-validation.
#' @return returns a list object with following items. 
#' \item{pred}{predictions for class labels on the new dataset}
#' \item{lik}{likelihood of each class on the new dataset}
sGda<-function(data=NULL,data.new=NULL,lam=0){
  
  group<-unique(colnames(data))
  prior<-rep(1,length(group))/length(group)
  name=colnames(data)
  nr<-nrow(data)
  data.list<-lapply(group,function(x){t(data[,grep(x,name)])})
  var=lapply(data.list,function(x){pdsoft(s=cov(x)+0.0001,lam=lam)$sigma})
  mu<-sapply(group,function(x){apply(data[,grep(x,name)],1,mean)})
  res=apply(data.new,2,function(x){
    sapply(1:length(group),function(t){dmvnorm(x, mu[,t],var[[t]], log=T)})
  })				
  res[which(is.nan(res),arr.ind=T)]=log(0)				 			
  rownames(res)<-group
  colnames(res)<-colnames(data.new)
  pred<-group[max.col(t(res))]
  return(list(pred=pred,lik=res))
}
#' Blockwise classifiers 
#' @param data data matrix with column names being the class labels and row names being the genes.
#' @param data.new the new data needs to be predicted.
#' @param len block size
#' @param times number of blocks
#' @param lam a sequence of \code{lambda}'s from cross-validation.
#' @return returns a list object with following items. 
#' \item{cv.error}{cross-validation errors for each block}
#' \item{pred}{predictions for class labels on the new dataset}
#' \item{lik}{likelihood of each class on the new dataset}
simpleAGG3<-function(data=NULL,data.new=NULL,len=100,times=100,lam=seq(0,0.1,length=10)){
  name=colnames(data.new)
  rate<-NULL
  pred<-NULL
  lik=NULL
  group<-unique(colnames(data))
  p=nrow(data)
  for(i in 1:times){
    seq<- (1:len)+len*(i-1)
    dat.sel=data[seq,]
    dat.sel.ts=data.new[seq,]
    rate.1<-NULL
    for(j in lam){
      rate.1<-c(rate.1,sGdaCV2(data=dat.sel,lam=j)$cv.error)
    }
    lamj<-lam[which.min(rate.1)]
    rate.1<-min(rate.1)
    res<-sGda(data=dat.sel,data.new=dat.sel.ts,lam=lamj)
    rate<-c(rate,rate.1)
    pred<-rbind(pred,res$pred)
    lik<-rbind(lik,res$lik)
  }
  
  return(list(cv.error=rate,pred=pred,lik=lik))
}

#' Spase Quadratic Discriminant Analysis
#' @param train.data data matrix with column names being the class labels and row names being the genes.
#' @param test.data the new data needs to be predicted.
#' @param len block size
#' @param lams a sequence of \code{lambda}'s from cross-validation.
#' @param presel pre-selection indicator.
#' @param prelam pre-selection sparisty parameter, only used when presel=T.
#' @param margin error margin for pre-selection, only used when presel=T.
#' @return returns a list object with following items. 
#' \item{pred}{predictions for class labels on the test.data}
#' \item{p}{the number of blocks selected}
#' @references The application of sparse estimation of covariance matrix to quadratic discriminant analysis. Jiehuan Sun and Hongyu Zhao.
#' @export
#' @examples 
#'data(exampledata)
#' res<-sQDA(train.data[1:100,],test.data[1:100,],lams=0.2,presel=FALSE)
#' sum(res$pred!=colnames(test.data))/ncol(test.data)  ##prediction error
#' res$p ## number of blocks selected
#' res$pred ## predicted class labels on test.data

sQDA<-function(train.data=NULL,test.data=NULL,len=100,lams=seq(0.02,1,length=10),presel=T,prelam=0.2,margin=0.05){
  group<-unique(colnames(train.data))
  table<-sortgene(train.data)
  seq<-as.numeric(rownames(table))
  dat.sel=train.data[seq,]
  dat.sel.ts=test.data[seq,]
  times=as.integer(nrow(train.data)/len)
  k<-(1:times)
  if(presel){
  res<-simpleAGG3(data=dat.sel,data.new=dat.sel.ts,len=len,times=times,lam=prelam)
  k<-(1:times)[res$cv.error <= min(res$cv.error)+margin]
  }
  lik<-NULL
  cverror=NULL
  for(j in k){
    seq= (1:len) + len*(j-1)
    dat.sel2<-dat.sel[seq,]
    dat.sel.ts2=dat.sel.ts[seq,]
    res<-simpleAGG3(data=dat.sel2,data.new=dat.sel.ts2,len=len,times=1,lam=lams)
    lik<-rbind(lik,res$lik)
    cverror<-c(cverror,res$cv.error)
  }
  
  j<-(1:length(k))
  s<-sort(c(2*j,2*j-1))
  lik<-lik[s,]
  
  p=1
  if(nrow(lik)>2){
    p=nrow(lik)/2
    seq1=seq(1,nrow(lik),by=2)
    seq2=seq(2,nrow(lik),by=2)
    res <- rbind(apply(lik[seq1,],2,sum), apply(lik[seq2,],2,sum))
  }else{res<-lik}
  
  rownames(res)<-group
  colnames(res)<-colnames(test.data)
  pred<-group[max.col(t(res))]
  return(list(pred=pred,p=p))
}

