#' @export
bag <- function(formula,data,status,nBoot=10, minbucket = 40){
  
  if(is.factor(data$status)){
    stop("status variable must be numeric")
  }
  
  n <- nrow(data)
  pred.m <- pred.s <- pred.h <- vector("list",nBoot)
  mf <- model.frame(formula,data,na.action=NULL)
  uniq <- unique(model.response(mf))
  m <- length(uniq)
  trees <- vector("list", length = nBoot)
  
  for (i in 1:nBoot){
    
  samp<-sample.ds(data, all.vars(formula(mf)[1]))
  trees[[i]]<- DStree(formula,data=samp,status=status,control=list(minbucket=minbucket,cp=0)) 
  pred<-predict.DStree(trees[[i]],data)
  
  mf.b<- model.frame(formula,samp,na.action=NULL)
  tab.b <- table(factor(model.response(mf.b),levels=uniq))
  ind.NA <- which(tab.b==0)
  
  
  if(length(ind.NA)==0){
    pred.h[[i]] <- pred$Haz
  }else{
    pred.h[[i]] <- matrix(0,n,m)
    pred.h[[i]][,-ind.NA] <-  pred$Haz
    pred.h[[i]][,ind.NA] <-  NA
    
  }
  

  }
  

  Haz <- apply(simplify2array(pred.h), c(1,2),function(x) mean(x,na.rm=T))
  Surv1 <- matrix(NA,n,m+1)
  Surv1[,1] <- 1
  for(i in 1:m){
    Surv1[,i+1] <- Surv1[,i]*(1-Haz[,i])
  }
  Surv <- Surv1[,-1]
  rownames(Surv) <- rownames(Haz)
  
  Med1 <- apply(Surv, 1, function(x) rev(x[x>=0.5])[1])
  Med2 <- apply(Surv, 1, function(x) x[x<0.5][1])
  
  MedSurv<-rep(NA,n)
  m<-rep(NA,n)
  
  for (i in 1:n){
    
    if(is.na(Med1[i])){
      Med1[i] <- 1
      m[i] <- 0
      }else{
      ms <- which(Surv[i,]==Med1[i])
      m[i] <- ms[length(ms)]
      }
    
    if(is.na(Med2[i])){
      MedSurv[i] <- NA
    }else{
      MedSurv[i] <- round(m[i]+(Med1[i]-0.5)/Med1[i]-Med2[i],1)
    }
      
  }
  ret<-list(MedSurv=MedSurv,Surv=Surv,Haz=Haz,minbucket=minbucket,nBoot=nBoot,trees=trees)
  class(ret)<-"DStreebag"
  return(ret)
}


