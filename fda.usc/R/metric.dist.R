#no importa sino se divide por la suma de los pesos
euclidean<-function(u,weights=rep(1,length(u))) {  
  u<- sqrt(sum(weights*u^2,na.rm=TRUE))
  st<-attributes(u)
  attributes(u)<-st
  attr(u, "weigthts") <- weights
  u
}
manhattan<-function(u,weights=rep(1,length(u))) {  
  u<- sqrt(sum(weights*abs(u),na.rm=TRUE))
  st<-attributes(u)
  attributes(u)<-st
  attr(u, "weigthts") <- weights
  u
}
minkowski<-function(u,weights=rep(1,length(u)),p) {  
  u<- (sum(weights*(u^p),na.rm=TRUE))^(1/p)
  st<-attributes(u)
  attributes(u)<-st
  attr(u, "weigthts") <- weights
  u
}
maximum<-function(u,weights=rep(1,length(u))) {  
  u<- sqrt(weights*max(u,na.rm=TRUE))
  st<-attributes(u)
  attributes(u)<-st
  attr(u, "weigthts") <- weights
  u
}
###########################################################
###########################################################

#euclidean<-function (u)   sqrt(sum(u^2, na.rm = TRUE))
#manhattan<-function(u) sqrt(sum(abs(u),na.rm=TRUE))
#minkowski<-function(u,p) (sum(u,na.rm=TRUE)^p)^1/p
#no hace falta w
#maximum<-function(u) sqrt(max(u,na.rm=TRUE))
################################################################################
metric.ldata=function(lfdata,lfdataref=lfdata,metric,par.metric=list(),weights,method="euclidean") {
if (class(lfdata)=="list"){
 lenl<-length(lfdata)
 lenl2<-length(lfdataref)
 n<-nrow(lfdata[[1]])
 m<-nrow(lfdataref[[1]]) 
 mdist2<-matrix(0,n,n)
 amdist<-array(NA,dim=c(n,m,lenl))
 ldist<-mdist<-list()
 nam1<-names(lfdata)
 nam2<-names(lfdataref) 
 attr<-list()
 if (is.null(nam1)) {names(lfdata)<-nam1<-paste("var",1:lenl,sep="")}
 if (is.null(nam2)) {names(lfdataref)<-nam2<-paste("var",1:lenl2,sep="")} 
 # verificar que todos los nam1 estan en nam2
 # verificar que longitud de metric y par.metric es la misma q ldata
 if (missing(metric)){
   if (is.fdata(lfdata[[nam1[1]]])) metric<-rep("metric.lp",len=lenl)
   else    metric<-rep("metric.dist",len=lenl)   
 }
if (missing(weights)) weights<-rep(1,length=lenl)
# print("metric.ldata2")
 for (i in 1:lenl){
# print("metric.ldata3")
if (is.character(metric)) {
 if (is.null(par.metric[[nam1[i]]])) par.metric[[nam1[i]]]<-list()
  if (is.fdata(lfdata[[nam1[i]]])) { 
   par.metric[[nam1[i]]][["fdata1"]]<-lfdata[[nam1[i]]]
   par.metric[[nam1[i]]][["fdata2"]]<-lfdataref[[nam1[i]]]   
    
   mdist<-do.call(metric[i],par.metric[[nam1[i]]])
   }
  else {
   par.metric[[nam1[i]]][["x"]]<-lfdata[[nam1[i]]]
   par.metric[[nam1[i]]][["y"]]<-lfdataref[[nam1[i]]]   
   mdist<-do.call(metric[i],par.metric[[nam1[i]]])
   }
}
else {
  
  if (is.list(metric)) mdist<-metric[[nam1[i]]]
  
}
  ldist[[nam1[i]]]<-mdist
  amdist[,,i]<- mdist
  #attr[[nam1[i]]]<-attributes(mdist)     
 }
# print("sale1")
#print(Weights)
###for (k in len) sqrt(w[i])*
# print(weights)
mdist<-apply(amdist,1:2,method,weights=weights)
# print("sale2")
}
else stop("Error in lfdata argument")
# print("a")
# print(dim(amdist))
# print(dim(mdist))
# print(attributes(mdist))
# print(attr)
#attributes(mdist)<-atr
attr(mdist, "method") <-method
attr(mdist, "weights") <- weights
for (i in 1:lenl) attr(mdist, nam1[i]) <- attributes(ldist[[nam1[i]]])
return(mdist)
}


################################################################################
metric.dist<-function(x,y=NULL,method="euclidean",p=2,dscale=1,...){
#print("*******************************************************entra meeeeetric.dist")
if (is.vector(x)) x<-matrix(x,nrow=1)
else x<-as.matrix(x)
ynull<-is.null(y)
if (method=="mahalanobis"){
    if (ynull)   {
        y <- x
        vc <- var(x)
        }
    else {
    y <- as.matrix(y)
    vc <- var(rbind(x, y))
    }
    n <- nrow(x)
    m <- nrow(y)
    mdist<- matrix(0, n, m)
    for (i in 1:m) {
        mdist[, i] <- mahalanobis(x, y[i, ], cov = vc)
    } 
    mdist<-sqrt(mdist)
}
else{
 if (!ynull) {    
    if (is.vector(y)) y<-matrix(y,nrow=1) 
    n<-nrow(y)
    nn<-nrow(x)
    mdist<-as.matrix(dist(rbind(x,y) , method = method, diag = TRUE, upper = TRUE,p=p))[1:nn,(nn+1):(nn+n)] 
    }
 else   mdist<-as.matrix(dist(x, method = method, diag = TRUE, upper = TRUE,p=p))  
 }
 if (is.vector(mdist)) mdist<-matrix(mdist,nrow=nn)   
 	etiq1=rownames(x)
	etiq2=rownames(y)
# namesx<-rownames(x)
# if (ynull) dimnames(mdist) <- list(namesx,namesx)
# else dimnames(mdist) <- list(namesx, rownames(y))
 if (is.function(dscale)) {
   if (nrow(mdist)==ncol(mdist)) diag(mdist)<-NA ################# ojjjooo solo para matrices cuadradas
   dscale<-dscale(as.dist(mdist))
   if (nrow(mdist)==ncol(mdist)) diag(mdist)<-0   
 } 
 mdist<-mdist/dscale
 attr(mdist, "call") <- "metric.dist"
 attr(mdist, "par.metric") <- list(method =method,p=p,dscale=dscale) 
 rownames(mdist)<-etiq1
 colnames(mdist)<-etiq2
 return(mdist)
}

