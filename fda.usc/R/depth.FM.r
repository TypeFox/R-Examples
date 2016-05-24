################################################################################
################################################################################
# Univariate version
mdepth.FM1=function(x,xx=x,scale=FALSE){
isx<-is.vector(x)
isxx<-is.vector(xx)
if (!isx)  stop("x object is not a vector")
if (!isxx)  stop("xx object is not a vector")
n<-length(x)
m<-length(xx)
if (n!=m)  stop("The length of x is not the same of xx")
Fn<-list()
Fn=ecdf(xx)
d=1-abs(0.5-Fn(x))
if (scale){     d<-(d-.5)*2    }   
return(list("dep"=d,"Fn"=Fn))  #
}
################################################################################
################################################################################
# univariate version
mdepth.TD1=function(x,xx=x,xeps=1e-15,scale=FALSE){
#print("entra TD")
isx<-is.vector(x)
isxx<-is.vector(xx)
if (!isx)  stop("x object is not a vector")
if (!isxx)  stop("xx object is not a vector")  
n<-length(x)
m<-length(xx)
if (n!=m)  stop("The length of x is not the same of xx")
Fn=ecdf(xx)
d=pmin(Fn(x),(1-Fn(x-xeps)))
if (scale){ d<-d*2      }
return(list("dep"=d,"Fn"=Fn))  #
}
################################################################################
################################################################################
# univariate version
mdepth.Liu1=function(x,xx=x,xeps=1e-15,scale=FALSE){
isx<-is.vector(x)
isxx<-is.vector(xx)
if (!isx)  stop("x object is not a vector")
if (!isxx)  stop("xx object is not a vector")
n<-length(x)
m<-length(xx)
if (n!=m)  stop("The length of x is not the same of xx")
Fn=ecdf(xx)
d=Fn(x)*(1-Fn(x-xeps))
if (scale){  d<-d*4      }
return(list("dep"=d,"Fn"=Fn))  
}

################################################################################
# multivariate version
mdepth.FM=function(x,xx=x,scale=FALSE,dfunc="mdepth.TD1"){
n<-nrow(x)
m<-ncol(xx)
d<-matrix(NA,nrow=n,ncol=m)
Fn<-list()
for (i in 1:m)   {
    Fn[[i]]=ecdf(xx[,i])
    d[,i]=1-abs(0.5-Fn(x[,i])) 
    }
d<-apply(d,1,mean,na.rm=TRUE)[1:n]
return(list("dep"=d,"Fn"=Fn))  #
}
################################################################################
################################################################################
mdepth.Liu=function(x,xx=x,xeps=1e-15,scale=FALSE){
n<-nrow(x)
m<-ncol(x)
d<-matrix(NA,nrow=n,ncol=m)
Fn<-list() 
for (i in 1:m){
         Fn[[i]]=ecdf(xx[,i])
         d[,i]=Fn[[i]](x[,i])*(1-Fn[[i]](x[,i]-xeps))
    }
if (scale){   d<-d*4   }
d<-apply(d,1,mean,na.rm=TRUE)[1:n]
return(list("dep"=d,"Fn"=Fn))  
}     
################################################################################
################################################################################
# multivariate version
mdepth.TD=function(x,xx=x,xeps=1e-15,scale=FALSE){
#print("entra TD")
n<-nrow(x)
m<-ncol(x)
d<-matrix(NA,nrow=n,ncol=m)
Fn<-list() 
for (i in 1:m){
         Fn[[i]]=ecdf(xx[,i])
         d[,i]=pmin(Fn[[i]](x[,i]),(1-Fn[[i]](x[,i]-xeps)))
    }
if (scale){ d<-d*2      }
d<-apply(d,1,mean,na.rm=TRUE)[1:n]
return(list("dep"=d,"Fn"=Fn)) #"dep.ori"=d2,
}
################################################################################
FM1<-function(x,Fn,scale=FALSE) {
 d=1-abs(0.5-Fn(x))
 if (scale){     d<-(d-.5)*2    }  
  d
 }
Liu1<-function(x,Fn,xeps=1e-15,scale=FALSE) {
  d=Fn(x)*(1-Fn(x-xeps))
 if (scale){   d<-d*4   }
 d
}
TD1<-function(x,Fn,xeps=1e-15,scale=FALSE) {
 d=pmin(Fn(x),(1-Fn(x-xeps)))      #HS1D
 if (scale){ d<-d*2      }
 d
}
LD1<-function(x,xx=x,scale=TRUE) {
 mdepth.LD(matrix(x,ncol=1),matrix(xx,ncol=1),scale=scale)$dep
}      
MhD1 <- function(x,xx=x,scale=FALSE){
 if (is.vector(x)) {
              D<- (1+(x-(mean(xx)))^2/sd(xx)^2)
              }
 else{stop("x is no a vector") }
   ans <- 1/D
ans
}   

 
################################################################################
depth.FM=function(fdataobj,fdataori=fdataobj,trim=0.25,scale=FALSE,dfunc="FM1",par.dfunc=list(scale=TRUE),draw=FALSE){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (!is.fdata(fdataori)) fdataori=fdata(fdataori)
if (is.null(fdataobj))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
nms<-rownames(fdataobj$data)

#nas<-apply(fd-ataobj$data,1,count.na)
#if (any(nas))  {
#   fdataobj$data<-fdataobj$data[!nas,]
#   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
#   }             
data<-fdataobj[["data"]]
data2<-fdataori[["data"]]
n<-nrow(data)
m<-ncol(data)
m2<-ncol(data2)
n2<-nrow(data2)
d<-matrix(NA,nrow=n,ncol=m)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.list(dfunc)) dfunc<-dfunc$dfunc
t=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depth.FM median"
Fn<-list()
tr<-paste("FM.tr",trim*100,"\u0025",sep="")
  dtt <- diff(t)
  eps <- as.double(.Machine[[1]] * 100)
  inf <- dtt - eps
  sup <- dtt + eps
  if (all(dtt > inf) & all(dtt < sup)) {  equi = TRUE  }
  else equi = FALSE
for (i in 1:m)   {
if (dfunc %in% c("TD1","Liu1","FM1")){
         Fn[[i]]=ecdf(data2[,i])
		 par.dfunc$x=data[,i]
		 par.dfunc$Fn=Fn[[i]]
         d[,i]=do.call(dfunc,par.dfunc)
         }
else     {
		par.dfunc$x=data[,i]
		par.dfunc$xx=data2[,i]
			d[,i]=do.call(dfunc,par.dfunc)
		 }
    }
#d<-int.simpson(fdata(d,t,rtt),equi=equi)
d<-apply(d,1,mean)[1:n]
ans<-d
if (scale) { ans=d/max(d)}
names(ans)<-nms
k=which.max(ans)    
med=data[k,]                                                      
lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
if (n>1)  mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)
else mtrim<-data                                     
med<-fdata(med,t,rtt,names1)
mtrim<-fdata(mtrim,t,rtt,names2)
rownames(med$data)<-"FM.med"
rownames(mtrim$data)<-tr
if (draw){
   ind1<-!is.nan(ans)
   ans[is.nan(ans)]=NA
   cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
   plot(fdataobj[ind1, ], col =  gray(cgray[ind1]),lty=1, main = "FM Depth")
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),
   box.col=0)
 }
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,
"dep"=ans,"Fn"=Fn)))
} 

  
#####################################################
################################################################################
depth.FMp=function(lfdata,lfdataref=lfdata,trim=0.25,dfunc="mdepth.MhD",
par.dfunc=list(scale=FALSE),draw=FALSE,ask=FALSE,...){
  #   print("FMp")
#if (!is.list(lfdata)) stop("lfdata1 must be a list")
#if (!is.list(lfdataref)) stop("lfdata1 must be a list")
nam1<-names(lfdata)
nam2<-names(lfdata)
len1<-length(lfdata)
len2<-length(lfdataref)   
  m0<-nrow(lfdata[[1]])
  if (is.null(rownames(lfdata[[1]]$data)))  rownames(lfdata[[1]]$data)<-1:m0
 nms<-rownames(lfdata[[1]]$data)
 nas<-NULL
for (i in 1:len1) {
  nas<-c(nas,na.action(na.omit(lfdata[[i]])))
} 
nas<-  (nas)
nullans<-!is.null(nas) 
for (i in 1:len1) {
#  if (nullans) lfdata[[i]]<-lfdata[[i]][-nas]
#  if (nullans)  lfdataref[[i]]<-lfdataref[[i]][-nas]  
}
#comprovar is.fdatalist
if (is.null(nam1)) nam1<-paste("var",1:len1,sep="")
if (is.null(nam2)) nam2<-nam1
 fdataobj<-lfdata[[1]]
 fdataori<-lfdataref[[1]] 
 if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
n<-nrow(fdataobj)
m<-nrow(fdataori)
p<-ncol(fdataori)
tt<-fdataori$argvals
rtt<-fdataori$rangeval
xref<-matrix(NA,m,len1)
x0<-matrix(NA,n,len1)
#depth<-paste("mdepth.",dfunc,sep="")
d<-array(NA,dim=c(n,p)) 
for (idat in 1:len1) {
   if (any(lfdata[[idat]]$argvals!=tt))  stop("Incorrect argvals in the fdata objects")}
# verificar dimensiones lfdata y mismos argvals/rangeval
for (iti in 1:p) {
#  lfdata[[idat]]
  for (idat in 1:len1) { 
   x0[,idat]<-lfdata[[idat]]$data[,iti]  
   xref[,idat]<-lfdataref[[idat]]$data[,iti]    
  } 
  par.dfunc$x<-x0
  par.dfunc$xx<-xref
  d[,iti]<-do.call(dfunc,par.dfunc)$dep 
 }   
#print(d)         
data<-fdataobj[["data"]]
data2<-fdataori[["data"]]
n<-nrow(data)
m<-ncol(data)
m2<-ncol(data2)
n2<-nrow(data2)
names1<-names2<-names<-lfdata[[1]][["names"]]
names1$main<-"depth.FM median"
 tr<-paste("FM.tr",trim*100,"\u0025",sep="")
  dtt <- diff(tt)
  eps <- as.double(.Machine[[1]] * 100)
  inf <- dtt - eps
  sup <- dtt + eps
  if (all(dtt > inf) & all(dtt < sup)) {
  equi = TRUE
  }
  else equi = FALSE
#d2<-int.simpson(fdata(d,tt,rtt),equi=equi)/n
ans<-apply(d,1,mean)[1:n]
#if (nullans) {
#  ans1<-rep(NA,len=m0)
 #ans1[-nas] <-ans 
# ans<-ans1    
# }       
names(ans)<-nms    
k=which.max(ans)    
lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
if (draw){
  mf=5
  if (len1>4) ask=TRUE
  if (ask) {par(mfrow = c(1, 1))
            dev.interactive()
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))}
   else{    mf<-switch(len1,
   "1"={c(1,1)},
   "2"={c(1,2)},
   "3"={c(1,3)},
   "4"={c(2,2)})            
            par(mfrow =mf)                    }          
 for (idat in 1:len1) {
   data<-lfdata[[idat]]$data
   med<-data[k,]    
  
   mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)
   med<-fdata(med,tt,rtt,names1)
   mtrim<-fdata(mtrim,tt,rtt,names2)
   rownames(med$data)<-"FMp.med"
   rownames(mtrim$data)<-tr
   ind1<-!is.nan(ans)
   ans[is.nan(ans)]=NA
   cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
   plot(lfdata[[idat]][ind1, ], col =  gray(cgray[ind1]),lty=1, main = paste(nam1[idat]," FMp Depth",sep=""))
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
}
# print("sale FMp")
return(invisible(list("lmed"=k,"ltrim"=lista,"dep"=ans,"par.dfunc"=par.dfunc)))
}
