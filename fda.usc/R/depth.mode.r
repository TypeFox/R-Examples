################################################################################
# dscale puede ser una funcion o un valor
################################################################################
depth.modep=function(lfdata,lfdataref=lfdata,h=NULL,metric,
                     par.metric=list(),method="euclidean",
                     scale=FALSE,trim=0.25,draw=FALSE,ask=FALSE) 
{
 equal<-identical(lfdata,lfdataref)
 if (class(lfdata)=="list"){
  lenl <- length(lfdata)
  lenl2 <- length(lfdataref)
  m0 <- nrow(lfdata[[1]])
  if (is.null(rownames(lfdata[[1]]$data))) 
    rownames(lfdata[[1]]$data) <- 1:m0
  nms <- rownames(lfdata[[1]]$data)
  nas <- NULL
  for (i in 1:lenl) {
    nas <- c(nas, na.action(na.omit(lfdata[[i]])))
  }
  nas <- unique(nas)
  nullans <- !is.null(nas)
  nam1<-names(lfdata)
  if (missing(metric)){
    if (is.fdata(lfdata[[nam1[1]]])) metric<-rep("metric.lp",len=lenl)
    else    metric<-rep("metric.dist",len=lenl)   
  }
  mdist2<-metric.ldata(lfdata=lfdata,lfdataref=lfdataref,metric=metric,par.metric=par.metric,method=method)
  if (equal) mdist<-mdist2
  else mdist<-metric.ldata(lfdata=lfdataref,lfdataref=lfdataref,metric=metric,par.metric=par.metric,method=method)
  m<-m0
  n<-nrow(lfdata[[1]])
}
else stop("Error in lfdata argument")
#attr(mdist, "par.metric") <-attr(atr, "par.metric")                                  
#attr(mdist, "call") <-attr(atr, "call")                                  
#attr(mdist, "method") <-method                                    
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
#if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
#h="h=0.15"
if (is.null(h))   {
  h<-0.1
  hq2=quantile(mdist,probs=h,na.rm=TRUE)
#print("es nulo h")  
}
else {
  #cat("no es nulo h ",h,"\n")    
  if (is.numeric(h))    hq2<-h  
  else hq2=quantile(mdist,probs=as.numeric(substring(h,first=3)),na.rm=TRUE)
}
#print(hq2)
class(mdist2)<-class(mdist)<-c("matrix","fdist")
ans<-Ker.norm(mdist2/hq2)    #### en (nxm) matrix (new X ref)
#print(ans)
ans<-apply(ans,1,sum,na.rm=TRUE)                                    
#print(ans)
#print(ans[1:3])

if (scale)   {
   mn<-min(ans,na.rm=TRUE)
   mx<-max(ans,na.rm=TRUE)
#   scl<-mx-mn
#   ans=as.vector(scale(ans,center=mn,scale=scl))
#   ans2=as.vector(scale(ans2,center=mn,scale=scl))
   ans=as.vector(ans/mx)   
}
 k=which.max(ans)
 lista=which(ans>=quantile(ans,probs=trim,na.rm=TRUE))
if (nullans) {
 ans1<-rep(NA,len=m0)
 ans1[-nas] <-ans 
 ans<-ans1    
 }       
names(ans)<-nms  
if (draw){
  mf=5
  if (lenl>4) ask=TRUE
  if (ask) {par(mfrow = c(1, 1))
            dev.interactive()
            oask <- devAskNewPage(TRUE)
            on.exit(devAskNewPage(oask))}
   else{    mf<-switch(lenl,
   "1"={c(1,1)},
   "2"={c(1,2)},
   "3"={c(1,3)},
   "4"={c(2,2)})            
    par(mfrow =mf)                    
   }
 names1<-names2<-names<-lfdata[[1]][["names"]]
 names1$main<-"depth.modep median"
 tr<-paste("mode.tr",trim*100,"\u0025",sep="")                     
 for (idat in 1:lenl) {
   data<-lfdata[[idat]]$data
   tt<-lfdata[[idat]]$argvals
   rtt<-lfdata[[idat]]$rangeval
   med<-data[k,]          
   if (n>1) mtrim=apply(data[lista,],2,mean,na.rm=TRUE)
   else mtrim=data[lista,]        
   med<-fdata(med,tt,rtt,names1)
   mtrim<-fdata(mtrim,tt,rtt,names2)
   rownames(med$data)<-"modep.med"
   rownames(mtrim$data)<-tr
   ind1<-!is.nan(ans)
   ans[is.nan(ans)]=NA     
   cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
   plot(lfdata[[idat]][ind1, ], col =  gray(cgray[ind1]),lty=1, main = paste(nam1[idat]," modep Depth",sep=""))
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
}
# print("fin modep")
return(invisible(list("lmed"=k,"ltrim"=lista,"dep"=ans,"metric"=metric,
"par.metric"=par.metric,"mdist"=mdist,"hq"=hq2)))
} 






























##################################################################################
depth.mode=function(fdataobj,fdataori=fdataobj,trim=0.25,metric=metric.lp,h=NULL,scale=FALSE,
draw=FALSE,...){    
if (is.fdata(fdataobj)) {
 fdat<-TRUE
# nas<-apply(fdataobj$data,1,count.na)
#if (any(nas))  {
#   fdataobj$data<-fdataobj$data[!nas,]
#   cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
#   }
 if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
 fdataobj<-na.omit.fdata(fdataobj)
 fdataori<-na.omit.fdata(fdataori) 
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas) 
 data<-fdataobj[["data"]]
 data2<-fdataori[["data"]]
 names1<-names2<-names<-fdataobj[["names"]]
 names1$main<-"depth.mode median"
 names2$main<-paste("depth.mode trim ",trim*100,"\u0025",sep="")
 tt=fdataobj[["argvals"]]
 rtt<-fdataobj[["rangeval"]]
#print("is fdata") 
}
else { stop("no fdata class object")
        data<-fdataobj
        data2<-fdataori
        fdat<-FALSE   
  }
n<-nrow(data)
m<-ncol(data)
m2<-ncol(data2)
n2<-nrow(data2)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.matrix(metric)) mdist=metric                    
else  mdist=metric(fdataori,fdataori,...)
class(mdist)<-"matrix"

if (n==n2 & m==m2) {
  equal<-all(fdataobj$data==fdataori$data)
  if (equal) mdist2<-mdist
  else mdist2<-metric(fdataobj,fdataori,...)}
else  mdist2<-metric(fdataobj,fdataori,...)
if (is.null(h))   {
  h<-0.15
  hq2=quantile(mdist,probs=h,na.rm=TRUE)
  }
else {
  if (is.numeric(h))    hq2<-h  
  else hq2=quantile(mdist,probs=as.numeric(substring(h,first=3)),na.rm=TRUE)
}
class(mdist)<-class(mdist2)<-c("matrix","fdist")
ans<-Ker.norm(mdist2/hq2)    ####
ans<-apply(ans,1,sum,na.rm=TRUE)                                    
if (scale)   {
#   mdist2<-metric(fdataori,fdataori,...)
#   ans2<-Ker.norm(mdist2/hq2)    ####
#   ans2<-apply(ans2,1,sum,na.rm=TRUE)                                       
   mn<-min(ans,na.rm=TRUE)
   mx<-max(ans,na.rm=TRUE)
#   scl<-mx-mn
#   ans=as.vector(scale(ans,center=mn,scale=scl))
#   ans2=as.vector(scale(ans2,center=mn,scale=scl))
   ans=as.vector(ans/mx)   
}    
k=which.max(ans)
med=data[k,]
lista=which(ans>=quantile(ans,probs=trim,na.rm=T))
if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-ans 
        ans<-ans1      
        }
 names(ans)<-nms   
if (length(lista)==1) {
  mtrim<-data[lista,]
  if (draw) {draw=FALSE;warning("The plot is not shown")}
  }
else mtrim=apply(fdataobj[lista]$data,2,mean,na.rm=TRUE)
tr<-paste("mode.tr",trim*100,"\u0025",sep="")
if (fdat) {
med<-fdata(med,tt,rtt,names1)
mtrim<-fdata(mtrim,tt,rtt,names2)
rownames(med$data)<-"mode.med"
rownames(mtrim$data)<-tr
if (draw){
    if (!scale){
     mn<-min(ans,na.rm=TRUE)
     mx<-max(ans,na.rm=TRUE)
     scl<-mx-mn
    }     
   ind1<-!is.nan(ans)
   ans[is.nan(ans)]=NA
   cgray=1-(ans-mn)/(scl)
   plot(fdataori, col = gray(.9),lty=1, main = "mode Depth")
   lines(fdataobj[ind1, ], col = gray(cgray[ind1]))
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
 }
out<-list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,
"dep"=ans,"hq"=hq2) 
if (scale) out$dscale=mx
return(invisible(out))
}
################################################################################
mdepth.LD=function(x,xx=x,metric=metric.dist,h=NULL,scale=FALSE,...){    
 if (is.vector(x)){
   if (all(xx==x)) x<-xx<-matrix(x,ncol=1)#stop("One of x or xx must be a matrix")
   else   {
 	m2=ncol(xx)
	if (length(x)!=m2) stop("Length of x does not match with dimension of xx")
	x=matrix(x,ncol=m2)
   }
 }    
		m2=ncol(xx)
		n2=nrow(xx)
		n<-nrow(x)
		m<-ncol(x)
if (is.null(rownames(x)))  rownames(x)<-1:nrow(x)
 nms<-rownames(x)
 x<-na.omit(x)  
 xx<-na.omit(xx)
 nas<-na.action(x)
 nullans<-!is.null(nas) 
        d <- ncol(x)
if (is.null(n) && is.null(m)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
if (is.matrix(metric)) {mdist=metric}
else {  mdist=metric(xx,xx,...)  }
class(mdist)<-"matrix"

if (n==n2 & m==m2) {
  equal<-all(x==xx)
  if (equal) mdist2<-mdist
  else mdist2<-metric(x,xx,...)}
else  mdist2<-metric(x,xx,...)
if (is.null(h))   {
  h<-0.15
  hq2=quantile(mdist,probs=h,na.rm=TRUE)
}
else {
  if (is.numeric(h))    hq2<-h  
  else hq2=quantile(mdist,probs=as.numeric(substring(h,first=3)),na.rm=TRUE)  
}
class(mdist)<-class(mdist2)<-c("matrix")
ans<-Ker.norm(mdist2/hq2)    ####
ans<-apply(ans,1,sum,na.rm=TRUE)                                    
mx = scale
if (scale)   {
#  mn<-min(ans,na.rm=TRUE)
  mx<-max(ans,na.rm=TRUE)
  ans=as.vector(ans/mx)   
}                                
if  (nullans){
        ans1<-rep(NA,len=n)
        ans1[-nas] <-ans 
        ans<-ans1      
        }
names(ans)<-nms   
out<-list("dep" = ans,"hq"=hq2,dscale=mx)                                                                  
return(invisible(out))
}
##########################################################################
##########################################################################
#out.p[[icls]]=classif.DD(yy,lxx,depth=dep.int[idep],classif=cls[icls])                  

#d<-depth.modep(lxx,lx.test)

