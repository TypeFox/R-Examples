depth.RP<-function(fdataobj,fdataori=fdataobj,trim=0.25,nproj=50,proj="vexponential",
dfunc="TD1",par.dfunc=list(),scale=FALSE,draw=FALSE,...){
if (!is.fdata(fdataobj)) fdataobj=fdata(fdataobj)
if (!is.fdata(fdataori)) fdataori=fdata(fdataori)

# nas<-apply(fdataobj$data,1,count.na)
# if (any(nas))  {
#    fdataobj$data<-fdataobj$data[!nas,]
#    cat("Warning: ",sum(nas)," curves with NA are not used in the calculations \n")
#    }

 if (is.null(rownames(fdataobj$data)))  rownames(fdataobj$data)<-1:nrow(fdataobj$data)
 nms<-rownames(fdataobj$data)
 m0<-nrow(fdataobj)
 fdataobj<-na.omit.fdata(fdataobj)
 fdataori<-na.omit.fdata(fdataori)
 nas<-na.action(fdataobj)
 nullans<-!is.null(nas)


data<-fdataobj[["data"]]
data2<-fdataori[["data"]]
n<-nrow(data)
m<-ncol(data)
m2<-ncol(data2)
n2<-nrow(data2)
if (is.null(n) && is.null(m))  stop("ERROR IN THE DATA DIMENSIONS")
if (is.null(m) && is.null(m2)) stop("ERROR IN THE DATA DIMENSIONS")
t=fdataobj[["argvals"]]
rtt<-fdataobj[["rangeval"]]
names1<-names2<-names<-fdataobj[["names"]]
names1$main<-"depth.RP median"
names2$main<-paste("RP trim",trim*100,"\u0025",sep="")
 dep=rep(0.0,n)
 dep2=rep(0.0,n2)
 #### new
 if (is.fdata(proj)) {
  nproj<-nrow(proj)
  if (fdataobj$argvals!=proj$argvals || ncol(fdataobj)!=ncol(proj)) stop("Error in proj dimension")
  z<-proj
  }
else {
z<-rproc2fdata(nproj,t,sigma=proj,norm=TRUE,...)
}
##
# z<-matrix(NA,nproj,m)

# integrated<-FALSE
# if (integrated){
#   dtt <- diff(tt)
#   eps <- as.double(.Machine[[1]] * 100)
#   inf <- dtt - eps
#   sup <- dtt + eps
#   if (all(dtt > inf) & all(dtt < sup)) {
#   equi = TRUE
#   }
#   else equi = FALSE
#     valor-inprod.fdata(fdataobj,z)
#     valor2<-inprod.fdata(fdataori,z)
# }
if (dfunc %in% c("Liu1","TD1","FM1")) {
Fn<-list()
 for (j in 1:nproj){
        valor=data%*%z$data[j,]
        valor2=data2%*%z$data[j,]
        Fn[[j]]=ecdf(valor2)
         par.dfunc$x<-valor
         par.dfunc$Fn<-Fn[[j]]
         dp<-do.call(dfunc,par.dfunc)
         dep<-dep+dp
 }
  dep=dep/nproj
  }
else if (dfunc %in% c("MhD1","LD1")) {
 for (j in 1:nproj){
        valor=data%*%z$data[j,]    #HACER INPROD
        valor2=data2%*%z$data[j,]
         par.dfunc$x<-drop(valor)
         par.dfunc$xx<-drop(valor2)
         dp<-do.call(dfunc,par.dfunc)
         dep<-dep+dp
 }
  dep=dep/nproj
  }
  if (scale){    dep<-dep/max(dep)    }
   k=which.max(dep)
   med=data[k,,drop=FALSE]
   nl=length(trim)
#   mtrim=matrix(NA,nrow=nl,ncol=m)
#   for (j in 1:length(trim)) {
#                    lista=which(dep2>=quantile(dep2,probs=trim[j],na.rm=TRUE))
#                    if (length(lista)==1) {
#                      mtrim[j,]<-data2[lista,]
#                      if (draw) {draw=FALSE;warning("The plot is not shown")}
#                    }
#                    else mtrim[j,]=apply(data2[lista,],2,mean,na.rm=TRUE)
#                       }
   lista=which(dep>=quantile(dep,probs=trim,na.rm=TRUE))
   mtrim=apply(data[lista,,drop=FALSE],2,mean,na.rm=TRUE)

if (nullans) {
        ans1<-rep(NA,len=m0)
        ans1[-nas] <-dep
        dep<-ans1
        }
 names(dep)<-nms

#else mtrim=data[lista,,drop=FALSE]
   tr<-paste("RP.tr",trim*100,"\u0025",sep="")
   med<-fdata(med,t,rtt,names1)
   mtrim<-fdata(mtrim,t,rtt,names2)
   rownames(med$data)<-"RP.med"
   rownames(mtrim$data)<-tr
if (draw){
   ans<-dep
   ind1<-!is.na(ans)
   cgray=1-(ans-min(ans,na.rm=TRUE))/(max(ans,na.rm=TRUE)-min(ans,na.rm=TRUE))
   plot(fdataori, col = gray(.9),lty=1, main = "RP Depth")
   lines(fdataobj[ind1, ], col = gray(cgray[ind1]))
   lines(mtrim,lwd=2,col="yellow")
   lines(med,col="red",lwd=2)
   legend("topleft",legend=c(tr,"Median"),lwd=2,col=c("yellow","red"),box.col=0)
 }
return(invisible(list("median"=med,"lmed"=k,"mtrim"=mtrim,"ltrim"=lista,
"dep"=dep,"proj"=z,dfunc=dfunc,par.dfunc=par.dfunc)))
}

