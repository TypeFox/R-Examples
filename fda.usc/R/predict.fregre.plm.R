predict.fregre.plm<-function(object,newx=NULL,...){
 if (is.null(object)) stop("No fregre.plm object entered")
 if (is.null(newx)) {
    cat("No newx entered \n")
    newx=object$fdataobj # los XX
    }
 data=newx
 formula=object$formula
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
 vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
# vnf=c(vnf2,vint) revisar sino no va bien vnf
 off<-attr(tf,"offset")
 kterms=1
 z=list()
 lenvnf=length(vnf)
 if (lenvnf>0) {
# cat(" Non functional variables: ",vnf,"\n")
 for ( i in 1:length(vnf)){
    if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
 if   (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
 }
 y=object$y
 n=nrow(y)
 if (length(vfunc)>0) {
  if (!is.fdata(data[[vfunc[1]]])) fdataobj=fdata(data[[vfunc[1]]],newx[["argvals"]],newx[["rangeval"]])
  else fdataobj=data[[vfunc[1]]]
  x.fd<-fdataobj[["data"]]
#  if (class(data[[vfunc[1]]])[1]=="fd")   	 x.fd=t(data[[vfunc[1]]]$coefs)
#  else    	 x.fd=data[[vfunc[1]]]
#  if (is.data.frame(x.fd)) x.fd=as.matrix(x.fd)
#  cat(" ",class(data[[vfunc[1]]])[1]," object: ",vfunc[1],"\n")
  I=diag(1,ncol=nrow(x.fd),nrow=nrow(x.fd))
  XX=as.matrix(data[["df"]][,vnf])
  colnames(XX)=vnf2
   xd=object$fdataobj
   x=object$XX
   h=object$h.opt
   betah=object$betah
   par.metric<-attr(object$mdist,"par.metric")
   par.metric[["fdata2"]]<-xd
   par.metric[["fdata1"]]<-fdataobj
#   parm<-attr(object$mdist,"par.metric")
#par.metric[3:(2+length(parm))]<-attr(object$mdist,"par.metric")
#names(par.metric[3:(2+length(parm))])<-names(attr(object$mdist,"par.metric")[3:(2+length(parm))])
   a1<-attr(object$mdist,"call")
   nmdist <- do.call(a1,par.metric)

#    xh=wh%*%object$XX
    nw=object$type.S(nmdist,h,object$Ker,cv=FALSE)
    nmh=nw%*%(y-x%*%betah)
#    mh=ww%*%(Y[j,]-XX%*%betah)
    yp=XX%*%betah+nmh
}
else {

   XX=data[[1]][,c(vnf2)]
   cat("Warning: predict.lm  done, non functional data in the formula \n")
#   if (is.null(newx)) yp=predict.lm(object,...)
#   else  {
#         if (!is.data.frame(XX)) XX=data.frame(XX)
#           print(names(XX))
#           yp=predict.lm(object,XX,...)
#          }
  }
}
else {
 cat("Warning: predict fregre.np object, only functional data in the formula \n")
# print(names(object))
# yp=predict.fregre.fd(object,data[[vfunc[1]]],...)
 }
 yp
}

