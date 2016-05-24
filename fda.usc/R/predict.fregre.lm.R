predict.fregre.lm<-function(object,newx=NULL,type="response",se.fit=FALSE,scale = NULL,df=df,
    interval = "none", level = 0.95,weights = 1, pred.var = res.var/weights, ...){
 if (is.null(object)) stop("No fregre.lm object entered")
 if (is.null(newx)) {
    yp=predict.lm(object,type=type,se.fit=se.fit,interval=interval,level=level,weights=weights,pred.var=pred.var,df=df,scale=scale,...)
    print("No newx entered")
    return(yp)
    }
 else {
 data=newx
 basis.x=object$basis.x
 basis.b=object$basis.b
 formula=object$formula.ini
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
# vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
 #vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 beta.l=list()
 kterms=1
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 if   (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
if (length(vnf)>0) {
 first=FALSE
# print(paste("no functional variable",vnf[i]))
   for ( i in 1:length(vnf)){
#     print(paste("Non functional covariate:",vnf[i]))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
#     else pf <- paste(pf, terms[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (attr(tf,"intercept")==0) {
#     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
   mf<-as.data.frame(model.matrix(formula(pf),data$df))
   vnf2<-names(mf)[-1]   
    for ( i in 1:length(vnf2))  pf<-paste(pf, "+", vnf2[i], sep = "")
    XX <- mf 
}
else {
    pf2<-paste(pf, "1", sep = "")
    XX <- data.frame(model.matrix(formula(pf2),data$df))
    first=TRUE
}
if (length(vnf)>0) {
  spm<-matrix(object$coefficients[names(XX)],ncol=1)
 yp<-as.matrix(XX)%*%spm
 
 }
else yp<-object$coefficients[1]*rep(1,len=nrow(newx[[vfunc[1]]])) 
if (length(vfunc)>0)  {
#   yp2<-a1 <- object$coefficients[1] * rep(1, len = nrow(data[[vfunc[1]]]))
   for (i in 1:length(vfunc)) {

   if(class(data[[vfunc[i]]])[1]=="fdata")  {
     fdataobj<-data[[vfunc[i]]]
      x.fd<-fdataobj[["data"]]
      tt<-fdataobj[["argvals"]]
      rtt<-fdataobj[["rangeval"]]
      if (!object$basis.x[[vfunc[i]]]$type=="pc"&!object$basis.x[[vfunc[i]]]$type=="pls") { 
 	      	x.fd = Data2fd(argvals = tt, y = t(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]]$data),
                      basisobj = basis.x[[vfunc[i]]],fdnames=rownames(x.fd))
	    	  r=x.fd[[2]][[3]]
          J<-object$JJ[[vfunc[i]]]
          Z = t(x.fd$coefs) %*% J
          colnames(Z) = colnames(J)
      }
      else { 
          name.coef<-paste(vfunc[i], ".",rownames(object$basis.x[[vfunc[i]]]$basis$data),sep ="")
          newXcen<-fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]]                  
                      if (object$basis.x[[vfunc[i]]]$type == "pls") {
                       if (object$basis.x[[vfunc[i]]]$norm)  {
                         sd.X <- sqrt(apply(object$basis.x[[vfunc[i]]]$fdataobj$data, 2, var))
                         newXcen$data<- newXcen$data/(rep(1, nrow(newXcen)) %*% t(sd.X))
                        }
                      } 
                    Z<- inprod.fdata(newXcen,object$vs.list[[vfunc[i]]]) 
                    
#          Z<- inprod.fdata(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]],object$vs.list[[vfunc[i]]])
          colnames(Z)<-name.coef
#         object$beta.l[[vfunc[i]]]$data <- matrix(object$beta.l[[vfunc[i]]]$data,nrow = 1)
#         b1 <- inprod.fdata(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]],object$beta.l[[vfunc[i]]])
#         yp2<-yp2+b1
      }
       if (first) {    XX=Z;              first=FALSE         }
       else XX = cbind(XX, Z)
      }
      else {
          if(class(data[[vfunc[i]]])[1]=="fd")  {
             if (class(object$basis.x[[vfunc[i]]])!="pca.fd") {
             x.fd<-fdataobj<-data[[vfunc[i]]]
 	    	     r=x.fd[[2]][[3]]
             J<-object$JJ[[vfunc[i]]]
             x.fd$coefs<-x.fd$coefs-object$mean[[vfunc[i]]]$coefs[,1]
             Z = t(x.fd$coefs) %*% J
             colnames(Z) = colnames(J)
             }
             else {
                       name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(object$basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
                       data[[vfunc[i]]]$coefs<- sweep(data[[vfunc[i]]]$coefs,1,(object$basis.x[[vfunc[i]]]$meanfd$coefs),FUN="-")
                       fd.cen<-data[[vfunc[i]]]
                       #          fd.cen<-data[[vfunc[i]]]-object$basis.x[[vfunc[i]]]$meanfd # de las CP basi
                       Z<- inprod(fd.cen,object$basis.x[[vfunc[i]]]$harmonics)
                       colnames(Z)<-name.coef[[vfunc[i]]]
                   }
              if (first) {    XX=Z;              first=FALSE         }
             else XX = cbind(XX, Z)
           }
          else stop("Please, enter functional covariate")
       }  }
       }
 nn<-nrow(XX)  
 if (!is.data.frame(XX)) XX=data.frame(XX)
    if (!object$rn)   return(predict.lm(object=object,newdata=XX,type=type,
    se.fit=se.fit,interval=interval,level=level,weights=weights,pred.var=pred.var,df=df,scale=scale,,...))  #cambiarlo
else {  
  for (i in 1:length(vfunc)){
  if (object$call[[1]]=="fregre.pls")  return(predict.lm(object=object,newdata=XX,type=type,se.fit=se.fit,...))
  if (object$basis.x[[vfunc[i]]]$type=="pc") {
   object$beta.l[[vfunc[i]]]$data<-matrix(object$beta.l[[vfunc[i]]]$data,nrow=1)
   b1<-inprod.fdata(fdata.cen(newx[[vfunc[i]]],object$mean.list[[vfunc[i]]])[[1]],object$beta.l[[vfunc[i]]])
   yp<-yp+b1                  
   }
   else{
    xcen<-fdata.cen(newx[[vfunc[i]]],object$mean.list[[vfunc[i]]])[[1]]
    x.fd=Data2fd(argvals=tt,y=t(xcen$data),basisobj=object$basis.x[[vfunc[i]]])
    C=t(x.fd$coefs)
    cnames<-colnames(object$JJ[[vfunc[i]]] )
    b.est<-matrix(object$coefficients[cnames],ncol=1)
    b1<- C%*%object$JJ[[vfunc[i]]]%*%b.est
    yp<-yp+b1  
   }
  } 
  XX2<-as.matrix(cbind(rep(1,len=nn),XX) )
  predictor<-drop(yp)
  if (se.fit || interval != "none") {   
    ip<-rowSums((XX2 %*%object$Vp*XX2))   
    res.var<-object$sr2    
    df<-object$df.residual
    if (interval != "none") {
        tfrac <- qt((1 - level)/2, df)
        hwid <- tfrac * switch(interval, confidence = sqrt(ip),prediction = sqrt(ip + pred.var))    
        predictor <- cbind(predictor, predictor + hwid %o%c(1, -1))
        colnames(predictor) <- c("fit", "lwr", "upr")
    }
  }
  if (se.fit)   {
        se <- sqrt(ip)
        return(list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var)))
        }
else return(predictor)          
 }    
 }
return(drop(yp))  
 }

 

