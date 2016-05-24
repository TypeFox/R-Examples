predict.fregre.glm<-function(object,newx=NULL,type="response",...){
 if (is.null(object)) stop("No fregre.glm object entered")
 if (is.null(newx)) {
    yp=predict.glm(object,type=type,...)
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
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
if   (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
##########
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
# vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
# vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 beta.l=list()
 kterms=1
if (length(vnf)>0) {
 first=FALSE
 XX=data.frame(data[["df"]][,c(vnf)])
 names(XX)=vnf
 for ( i in 1:length(vnf)){
# print(paste("no functional variable",vnf[i]))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (attr(tf,"intercept")==0) {
     print("No intecept")
     pf<- paste(pf,-1,sep="")
     }
}
else  first=TRUE
if (length(vfunc)>0)  {
#   yp2<-a1 <- object$coefficients[1] * rep(1, len = nrow(data[[vfunc[1]]]))
   for (i in 1:length(vfunc)) {
   if(class(data[[vfunc[i]]])[1]=="fdata")  {
      fdataobj<-data[[vfunc[i]]]
      x.fd<-fdataobj[["data"]]
      if (nrow(x.fd)==1) rwn<-NULL
      else rwn<-rownames(x.fd)
      tt<-fdataobj[["argvals"]]
      rtt<-fdataobj[["rangeval"]]
      if (!object$basis.x[[vfunc[i]]]$type=="pc"&!object$basis.x[[vfunc[i]]]$type=="pls") {
 	  	x.fd = Data2fd(argvals = tt, y = t(fdata.cen(fdataobj,object$mean[[vfunc[i]]])[[1]]$data),
                      basisobj = basis.x[[vfunc[i]]],fdnames=rwn)
	    	r=x.fd[[2]][[3]]
        J<-object$JJ[[vfunc[i]]]
        Z = t(x.fd$coefs) %*% J
#        colnames(Z) = paste(vfunc[i], ".",colnames(J), sep = "")
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
if (first) return(rep(object$coefficient,length=nrow(newx[[1]])) )        
if (!is.data.frame(XX)) XX=data.frame(XX)         
 yp=predict.glm(object=object,newdata=XX,type=type,x=TRUE,y=TRUE,...)
return(yp)
}
}


