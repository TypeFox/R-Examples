#################################################################
#################################################################
summary.fdata.comp=function(object,y=NULL,biplot=TRUE,corplot=FALSE,...) {
if (inherits(object, "fdata.comp"))         {
   a1=TRUE
   pr.com<-object
   if (is.null(y)) {
      if (object$call[[1]]=="fdata2pls" | object$call[[1]]=="fdata2ppls")
       y<-object$y
#                    else stop("Argument y no introduced")
                    }
   }
else if (inherits(object, "fregre.fd"))     {
   a1=FALSE
   pr.com<-object$fdata.comp
   y<-object$y
     }
else stop("Error in input data")
#  pr.com=object
  out2=pr.com$x
#  object<-object[["data"]]
#  if (nrow(object) != (length(y))) stop("ERROR IN THE DATA DIMENSIONS")
 l<-object$l
 le=length(l)
 rotation=aperm(pr.com$rotation$data)
 p=pr.com$x
if (object$call[[1]]=="fdata2pls" | object$call[[1]]=="fdata2ppls") {
 var.1<-apply(p, 2, var)
 pr.x2= var.1/sum(var.1)
}
else {
 d<-object$d
 if (is.null(object$rn)) rn<-0
 pr.x2<-(d^2+rn)/sum(d^2+rn)
 names(pr.x2)<-colnames(out2)[l]
 }
# print(pr.x2)
 C<-match.call()
 lenC=length(C)
 cor.y.pc=rep(NA,le)
 xxx=cbind(y,pr.com$x)
 cor.y.pc=round(cor(xxx[,c(1,l+1)]),3)[1,-1]
 types<-colnames(pr.com$x)
 cat("\n     - SUMMARY:  ",object$call[[1]]," object   -\n")
 if (object$call[[1]]=="fdata2pc" | object$call[[1]]=="fdata2ppc") {
   cat("\n-With",le," components are explained ",round(sum(pr.x2[l])*100
 ,2),"%\n of the variability of explicative variables.\n \n-Variability for each component (%):\n")
# print(round(pr.x[l] * 100, 2))
print(round(pr.x2[l] * 100, 2))
}
  if (corplot){
   if (is.null(y)) stop("Argument y no introduced")
   names(cor.y.pc)=paste("cor(y,",types[l],")",sep="")
   cat("\n-Correlations:\n")
   print(cor.y.pc)
   j=1
   dev.new()
     while (j<=lenC) {
       if (names(C)[j]=="ask") {ask=C[[j]];j=lenC +1}
          else {j=j+1;ask=FALSE}
     }
     if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          }
     else   par(mfrow=c(ceiling(le/2),2))
    for (i in 1:le)   plot(pr.com$x[,l[i]],y,main=paste("cor=",
        round(cor.y.pc[i],3)),xlab=colnames(pr.com$x)[l[i]],ylab="y",...)
  }
 if (biplot){
  j=1
  while (j<=lenC) {
        if (names(C)[j]=="ask") {
           ask=C[[j]]
           j=lenC +1
           }
        else { j=j+1; ask=FALSE  }
  }
   dev.new()
   if (ask) {
          par(mfrow=c(1,1))
          dev.interactive()
          oask <- devAskNewPage(TRUE)
          on.exit(devAskNewPage(oask))
          for (i in 1:le) {
          ts.plot(rotation[,l[i]],ylab=c("loadings",l[i],sep=""),
      main=c(paste("components",l[i],"- Expl. Var. ",round(pr.x2[l[i]] * 100, 2),"%",sep="")))
      if (i<le)
      for (j in (i+1):le) {

            if (nrow(out2)<50)   {
                        plot(p[,c(i,l[j])],main="BIPLOT",type="n")
                        text(p[,c(i,l[j])])#,rownames(out2))
                        }
             else                         plot(p[,c(i,l[j])],main="BIPLOT")
            if (nrow(out2)<50)      {
                           plot(p[,c(l[j],i)],main="BIPLOT",type="n")
                           text(p[,c(l[j],i)])#,rownames(out2))
               }
           else  plot(p[,c(l[j],i)],main="BIPLOT")
        } }  }
    else   {
    par(mfrow=c(le,le))
    for (i in 1:le) {
      par(mfg=c(i,i))
      ts.plot(rotation[,l[i]],ylab=c("loadings",l[i],sep=""),
       main=c(paste("Component",l[i],"- Expl. Var. ",round(pr.x2[l[i]] * 100, 2),"%",sep="")))
      if (i<le)
      for (j in (i+1):le) {
            par(mfg=c(i,j))
            if (nrow(out2)<50)     {
                       plot(p[,c(i,l[j])],main="BIPLOT")
                      text(p[,c(i,l[j])])
                      }
            else plot(p[,c(i,l[j])],main="BIPLOT")
            par(mfg=c(j,i))
            if (nrow(out2)<50)            {
               plot(p[,c(l[j],i)],main="BIPLOT",type="n")
               text(p[,c(l[j],i)])#,rownames(out2))
               }
            else plot(p[,c(l[j],i)],main="BIPLOT")
     }  }    }  }
#return(invisible(pr.com))
return(invisible(object))
}
#################################################################
#################################################################



