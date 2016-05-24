slider.bootstrap.lm.plot<-function(x,y=NULL,...)
{ 
 x.name<-deparse(substitute(x))
 y.name<-deparse(substitute(y))
 if(length(x)<2) return("Error: x is of length 0 or 1")
 if(!is.null(y)){ 
   if(length(y)<2) return("Error: y must be a vector")
   if(length(x)!=length(y)) 
     return("Error: x and y must have the same length")
   x<-cbind(x,y) 
 }
 if(!is.matrix(x) && !is.data.frame(x)){ 
   x<-cbind(seq(x),x) 
   y.name<-x.name; x.name<-"index"
 }
 if(is.null(y.name)){x.name<-colnames(x)[1]; y.name<-colnames(x)[2]}
 y<-x[,2]; x<-x[,1]

 args<-list(...)
 n<-length(x)
 ind<-order(x); x.orig<-x<-x[ind]; y.orig<-y<-y[ind]
 xx<-seq(min(x),max(x),length=100)         
 # plot(x,y,...)       
 do.call("plot",c(alist(x,y,bty="n"),args))
 abline(lm(y~x),lwd=5)
 refresh<-function(...){
   # plot(x,y,...); 
   do.call("plot",c(alist(x,y,bty="n"),args))
   abline(coefyx<-lm(y~x)$coef, lwd=3)
   polytype<-slider(no=1)
   form<-paste(paste(sep="","I(x^",1:polytype,")"),collapse="+")
   form<-as.formula(paste("y ~",form)); coef<-lm(form)$coef
   yy<-outer(xx,0:polytype,"^")%*%coef; lines(xx,yy,lwd=2)
   B<-slider(no=2); zz<-slider(no=3); set.seed(zz)
   result<-matrix(0,1+polytype,B)
   for(i in 1:B){
     index<-sample(1:n,n,replace=TRUE)
     x<-x.orig[index]; y<-y.orig[index]
     coef<-lm(form)$coef
     yy<-outer(xx,0:polytype,"^")%*%coef
     lines(xx,yy,lwd=2,col=i,lty=2)
     result[,i]<-coef
   }
   abline(coefyx, lwd=5)
   result<-t(result); 
   colnames(result)<-c("intercept",paste(sep="","beta: x^",1:polytype))
   print(summary(result))
 }
 slider(refresh,c("polynomial degree","number repetitions","random seed"),
        c(1,1,1),c(5,50,100),c(1,1,1),c(1,10,1))
 refresh()
 "ok"
}

