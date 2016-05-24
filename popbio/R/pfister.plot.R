pfister.plot<-function(A)
{
    n<-length(A)
    # error-checking
    if(class(A)!="list"){stop("A LIST of annual matrices is required")}
    ## matrix inputs should have same dimensions
    if(length(unique(lapply(A, dim)))>1){stop("Matrices have different dimensions")}
    if(n<2){stop("A list of TWO or more annual matrices is required input")}     
    #COLUMN NAMES
   col<-names(A)
   #NUMBER of stages
    x<-dim(A[[1]])[1]
    ## ROW names
   row<-paste("a", 1:x, rep(1:x,each=x), sep="")
   ## annual matrix elements
   vr<-data.frame(matrix(unlist(A), ncol=n, dimnames=list(row, col)))
   ## MEAN, var, and cv
   vr$mean<-apply(vr, 1, mean)
   vr$var<-apply(vr, 1, var)
   vr$cv<-vr$var^.5/vr$mean*100
   ## mean matrix
   meanA<- matrix(vr$mean, nrow=x)
   ## Sensitivities and elasticities of mean matrix
   eigA<-eigen.analysis(meanA)
   vr$sens<-as.vector(eigA$sensitivities)
   vr$elas<-as.vector(eigA$elasticities)
   #NON-ZERO elements
   vr1<-subset(vr, mean>0 & var>0)
   #### log-log PLOTS     
   op<-par(mfrow=c(1,2))
   plot(vr1$var, vr1$sens, xlab="Variance", ylab="Sensitivity", log="xy", pch=16, col="blue")
   plot(vr1$cv, vr1$elas, xlab="CV", ylab="Elasticity", log="xy", pch=16, col="blue")

   par(op)
   # output plot values 
   vr1[,(n+1):(n+5)]
}

