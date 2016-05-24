`cca` <-
function(x,y,xlab=colnames(x),ylab=colnames(y),xcenter=TRUE,ycenter=TRUE,xscale=FALSE,yscale=FALSE,standardize.scores=TRUE,use="complete.obs",na.rm=TRUE){
   #Perform the preliminaries
   if(is.data.frame(x))                  #Some routines require matrices...
     x<-as.matrix(x)
   if(is.data.frame(y))
     y<-as.matrix(y)
   #Center the data matrices, if needed
   if(xcenter||xscale)
      x<-scale(x,scale=xscale,center=xcenter)   
   if(ycenter||yscale)
      y<-scale(y,scale=yscale,center=ycenter) 
   if(is.null(dim(x)))                   #Make sure we have matrices
      x<-matrix(x,ncol=1)
   if(is.null(dim(y)))
      y<-matrix(y,ncol=1)
   #Find out how large these data matrices are
   nx<-dim(x)[2]      
   ny<-dim(y)[2]
   ncv<-min(nx,ny)
   cvlab<-paste("CV",1:ncv)
   o<-list()
   #Get covariance matrices
   cxx<-cov(x,use=use)
   cyy<-cov(y,use=use)
   cxy<-cov(x,y,use=use)
   cyx<-t(cxy)
   #Find the projections
   ey<-eigen(qr.solve(cyy,cyx)%*%qr.solve(cxx,cxy))
   ex<-list(values=ey$values,vectors=qr.solve(cxx,cxy)%*%(ey$vec))
   o$corr<-(ex$val[1:ncv])^0.5
   names(o$corr)<-cvlab
   o$corrsq<-o$corr^2                                 #Get the variance accounted for by each canonical variate
   names(o$corrsq)<-cvlab
   o$xcoef<-ex$vec[,1:ncv,drop=FALSE]
   rownames(o$xcoef)<-xlab
   colnames(o$xcoef)<-cvlab
   o$ycoef<-ey$vec[,1:ncv,drop=FALSE]
   rownames(o$ycoef)<-ylab
   colnames(o$ycoef)<-cvlab
   #Find the canonical variates (using the coefficients) and compute structural correlation information
   o$canvarx<-x%*%o$xcoef    #Construct the canonical variates
   rownames(o$canvarx)<-rownames(x)
   colnames(o$canvarx)<-cvlab
   o$canvary<-y%*%o$ycoef
   rownames(o$canvary)<-rownames(y)
   colnames(o$canvary)<-cvlab
   if(standardize.scores){         #If needed, standardize the scores/coefs
     sdx<-apply(o$canvarx,2,sd)
     sdy<-apply(o$canvary,2,sd)
     o$canvarx<-sweep(o$canvarx,2,sdx,"/")
     o$canvary<-sweep(o$canvary,2,sdy,"/")
     o$xcoef<-sweep(o$xcoef,2,sdx,"/")
     o$ycoef<-sweep(o$ycoef,2,sdy,"/")
   }
   o$xstructcorr<-cor(x,o$canvarx,use=use)   #Find structural correlations
   rownames(o$xstructcorr)<-xlab
   colnames(o$xstructcorr)<-cvlab
   o$ystructcorr<-cor(y,o$canvary,use=use)
   rownames(o$ystructcorr)<-ylab
   colnames(o$ystructcorr)<-cvlab
   o$xstructcorrsq<-o$xstructcorr^2         #Find variance explained by structural correlations
   rownames(o$xstructcorrsq)<-xlab
   colnames(o$xstructcorrsq)<-cvlab
   o$ystructcorrsq<-o$ystructcorr^2 
   rownames(o$ystructcorrsq)<-ylab
   colnames(o$ystructcorrsq)<-cvlab
   o$xcrosscorr<-cor(x,o$canvary,use=use)   #Find cross-correlations
   rownames(o$xcrosscorr)<-xlab
   colnames(o$xcrosscorr)<-cvlab
   o$ycrosscorr<-cor(y,o$canvarx,use=use)
   rownames(o$ycrosscorr)<-ylab
   colnames(o$ycrosscorr)<-cvlab
   o$xcrosscorrsq<-o$xcrosscorr^2   #Find variance exp. by cross-correlations
   rownames(o$xcrosscorrsq)<-xlab
   colnames(o$xcrosscorrsq)<-cvlab
   o$ycrosscorrsq<-o$ycrosscorr^2
   rownames(o$ycrosscorrsq)<-ylab
   colnames(o$ycrosscorrsq)<-cvlab
   o$xcancom<-apply(o$xstructcorrsq,1,sum,na.rm=na.rm)      #Find the canonical communalities (total var explained)
   names(o$xcancom)<-xlab
   o$ycancom<-apply(o$ystructcorrsq,1,sum,na.rm=na.rm)
   names(o$ycancom)<-ylab
   o$xcanvad<-apply(o$xstructcorrsq,2,mean,na.rm=na.rm)      #Find the canonical variate adequacies
   names(o$xcanvad)<-cvlab
   o$ycanvad<-apply(o$ystructcorrsq,2,mean,na.rm=na.rm)
   names(o$ycanvad)<-cvlab
   o$xvrd<-o$xcanvad*o$corrsq                             #Find the redundancy indices (Rd) for X|Y and Y|X
   names(o$xvrd)<-cvlab
   o$yvrd<-o$ycanvad*o$corrsq
   names(o$yvrd)<-cvlab
   o$xrd<-sum(o$xvrd,na.rm=na.rm)
   o$yrd<-sum(o$yvrd,na.rm=na.rm)
   bartvbase<--(NROW(x)-1-(nx+ny+1)/2)              #Bartlett's chisq stats
   o$chisq<-bartvbase*(sum(log(1-o$corr^2))-c(0,cumsum(log(1-o$corr^2))[-ncv]))
   o$df<-(nx+1-(1:ncv))*(ny+1-(1:ncv))
   names(o$chisq)<-cvlab
   names(o$df)<-cvlab
   o$xlab<-xlab                                     #Save labels, just in case
   o$ylab<-ylab
   #print(eigen(solve(cxx)%*%cxy%*%solve(cyy)%*%cyx))   
   #Return the results
   class(o)<-"cca"
   o
}

