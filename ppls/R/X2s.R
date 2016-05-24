`X2s` <-
function(X,Xtest=NULL,deg=3,nknot=NULL,reduce.knots=FALSE){


    ####################################################
    # internal function
    ###################################
    
    compute.knots<-function(x,nknot,deg){
        xl<-min(x);
        xr<-max(x);
	a<-xr-xl
	if (a==0){
		a<-1
	}
        xmin<-xl-a/100
        xmax<-xr+a/100
        dx<-(xmax-xmin)/(nknot-1)
        knots<-seq(xmin-deg*dx,xmax+deg*dx,by=dx)
        return(knots)
    }
    #############################
    
  if (is.null(Xtest)){

    Xtest<-X

  }

  if (is.vector(X)){

    X<-matrix(X,ncol=1)

  }


    if (is.vector(Xtest)){

        Xtest<-matrix(Xtest,ncol=1)

      }


    if (ncol(X)!=ncol(Xtest)){

      "error: X and Xtest must have the same number of columns"

    }

    p<-ncol(X)
  
    if (is.null(nknot)==TRUE){
        nknot<-rep(20,p)
    }

  
  
  nknots<- nknot -2*deg # transform the number of knots such that
                        # compute knots returns the correct amount
                        # of knots    
  
  for (j in 1:p){
    nknots[j]=max(nknots[j],deg+1) # there must be at least as many knots
                                    # as the order 
  }
  
  
    sizeZ<-vector(length=p)
  
  
  
    Z<-c()

  Ztest<-c()
  for (j in 1:p){
  if (reduce.knots==TRUE){
  too.many.knots=TRUE
  while ((too.many.knots==TRUE) & (nknots[j]>deg+1)){
  too.many.knots=FALSE
    knotsj<-compute.knots(X[,j],nknots[j],deg)

    #Zj<-spline.des(knots=knotsj,x=X[,j],ord=deg+1,outer.ok=TRUE)$design
    Zj<-splineDesign(knots=knotsj,x=X[,j],ord=deg+1,outer.ok=TRUE)

    #Zjtest<-spline.des(knots=knotsj,x=Xtest[,j],ord=deg+1,outer.ok=TRUE)$design
    Zjtest<-splineDesign(knots=knotsj,x=Xtest[,j],ord=deg+1,outer.ok=TRUE)
    if (min(apply(Zj^2,2,sum))==0){
  
    too.many.knots=TRUE
    nknots[j]=nknots[j]-1
    }
}
    if (nknots[j]==deg+1){
	nknots[j]=nknots[j]+1
        knotsj<-compute.knots(X[,j],nknots[j],deg)
	Zj<-splineDesign(knots=knotsj,x=X[,j],ord=deg+1,outer.ok=TRUE)

    #Zjtest<-spline.des(knots=knotsj,x=Xtest[,j],ord=deg+1,outer.ok=TRUE)$design
    Zjtest<-splineDesign(knots=knotsj,x=Xtest[,j],ord=deg+1,outer.ok=TRUE)
	}
    Z=cbind(Z,Zj)
    Ztest=cbind(Ztest,Zjtest)
    sizeZ[j]=ncol(Zj)
    }
  else{
  knotsj<-compute.knots(X[,j],nknots[j],deg)

    Zj<-spline.des(knots=knotsj,x=X[,j],ord=deg+1,outer.ok=TRUE)$design
    

    Zjtest<-spline.des(knots=knotsj,x=Xtest[,j],ord=deg+1,outer.ok=TRUE)$design
    Z=cbind(Z,Zj)
    Ztest=cbind(Ztest,Zjtest)
    sizeZ[j]=ncol(Zj)
  
  }
  }
  
  return(list(Z=Z,Ztest=Ztest,sizeZ=sizeZ))

}
