spl<-function(x, k=10, knots=NULL, type="LRTP"){

  if(is.null(knots)){
    knots<-quantile(as.numeric(x), 1:k/(k+1), na.rm=T)
  }

  if(type=="LRTP"){  # low-rank thin plate slpine
    Z<-outer(as.numeric(x), knots, function(x,z){abs(x-z)^3})
    B<-outer(knots, knots, function(x,z){abs(x-z)^3})
  }else{
    stop("sorry, only low-rank thin plate splines are implemented")
  }
  if(type=="cubic"){ # non-cyclic cubic spline
  }
  if(type=="cyc-cubic"){ # cyclic cubic spline
  }
  if(type=="res-cubic"){ # restricted cubic spline
  }

  Bsvd<-svd(B)
  nonzeros<-which(Bsvd$d>sqrt(.Machine$double.eps))

  Z<-Z%*%Bsvd$v[,nonzeros, drop=FALSE]%*%diag(sqrt(1/Bsvd$d[nonzeros]))%*%t(Bsvd$u[,nonzeros, drop=FALSE])

  return(Z)
} 


