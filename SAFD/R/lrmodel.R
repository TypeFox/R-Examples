lrmodel <-
function(XX, YY,theta=1/3){
 #given data for the regression XX and YY of trapezoidal fuzzy sets
 #nobs...sample size
 #in case the data is not trapezoidal it is automatically transformed to
 #trapezoidal one using the translator function
 
 kx<-length(XX)
 ky<-length(YY)

 if(kx!=ky){
   print("lists must have same length (i.e. input and output must have same sample size")
   }
 if(kx==ky){
  nobs<-kx
  for (i in 1:nobs){
   XX[[i]]<-translator(XX[[i]],2)
   YY[[i]]<-translator(YY[[i]],2)
  }

 #bind XX and YY in list to simplify check for compatibility
 ZZ<-vector("list",length=(2*nobs))
 ZZ[1:nobs]<-XX[1:nobs]
 ZZ[(nobs+1):(2*nobs)]<-YY[1:nobs]
 temp_sum<-Msum(ZZ)
 if(nrow(temp_sum)>1){
  #calculate feasible set Feas =[-amax0,bmax0]
  XXs <-rep(0,nobs)
  YYs <-rep(0,nobs)
  XXl <-rep(0,nobs)
  XXr <-rep(0,nobs)
  YYl <-rep(0,nobs)
  YYr <-rep(0,nobs)
  M<-data.frame(XXs=XXs,YYs=YYs,XXl=XXl,YYl=YYl,XXr=XXr,YYr=YYr)

  for(i in 1:nobs){
	 M$XXs[i] <-0.5*(XX[[i]]$x[3]-XX[[i]]$x[2])
	 M$YYs[i] <-0.5*(YY[[i]]$x[3]-YY[[i]]$x[2])
	 M$XXl[i] <-(XX[[i]]$x[2]-XX[[i]]$x[1])
	 M$YYl[i] <-(YY[[i]]$x[2]-YY[[i]]$x[1])
	 M$XXr[i] <-(XX[[i]]$x[4]-XX[[i]]$x[3])
	 M$YYr[i] <-(YY[[i]]$x[4]-YY[[i]]$x[3])
	}
	
	 M1<-subset(M,M$XXs>0&M$XXl>0&M$XXr>0) 
   if(nrow(M1)==0){
    Feas<-c(NA,NA)
    }
  
   if(nrow(M1)>0){
    M1$valneg1 <- M1$YYs/M1$XXs
    M1$valpos2 <- M1$YYl/M1$XXl
	  M1$valneg2 <- M1$YYr/M1$XXl
	  M1$valpos3 <- M1$YYr/M1$XXr
	  M1$valneg3 <- M1$YYl/M1$XXr
   
    bmax0 <-min(M1$valneg1,M1$valpos2,M1$valpos3)
    amax0 <-min(M1$valneg1,M1$valneg2,M1$valneg3)
	  Feas<-c(-amax0,bmax0)
	  }

  varX  <- Bvar(XX,theta)
  covXY  <- Bcov(XX,YY,theta)
  mXX<-list()
  for(i in 1:nobs){
   mXX[[i]]<-sc_mult(XX[[i]],-1)
  }
  covmXY <- Bcov(mXX, YY,theta)
  #cat(varX, covXY, covmXY, "\n")
  #calculate beta and gamma
   beta<-0
   if(covmXY >0){
    if(is.na(Feas[1])==TRUE){
    beta<-1
    }
    if(is.na(Feas[1])==FALSE){
 	 beta <- min(1,Feas[1] * varX/covmXY)
 	 }
 	}
   gamma<-0
   if (covXY >0){
    if(is.na(Feas[1])==TRUE){
     gamma<-1
    }
    if(is.na(Feas[1])==FALSE){
     gamma <- min(1,Feas[2] * varX/covXY)
 	 }
   }
  #cat("beta=", beta, "gamma=", gamma, "\n")
  #print(c(beta,gamma))
 #estimate a and B
  if (gamma==0|beta==0){
  	a<-gamma*covXY/varX - beta*covmXY/varX
    Best<- hukuhara(sc_mult(Mmean(XX),a),Mmean(YY),0)
    Result<-list(a=a,B=Best)
    }
  if (gamma!=0&beta!=0){
  	cond<-covmXY/covXY - (2*gamma-gamma^2)/(2*beta-beta^2)
    if(cond>0){
     a<-beta*covmXY/varX
     Best<- hukuhara(sc_mult(Mmean(XX),a),Mmean(YY),0)
     Result<-list(a=a,B=Best)
     }
    if(cond<0){
     a<-gamma*covXY/varX
     Best<- hukuhara(sc_mult(Mmean(XX),a),Mmean(YY),0)
     Result<-list(a=a,B=Best)
    }
    if(cond==0){
     a<-c(beta*covmXY/varX,gamma*covXY/varX)
     B<-vector("list",length=2)
     B[[1]]<- hukuhara(sc_mult(Mmean(XX),a[1]),Mmean(YY),0)
     B[[2]]<- hukuhara(sc_mult(Mmean(XX),a[2]),Mmean(YY),0)
     Result<-list(a=a,B=B)
     }
    }
    invisible(Result)
   }
 }
}
