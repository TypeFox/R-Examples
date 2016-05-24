bootstrap_lm_cov_latent <- function(X1,X2,param="multilogit",Psi,Be,Ga,B=100,fort=TRUE){
    
  mPsi = 0
  mBe = 0
  m2Psi = 0
  m2Be = 0
  if(param=="multilogit"){
  	mGa = 0
  	m2Ga = 0
  }else if(param=="difflogit"){
  	mGa = vector("list",2)
    m2Ga = vector("list",2)
    mGa[[1]] = matrix(0,dim(Ga[[1]]))
    mGa[[2]] = matrix(0,dim(Ga[[2]]))
    m2Ga[[1]] = matrix(0,dim(Ga[[1]]))
    m2Ga[[2]] = matrix(0,dim(Ga[[2]]))
  }

  dPsi = dim(Psi)
  if(length(dPsi)==1) k = 1
  else k = dPsi[2]
  for (b in 1:B) {
    cat("boostrap sample n. ",b,"\n")
    out = draw_lm_cov_latent(X1,X2,param,Psi,Be,Ga,fort=fort)
    Yb = out$Y
    out = est_lm_cov_latent(Yb,X1,X2,param=param,k=k,fort=fort)
    mPsi = mPsi + out$Psi/B
    mBe = mBe + out$Be/B
    if(param=="multilogit"){
   	 	mGa = mGa + out$Ga/B
    	m2Ga = m2Ga + out$Ga^2/B   	
    }else if(param=="difflogit"){
    	mGa[[1]] = mGa[[1]]+out$Ga[[1]]/B
    	mGa[[2]] = mGa[[2]]+out$Ga[[2]]/B
    	m2Ga[[1]] = m2Ga[[1]] + out$Ga[[1]]^2/B   	
      	m2Ga[[2]] = m2Ga[[2]] + out$Ga[[2]]^2/B   	
    }
    m2Psi = m2Psi + out$Psi^2/B
    m2Be = m2Be + out$Be^2/B
 
  }
  sePsi = sqrt(m2Psi - mPsi^2)
  seBe = sqrt(m2Be - mBe^2)
  if(param=="multilogit"){ 
  	seGa = sqrt(m2Ga - mGa^2)
  }else if(param=="difflogit"){
  	seGa = vector("list",2)
  	seGa[[1]] = sqrt(m2Ga[[1]] - mGa[[1]]^2)
  	seGa[[2]] = sqrt(m2Ga[[2]] - mGa[[2]]^2)	
  }
  out = list(mPsi = mPsi, mBe = mBe, mGa = mGa, sePsi = sePsi, 
             seBe = seBe, seGa = seGa)
  
}