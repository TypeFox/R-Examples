space.discretized<-function(xinit,yuima, env){


##:: initialize state variable
	sdeModel<-yuima@model
	
	modelstate <- sdeModel@solve.variable
	modeltime <- sdeModel@time.variable
	V0 <- sdeModel@drift
	V <- sdeModel@diffusion
	r.size <- sdeModel@noise.number
	d.size <- sdeModel@equation.number
	Terminal <- yuima@sampling@Terminal[1]
	n <- yuima@sampling@n[1]
#	dX <- xinit
  if(length(unique(as.character(xinit)))==1 &&
       is.numeric(tryCatch(eval(xinit[1],env),error=function(...) FALSE))){
    dX_dummy<-xinit[1]
    dummy.val<-eval(dX_dummy, env)
    if(length(dummy.val)==1){dummy.val<-rep(dummy.val,length(xinit))}
    for(i in 1:length(modelstate)){
      assign(modelstate[i],dummy.val[i] ,env)
    }
    dX<-vector(mode="numeric",length(dX_dummy))
  
    for(i in 1:length(xinit)){
      dX[i] <- dummy.val[i]
    }
  }else{
	  dX_dummy <- xinit
	  if(length(modelstate)==length(dX_dummy)){
	    for(i in 1:length(modelstate)) {
	      if(is.numeric(tryCatch(eval(dX_dummy[i],env),error=function(...) FALSE))){
	        assign(modelstate[i], eval(dX_dummy[i], env),env)
	      }else{
	        assign(modelstate[i], 0, env)
	      }
	    }
	  }else{ 
	    yuima.warn("the number of model states do not match the number of initial conditions")
	    return(NULL)
	  }
	# 20/11 we need a initial variable for X_0
	  dX<-vector(mode="numeric",length(dX_dummy))

	  for(i in 1:length(dX_dummy)){
	    dX[i] <- eval(dX_dummy[i], env)
	  }
  }
	
	
##:: set time step	
	delta <- Terminal/n 

	
	
	
##:: using Space-discretized Euler-Maruyama method
	


##:: function for approximation of function G
gfunc <- function(x){
  c0 <- 10
  c1 <- 10
  ret <- rep(0, length(x))
  idx <- which(x < 1/c0)
  ret[idx] <- 1
  
  idx <- which(1/c0 <= x)
  ret[idx] <- 1-pnorm(x[idx])
  for(i in 1:length(idx)){
    n <- 1:floor(c1/x[idx[i]])
    ret[idx[i]] <- 4 * (ret[idx[i]] - sum( pnorm((4*n+1)*x[idx[i]]) - pnorm((4*n-1)*x[idx[i]]) ))
  }
  
  idx <- which(1 < ret)
  ret[idx] <- 1
  return(ret)
}


   
    dxx <- 0.0001
    xx <- seq(0, 1.7, dxx)
    
    ##:: approximate function G(gg)
    gg <- gfunc(xx)
    appfunc <- suppressWarnings(approxfun(gg, xx))
    
    ##:: calculate inverse of G
    unif.a <- runif(n*2)
    inv.a <- pmin(qnorm(1 - unif.a/4), appfunc(unif.a), na.rm=TRUE)
    
    ##:: make random time steps
    ep <- sqrt(delta)
    dTW <- (ep/inv.a)^2
    time_idx <- cumsum(dTW) ##:: time index should be attached            
    div_sd <- min(which(time_idx > Terminal)) ##:: cut by time=1
    time_idx <- time_idx[1:div_sd]
    
    ##:: add diffusion term
    dTW <- rbind(dTW[1:div_sd],
                 t(matrix( (rbinom(div_sd*r.size, 1, 0.5)*2-1) * ep,
                          nrow=div_sd,
                          ncol=r.size)
                   )
                 )
    
    X_mat <- matrix(0, d.size, div_sd+1)              
    X_mat[,1] <- dX
    
    ##:: function to calculate coefficients of dTW
    p.b <- function(t, X=numeric(d.size)){
      ##:: assign names of variables
      for(i in 1:length(modelstate)){
        assign(modelstate[i], X[i], env)
      }
      assign(modeltime, t, env)
      tmp <- matrix(0, d.size, r.size+1)
      for(i in 1:d.size){
        tmp[i,1] <- eval(V0[i],env)
        for(j in 1:r.size){
          tmp[i,j+1] <- eval(V[[i]][j], env)
        }
      }
      return(tmp)
    }
    ##:: calcurate difference equation
    for(i in 1:div_sd){
      dX <- dX + p.b(t=time_idx[i], X=dX) %*% dTW[,i]
      X_mat[,i+1] <- dX
    }
    ##tsX <- ts(data=t(X_mat), deltat=delta , start=0)
    ##:: output zoo data
    zooX <- zoo(x=t(X_mat), order.by=c(0, time_idx))
    yuimaData <- setData(original.data=zooX)
	return(yuimaData)

}
 