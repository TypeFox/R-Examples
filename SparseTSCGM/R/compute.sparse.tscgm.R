

compute.sparse.tscgm <- function(X=X, Y=Y, lam1=lam1, lam2=lam2,
       optimality = c("NULL", "bic","bic_ext","bic_mod","aic","gic"),
       setting=setting)
{

  T <- dim(Y)[1]
  p <- dim(X)[2]
  n <- dim(Y)[3]
  q <- dim(Y)[2]
  xtyi <- array(NA, c(p,q,n))
  xtxi <- array(NA, c(p,p,n))
  ytyi <- array(NA, c(q,q,n))

  for(i in 1:n){
      XX <- X[,,i]
      YY <- Y[,,i]
      xtyi[,,i]=crossprod(XX,YY)
      xtxi[,,i]=crossprod(XX)
      ytyi[,,i]=crossprod(YY)

    }
  xty=apply(xtyi, c(1,2), sum)
  xtx=apply(xtxi, c(1,2), sum)
  yty=apply(ytyi, c(1,2), sum)
  xtxt=apply(xtxi, c(1,2), sum)/(n*T)

  nlam=(n*T)*lam2
  #starting value for matrix autoregressive coefficients
  old.B = qr.solve(xtx + nlam*diag(p), xty)
  if(!is.numeric(old.B))  old.B <- matrix(0,p,q)

  k=0
  mab =  sum(sum(abs(old.B)))
  while(1)
  {
    k=k+1
        ##Initial estimates for the precision matrix using Glasso
        samp.cov = ( yty - t(xty) %*% old.B - t(old.B) %*% xty + t(old.B) %*% xtx %*% old.B)/(n*T)
        #samp.cor=cov2cor(samp.cov)
        g.out=glasso(s=samp.cov, rho=lam1, thr=1.0e-4, maxit=1e4, penalize.diagonal=FALSE,
        approx=FALSE ) #, start="warm", w.init=old.om0i, wi.init=old.om0)
        old1.om=g.out$wi
	      old1.om.i=g.out$w
         if(!is.numeric(old1.om))  old1.om <- diag(q)
        # Calculating SCAD penalty
        wt <- matrix(NA, nrow=q,ncol=q)
        a <- 3.7                                       
        for(i in 1:q){
          for(j in 1:q){
           if(abs(old1.om[i,j]) <= lam1 ) wt[i,j] <- lam1
             else { if((lam1 <= abs(old1.om[i,j])) & (abs(old1.om[i,j]) <  a*lam1 )) {
             wt[i,j] <- ((a*lam1-abs(old1.om[i,j]))/((a-1))) }
          	   else wt[i,j] <-  0 }
        }}
        diag(wt) <- lam1
        ##Estimating the precision matrix using Glasso based on SCAD penalty
        g.out1=glasso(s=samp.cov, rho=wt, thr=1.0e-4, maxit=1e4, penalize.diagonal=FALSE,
        approx=FALSE, start="warm", w.init=old1.om.i, wi.init=old1.om)
	    old.om=g.out1$wi
	    old.om.i=g.out1$w
       if(!is.numeric(old.om) )  old.om <- diag(q)
      xtyom=(xty%*%old.om)
      warmstart=1
 	    if(k == 1) warmstart=0

      ##Estimating the auotregressive coeficients based on SCAD penalty
        wt1 <- matrix(NA, nrow=p,ncol=q)
        a <- 3.7
        for(i in 1:p){
               if(xtx[i,i] == 0) xtx[i,i] <- 1 
             
          for(j in 1:q){
                if(old.om[j,j] == 0) old.om[j,j] <- 1 
               
              if(abs(old.B[i,j]) <= lam2 ) wt1[i,j] <- lam2
              else {if((lam2 < abs(old.B[i,j])) & (abs(old.B[i,j]) <  a*lam2 )) {
                  wt1[i,j] <- ((a*lam2-abs(old.B[i,j]))/((a-1))) }
               else wt1[i,j] <-   0  }

          }}
        
        rho2 <- wt1*(n*T)
 	     B=rblasso(s=xtx, m=xtyom, om=old.om, nlam=rho2, tol=1e-5, sbols=mab, maxit=setting$maxit.in, warm=warmstart, B0=old.B)
  	   bdist = sum(sum(abs(B-old.B)))
  	   if(!is.numeric(B)) B <- matrix(0,p,q) 
   	   old.B=B
          if( (bdist < setting$tol.out*mab) | (k > setting$maxit.out))
    	      break
     cat("Outer iterations: ", k, "\n")
     cat("lambda1=",lam1, "\n")
     cat("lambda2=",lam2, "\n")
   }
    
  if(setting$silent ==FALSE) cat("Total outer iterations for tscgm : ", k, "\n")
  
  return(list(gamma=old.B, theta=old.om))
}

