# This function optimizes the parameters of 
# a H component mixture of k-variate student-t densities
# using IS weighted EM algorithm, given N draws from 
# the initial H component mixture
#
# inputs:
#   theta        : [Nxk matrix] of N parameter draws, each with k parameters
#   mit          : [list] initial mixture density with H (>=1) components.
#                Following components must be provided:
#      p         : [H vector] probability of each component
#      mu        : [Hxd matrix] modes (in row) of each component
#      Sigma     : [Hx(d*d)  matrix] scales (in row, vectorized) of each component
#      df        : [H vector] degree of freedom for each component
#    w           : [N vector] of IS weights
#    control     : [list] of optimization parameters. Following arguments can be provided:
#      maxit.EM  : [double>0] max. number of iterations for the EM algorithm, default: 1000
#      tol.EM    : [double>0] tolerance for EM steps' convergence, default: 0.001
#      trace.EM  : [logical] tracing information on the fitting procedure of IS_EM(logical). 
#                (Default: FALSE; i.e., no tracing information)
#      optim.df  : [logical] should df be optimized? default: TRUE, df parameters are optimized in the program
#                (rest of the inputs relevant only if df are optimized - df=NULL)
#      inter.df  : [2-vector] range of search values for df optimization, default: c(0,50)
#      tol.df    : [double>0] tolerance for degree of freedom optimization step (default 0.0001220703)
#      maxit.df  : [integer>0] maximum number of iterations for degree of freedom optimization (default 1000)
#      trace.df  : [logical] tracing information on the fitting procedure of degrees of freedom (logical). 
#               (Default: FALSE, i.e., no tracing information.
# outputs: [list] with the following components:
#   mit        : [list] containing updated mixture density (has the same structure as input 'mit')
#   summary.EM : [list] containing optimization messages
#
# note: If one t component has non-pds scale or too small mixture probability, the component is removed.
# author : Nalan Basturk
# date   : 20120912

'fn.optimt' <- function(theta,mit,w,control=list()){
   # set control parameters
   con <- list(maxit.EM = 1000, tol.EM = 0.001,trace.EM=FALSE,
   optim.df=TRUE,inter.df=c(0.01,50),tol.df=.Machine$double.eps^0.25,
   maxit.df=1000,tol.pr=0,trace.df=TRUE)
   con[names(control)] <- control
      
   N <- nrow(theta)   # number of parameter draws
   k <- ncol(theta)   # number of parameters to optimize
      
   # ISEM step until convergence
    conv       <- 0;  # convergence indicator
    iter       <- 0;  # number of current iteration
    loglik.old <- -Inf # previous log-likelihood from ISEM
    mit.old    <- mit; # old mixture component     
    H          <- length(mit.old$p)
    while(conv==0 & iter < con$maxit.EM){
     shrinkmit = 1 # indicator for 'shrinking' the Mit density  
     # apply/repeat optimization if mixture component is 'shrank'
     while(shrinkmit){
       iter = iter + 1
       # IS-EM iteration
       out.ISEM  <- fn.ISEM(theta=theta,mit=mit.old,w,control=con)
       mit.new   =  out.ISEM$mit
       # shrink mit density only if H>1
       if(H > 1){
         tmp       <- fn.shrink.mit(mit.new,tol.pr=con$tol.pr) 
         shrinkmit =  tmp$crash
         mit.old   <- tmp$mit
       }else{
         shrinkmit = 0
       } 
       if(shrinkmit)
         iter = 0
     }  
     H <- length(mit.new$p) # number of mixture components

     # new (IS-weighted) loglikelihood of the mixture    
     tmp <- dmvgt(theta,mit=mit.new,log=TRUE)
     loglik.new <- sum(w*tmp)

     # convergence critera
     if(iter >1 & abs((loglik.new-loglik.old)/loglik.old) <= con$tol.EM)
       conv = 1
     loglik.old = loglik.new
   }
   mit=mit.new
   if(iter==con$maxit.EM)
     conv= 2 # maximum number of iterations reached

   summary.EM <- data.frame(iter, conv)
   names(summary.EM) <- c("EM iter","EM conv")
   
   return(list(mit=mit,summary.EM=summary.EM))
}
fn.ISEM<-function(theta,mit,w,control){
   N <- nrow(theta) # number of parameter draws
   k <- ncol(theta) # number of parameters to optimize
   H <- length(mit$p) # number of mixture components
   # auxiliary matrices
   dfmat  <- matrix(mit$df,N,H,byrow=TRUE) # NxH matrix with df in rows
   w.mat  <- matrix(w,N,H) # NxH matrix with IS weights in columns
   z      <- matrix(nrow=N,ncol=H)  # z <- tilde z <- E(z|theta) expected student-t indicator
   z.wg   <- matrix(nrow=N,ncol=H)  # z.wg<- tilde (z/w) <- E(z/w|theta) weighted student-t indicator
   wg.ln  <- matrix(nrow=N,ncol=H)  # wg.ln <- tilde (log(w)) <- E(log(w|theta)) expected student-t indicator
   delta  <- matrix(nrow=N,ncol=H)  # delta <- tilde (1/w) <- E(1/w) expected (inverse) IG draw
   rho    <- matrix(nrow=N,ncol=H)  # (theta-mu)'Sigma(theta-mu)

   # 1/2 EXPECTATION STEP: updates conditional expectations given last parameters 
   for(h in 1:H){
      mit.h <- list(p=1,mu=t(matrix(mit$mu[h,])),
                    Sigma = t(matrix(mit$Sigma[h,])),
                    df = mit$df[h])
      z[,h] <- exp(log(mit$p[h]) + dmvgt(theta,mit.h,log=TRUE))     
   } 
   z = z / matrix(rowSums(z),N,H)  # normalize indicator probability
   # calculate digamma functions for equation (11)
   psi1 <- digamma((k+mit$df)/2)
   psi2 <- digamma(mit$df/2)
   for(h in 1:H){
      # calc rho for eq (10) #
      mu.mat  <- matrix(mit$mu[h,],N,k,byrow=TRUE)
      tmp     <- chol(chol2inv(chol(matrix(mit$Sigma[h,],k,k)))) 
      tmp      = tmp %*% t(theta - mu.mat)  
      rho[,h]  = colSums(tmp^2)
      # calc wg.ln <- tilde (log(w)) <- E(log(w|theta)) eq(11)#
      tmp.p   <- cbind(z[,h],(1-z[,h]))
      tmp.m   <- cbind( log((rho[,h]+mit$df[h])/2),log(mit$df[h]/2))
      tmp.psi <- matrix(c(psi1[h],psi2[h]),N,2,byrow=TRUE)
      tmp.m    = tmp.m - tmp.psi
      wg.ln[,h]= rowSums(tmp.m * tmp.p)
   }
   # calc weighted membership eq(10) #
   z.wg =  z * (k + dfmat) / (dfmat + rho)
   # calc delta eq(12) using equality in (10)#
   delta <- z.wg  + (1-z)

   # 2/2 MAXIMIZATION STEP: updates mixture parameters using IS weights 
   Sigma <- matrix(0,nrow = H, ncol=(k*k))
   mu <- matrix(0,nrow = H, ncol=k)
   for(h in 1:H){
      tmp.wg    = w * z.wg[,h] / sum(w * z.wg[,h])
      tmp       = matrix(tmp.wg,N,k)
      mu[h,]    <- colSums(tmp*theta)      
      tmp.theta <- theta - matrix(mu[h,],N,k,byrow=TRUE)
      tmpSigmas <- t(apply(tmp.theta,1,function(x)(c(x %o% x))))
      tmp.wg    <- w * z.wg[,h] / sum(w * z[,h])
      tmp.wg    = matrix(tmp.wg,N,(k^2))
      Sigma[h,] = colSums(tmpSigmas * tmp.wg)
   }
   # update eta, mixture probability
   eta <- colSums(w.mat * z) / sum(w)
   # update degrees of freedom if optimized
   df = mit$df # df not altered if not optimized
   if(control$optim.df) {
      # calculate weigted expectations in equation (16)
      w.exp <- colSums(w.mat * (wg.ln + delta)) / sum(w)
      # (increasing) objective function
      fn_sub <- function(df,w.exp){
         r = - digamma(df/2) + log(df/2) + 1 - w.exp
         (r)
      }
      df.new <- df
      for (h in 1:H){
      # check if search interval is correct
      tmp <- as.matrix(control$inter.df)
      f.bounds <- apply(tmp, 1, FUN=fn_sub,w.exp[h])
      if(all(f.bounds>0)){
        df.new[h] = control$inter.df[2]  
      }else if(all (f.bounds<0)){
        df.new[h] = control$inter.df[1]
      }else{
        optim.df <- uniroot(f=fn_sub,interval=control$inter.df,tol=control$tol.df,maxiter=control$maxit.df,
                            w.exp=w.exp[h])
        df.new[h] <- optim.df$root
      }
      }
   df=df.new
   }
   mit <- list(p=eta,mu=mu,Sigma=Sigma,df=df)
   return(list(mit=mit))
}


