rnum <- function(x0, t, f, p, dt, n){

 # Numerical solution of random model
 # arguments: x0 init cond
 #            t times for output
 #            f model
 #            p parameters
 #            dt time step
 #            n number of realizations
 # arrays to store results
    nt <- length(t); nx <- length(x0)
    ns <- floor((t[2]-t[1])/dt)
    X <- structure(c(1:(nt*nx*n)),dim=c(nt,nx,n));X[,,]<-NA
    Xm <- matrix(nrow=nt,ncol=nx); Xs <- Xm
 # loops
   # realizations
   for(k in 1:n){
 # first value of X is initial condition 
    Xt <- x0; tt <- t[1]
    for(i in 1:nt){
     X[i,,k] <- Xt
     for(j in 1:ns){
      tt <- tt + dt;Xt <- Xt + dt*f(tt,p,Xt)
      # X assumed positive
      if(Xt<0) Xt <- 0
     }
    } # end of one realization
   } # end of run

 # statistics for the run
      for(i in 1:nt){
       for(j in 1:nx){
        Xm[i,j] <- mean(X[i,j,])
        Xs[i,j] <- sd(X[i,j,])
       }
      }
 return(cbind(Xm,Xs))
} # end of function


