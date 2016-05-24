ramos <- function(x0, t, f, p, dt){

 # Ramos nonstandard explicit integration algorithm (EIA)
 # arguments: x0 initial condition
 #            t times for output
 #            f model
 #            df derivative of model
 #            p parameters
 #            dt time step
 # arrays to store results
    nt <- length(t); nx <- length(x0)
    ns <- floor((t[2]-t[1])/dt)
    X <- matrix(nrow=nt,ncol=nx)
 # first value of X is initial condition 
    Xt <- x0; tt <- t[1]
 # loops
   for(i in 1:nt){
     X[i,] <- Xt
     for(j in 1:ns){
      tt <- tt + dt
      num <- 2*dt*f(tt,p,Xt)[1]
      den <- 2-dt*f(tt,p,Xt)[2]
      Xt <- Xt + num/den
      # X assumed positive
      for(k in 1:nx) if(Xt[k]<0) Xt[k] <- 0
     }
   } # end of run
 return(X)
} # end of function


