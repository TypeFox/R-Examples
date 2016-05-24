RK4D <- function(x0, t, f, p, dt, g, tz){

 # RK4 numerical solution of model
 # arguments: x0 init cond
 #            t times for output 
 #            f continuous model
 #            p parameters
 #            g discontinuous model
 #            tz discontinuities
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
      k1<- f(tt,p,Xt)
      k2<- f(tt+dt/2,p,(Xt+k1*dt/2))
      k3<- f(tt+dt/2,p,(Xt+k2*dt/2))
      k4<- f(tt+dt,p,(Xt+k3*dt))
      kavg <- (dt/6)*(k1+ 2*k2+ 2*k3+ k4)
      tt <- tt + dt; Xt <- Xt + kavg
      Xt <- Xt + g(tt,p,Xt,tz)
     }
  } # end of run
 return(X)
} # end of function

