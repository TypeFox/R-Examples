vectorsetup <- 
  function(obj) 
# obj is an inputlist class object yielding all the input parameters    
    {

# setup of vectors for the mean and the two covariance factors of the diffusion process  

dt2 <- 0.5*obj$deltat
ntime <- floor((obj$Tfin - obj$t0)/obj$deltat)  

time.points <- seq(obj$t0,obj$Tfin,obj$deltat)

aa <-a(time.points)
app <- (aa[1:ntime]+aa[2:(ntime+1)])*dt2
app <- cumsum(app)
rm(aa)

vp <- exp(app)  
vp <- c(1.0, vp)
#
aa<-b(time.points)/vp
mp <- (aa[1:ntime]+aa[2:(ntime+1)])*dt2
mp <- c (0, cumsum(mp))
rm(aa)

mp <- vp*(obj$x0 + mp) 


aa <- cc(time.points)/(vp^2)
up <- (aa[1:ntime]+aa[2:(ntime+1)])*dt2
up <- c(0, cumsum(up))
rm(aa)

up <- up*vp 
answer <- cbind(mp,up,vp)
return(answer)
}
