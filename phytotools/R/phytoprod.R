
phytoprod <- function(PE, Ein, kpar,
                      cz=matrix(data=c(1,1),ncol=2),
                      zmax=NA){
  
  ##########################
  #Initialize Time and Depth
  
  #Calculate number of days to simulate 
  dayofyear <- unique(floor(Ein[,1]))
  dayofyear <- dayofyear[1:(length(dayofyear)-1)]
  n         <- length(dayofyear)
  
  #Determine time steps
  ts <- length(Ein[,1])
  
  #Calculate intervals per day
  interval <- (ts-1)/n
  
  #Determine euphotic depth as 0.5% of light penetration (-ln(0.005) = 5.3) 
  zeu <- ceiling(5.3/kpar)
  
  #If zmax is not specified or deeper than zeu, then set zmax = zeu
  if (is.na(zmax) | (is.finite(zmax) & zmax > zeu)) {
    zmax <- zeu              
  }
  
  #Setup all arrays
  zseq <- seq(0,zmax, length.out = 50)  #Depths over which to calculate E and P
  E    <- array(NA, c(ts,50))           #E is irradiance through time and depth
  P    <- array(NA, c(ts,50))           #P is productivity through time and depth
  Pz   <- rep(NA, ts)                   #Depth integrated P
  PP   <- rep(NA, n)                    #Daily areal productivity
  
  #####################################
  #Calculate E(z,t) = E0 * exp(-kpar*z)

  for (t in 1:ts){
    E[t,] <- Ein[t,2] * exp(-1*kpar*zseq) 
  }
  
  #################################
  #Calculate model dependent P(t,z)
  
  if (PE$model=="Webb")
  {P <- PE$alpha[1]*PE$ek[1]*(1-exp(-1*E/PE$ek[1]))}
  
  if (PE$model=="JP")  
  {P <- PE$alpha[1]*PE$ek[1]*tanh(E/PE$ek[1])}
  
  if (PE$model=="PGH")  
  {P <- PE$ps[1]*(1-exp(-1*PE$alpha[1]*E/PE$ps[1]))*exp(-1*PE$beta[1]*E/PE$ps[1])}
  
  if (PE$model=="EP") 
  {P <- E/((1/(PE$alpha[1]*PE$eopt[1]^2))*E^2+(1/PE$ps[1]-2/(PE$alpha[1]*PE$eopt[1]))*E+(1/PE$alpha[1]))}

  #########################################################
  #Interpolate biomass through depth, calculate P x biomass 
  
  if (length(cz[,2])>1){
    c.approx <- approx(x=cz[,2], y=cz[,1], xout=zseq, rule=2)$y
  }else{
    c.approx <- rep(cz[,2], 50)
  }
  
  #Multiply P by to depth dependent concentration
  P        <- sweep(P, MARGIN=2, c.approx, `*`)

  #################################################
  #Unit Conversion if PE data is quantum efficiency

  if (PE$normalize==T){
    P <- P * 1e-3 * 3600 #convert umol m-3 s-1 to mmol m-3 hr-1
  }
  
  ###################################
  #Integrate phytoplankton production 
  
  #Integrate through depth 
  for (t in 1:ts){
    fn       <- splinefun(zseq, P[t,])                            
    Pz[t]    <- integrate(fn, 0, zmax)$value
  }
  
  #Integrate each day
  for (i in 1:n){    #Loop through n days
    index    <- seq((i-1) * interval+1, (i*interval+1), by=1)
    fn       <- splinefun(seq(0,1,length.out=interval+1), Pz[index])                            
    PP[i]    <- round(integrate(fn,0,1)$value*24, digits=3)
  }
  
  #####################
  #Return relevant data
  
  return(list(PP=cbind(dayofyear,PP), z=zseq, t=Ein[,1], P=P))

}
