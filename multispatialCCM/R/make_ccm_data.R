make_ccm_data <-
function(sp_sd=0.125, obs_sd=0.025, Sstr=0.375, times=10, burnin=100, number_of_chains=20, seednum=2718) {
  #Prepare vectors to create simulated data
  #sp_sd, Standard deviation used to add process noise
  #obs_sd, Standard deviatoin used to add observation error
  #Sstr, Forcing strength of process
  
  #times, How many sequential observations in each plot?
  #burnin, Burnin time before starting the experiment.
  #Is used to remove correlation among plots that occurs because of starting conditions
  
  #number_of_chains, Number of plots in spatially replicated data
  
  #Create emtpy vectors
  Accm<-NULL
  Bccm<-NULL
  time_ccm<-NULL
  
  #Set seed for generating random numbers. Remove this if you want
  #distinct results each time you run the script
  set.seed(seednum)
  
  #The below loop will populate Accm and Bccm. Accm is forced by Bccm based on
  #the forcing strength SStr. Bccm is not forced by Accm.
  #The output for all plots is stored as a single long vector.
  #Plots are separated by a single "NA" record
  for(i in 1:number_of_chains) { #Build plots
    
    #Add process error to species intrinsic growth rates
    #This will cause species to have slightly different
    #growth rates in different plots
    rx<-3.8+rnorm(1,0,sp_sd)
    ry<-3.6+rnorm(1,0,sp_sd)
    
    #Simulate competition
    #Note - y influences x here, but x does not influence y
    x<-numeric(length(times)+burnin)
    y<-numeric(length(times)+burnin)
    x[1]<-runif(1); y[1]<-runif(1)
    for(ii in 2:max(times+burnin)) {
      x[ii]<-x[ii-1]*(rx-rx*x[ii-1]-Sstr*y[ii-1])
      y[ii]<-y[ii-1]*(ry-ry*y[ii-1])
      
      #Reflecting boundaries for population >1 or <0
      while(x[ii]<0|x[ii]>1) {
        if(x[ii]<0) {x[ii]<-(-x[ii])}
        if(x[ii]>1) {x[ii]<-1+(1-x[ii])}
      }
      
      while(y[ii]<0|y[ii]>1) {
        if(y[ii]<0) {y[ii]<-(-y[ii])}
        if(y[ii]>1) {y[ii]<-1+(1-y[ii])}
      }        
    }
    
    #Remove data from burnin time
    x<-x[-c(1:burnin)]
    y<-y[-c(1:burnin)]
    
    #Add multiplicative observation error
    x<-x*rlnorm(length(x), 0, obs_sd)
    y<-y*rlnorm(length(y), 0, obs_sd)
    
    #Save variables
    Accm<-c(Accm, NA, x)
    Bccm<-c(Bccm, NA, y)
    time_ccm<-c(time_ccm, NA, 1:times)
    
  } # End plot building
  return(list(Accm=Accm, Bccm=Bccm, time_ccm=time_ccm))
}
