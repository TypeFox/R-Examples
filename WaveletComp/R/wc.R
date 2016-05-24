wc <-
function(x, y, start = 1,
         dt = 1, dj=1/20, 
         lowerPeriod = 2*dt, upperPeriod = floor(length(x)/3)*dt,
         window.type.t=1, window.type.s=1, window.size.t=5, window.size.s=1/4,
         make.pval = T,
         method = "white.noise",
         params = NULL,
         n.sim = 100, save.sim = F){

  ###############################################################################
  ## Call function WaveletCoherency
  ## Retrieve the wavelet coherency, coherence and cross-wavelet transform and power, 
  ## as well as individual wavelet transforms, power, phases, amplitudes
  ###############################################################################       
  # computation of wavelet coherence of x over y
  
  WC =  WaveletCoherency(x, y, dt = dt, dj = dj, 
                         lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
                         window.type.t=window.type.t, window.type.s=window.type.s, window.size.t=window.size.t, window.size.s=window.size.s)                          
       
  # cross-wavelet and coherency results  
  
  Wave.xy  = WC$Wave.xy
  sWave.xy = WC$sWave.xy
  
  Angle  = Arg(Wave.xy)
  sAngle = Arg(sWave.xy)
  
  Power.xy = WC$Power.xy
  Power.xy.avg = rowMeans(Power.xy)  
   
  Coherency = WC$Coherency
  Coherence = WC$Coherence  
  Coherence.avg = rowMeans(Coherence)  
  
  # individual wavelet results
  
  Wave.x = WC$Wave.x
  Wave.y = WC$Wave.y
  Phase.x = WC$Phase.x
  Phase.y = WC$Phase.y
  Ampl.x = WC$Ampl.x
  Ampl.y = WC$Ampl.y 
  
  Power.x = WC$Power.x
  Power.y = WC$Power.y
  Power.x.avg = rowMeans(Power.x)
  Power.y.avg = rowMeans(Power.y)
  
  sPower.x = WC$sPower.x
  sPower.y = WC$sPower.y
  
  # parameters
  
  Period = WC$Period
  Scale  = WC$Scale
  nr  = WC$nr
  nc  = WC$nc 
  
  rm(WC)
  
  ###############################################################################
  ## Compute p values for significance check
  ###############################################################################
   
  Coherence.pval = NULL 
  Power.xy.pval  = NULL
  Power.x.pval   = NULL
  Power.y.pval   = NULL
  
  Coherence.avg.pval = NULL
  Power.xy.avg.pval  = NULL 
  Power.x.avg.pval   = NULL
  Power.y.avg.pval   = NULL
  
  series.sim = NULL
  
  x.sim = NULL
  y.sim = NULL
  
  if (make.pval == T) {
  
      Coherence.pval     =  matrix(0, nrow = nr, ncol = nc)
      Power.xy.pval      =  matrix(0, nrow = nr, ncol = nc)
      Power.x.pval       =  matrix(0, nrow = nr, ncol = nc)
      Power.y.pval       =  matrix(0, nrow = nr, ncol = nc)
      
      Coherence.avg.pval =  rep(0, nr)
      Power.xy.avg.pval  =  rep(0, nr)
      Power.x.avg.pval   =  rep(0, nr)
      Power.y.avg.pval   =  rep(0, nr)          
      
      if (save.sim == T) { 
          series.sim$x = matrix(NA, nrow=nc, ncol=n.sim)
          series.sim$y = matrix(NA, nrow=nc, ncol=n.sim)
      }    
  
      pbar = txtProgressBar(min = 0, max = n.sim, style = 3) # create a progress bar
      for(ind.sim in 1:n.sim){
      
          x.s = SurrogateData(x, method = method, params = params)
          y.s = SurrogateData(y, method = method, params = params)
          
          if (save.sim == T) { 
              series.sim$x[,ind.sim] = x.s
              series.sim$y[,ind.sim] = y.s
          }    
          
          WC.s = WaveletCoherency(x.s, y.s, dt = dt, dj = dj,
                                  lowerPeriod = lowerPeriod, upperPeriod = upperPeriod,
                                  window.type.t=window.type.t, window.type.s=window.type.s, window.size.t=window.size.t, window.size.s=window.size.s)   
                                  
          Coherence.s     = WC.s$Coherence
          Power.xy.s      = WC.s$Power.xy
          Power.x.s       = WC.s$Power.x
          Power.y.s       = WC.s$Power.y
          
          Coherence.avg.s = rowMeans(Coherence.s)
          Power.xy.avg.s  = rowMeans(Power.xy.s)
          Power.x.avg.s   = rowMeans(Power.x.s)
          Power.y.avg.s   = rowMeans(Power.y.s)
          
          rm(WC.s)
          
          Coherence.pval[which(abs(Coherence.s) >= abs(Coherence))] = Coherence.pval[which(abs(Coherence.s) >= abs(Coherence))] + 1
          Power.xy.pval[which(abs(Power.xy.s) >= abs(Power.xy))] = Power.xy.pval[which(abs(Power.xy.s) >= abs(Power.xy))] + 1
          Power.x.pval[which(abs(Power.x.s) >= abs(Power.x))] = Power.x.pval[which(abs(Power.x.s) >= abs(Power.x))] + 1
          Power.y.pval[which(abs(Power.y.s) >= abs(Power.y))] = Power.y.pval[which(abs(Power.y.s) >= abs(Power.y))] + 1
          Coherence.avg.pval[which(abs(Coherence.avg.s) >= abs(Coherence.avg))] = Coherence.avg.pval[which(abs(Coherence.avg.s) >= abs(Coherence.avg))] + 1
          Power.xy.avg.pval[which(abs(Power.xy.avg.s) >= abs(Power.xy.avg))] = Power.xy.avg.pval[which(abs(Power.xy.avg.s) >= abs(Power.xy.avg))] + 1
          Power.x.avg.pval[which(abs(Power.x.avg.s) >= abs(Power.x.avg))] = Power.x.avg.pval[which(abs(Power.x.avg.s) >= abs(Power.x.avg))] + 1
          Power.y.avg.pval[which(abs(Power.y.avg.s) >= abs(Power.y.avg))] = Power.y.avg.pval[which(abs(Power.y.avg.s) >= abs(Power.y.avg))] + 1
          
          setTxtProgressBar(pbar, ind.sim) # set progress bar
      }
      close(pbar) # close progress bar
      
      # p-values

      Coherence.pval     = Coherence.pval/n.sim
      Power.xy.pval      = Power.xy.pval/n.sim 
      Power.x.pval       = Power.x.pval/n.sim
      Power.y.pval       = Power.y.pval/n.sim
      
      Coherence.avg.pval = Coherence.avg.pval/n.sim
      Power.xy.avg.pval  = Power.xy.avg.pval/n.sim
      Power.x.avg.pval   = Power.x.avg.pval/n.sim
      Power.y.avg.pval   = Power.y.avg.pval/n.sim
        
  }    
  

  ###############################################################################
  ## Compute the cone of influence COI
  ###############################################################################

  coi <- COI(start = start, dt = dt, nc = nc, nr = nr, Period = Period)

  ###############################################################################
  ## Prepare the output
  ###############################################################################
    
  output <- list(Wave.xy = Wave.xy, Angle = Angle,
                 sWave.xy = sWave.xy, sAngle = sAngle,     
                 Power.xy = Power.xy, Power.xy.avg = Power.xy.avg, 
                 Power.xy.pval = Power.xy.pval, Power.xy.avg.pval = Power.xy.avg.pval, 
                 Coherency = Coherency,
                 Coherence = Coherence, Coherence.avg = Coherence.avg,
                 Coherence.pval = Coherence.pval, Coherence.avg.pval = Coherence.avg.pval,
                 Wave.x = Wave.x, Wave.y = Wave.y,
                 Phase.x = Phase.x, Phase.y = Phase.y,
                 Ampl.x = Ampl.x, Ampl.y = Ampl.y,             
                 Power.x = Power.x, Power.y = Power.y, 
                 Power.x.avg = Power.x.avg, Power.y.avg = Power.y.avg,
                 Power.x.pval = Power.x.pval, Power.y.pval = Power.y.pval, 
                 Power.x.avg.pval = Power.x.avg.pval, Power.y.avg.pval = Power.y.avg.pval,
                 sPower.x = sPower.x, sPower.y = sPower.y,
                 Period = Period, Scale = Scale,
                 coi.1 = coi$x, coi.2 = coi$y, 
                 nc = nc, nr = nr,
                 axis.1 = coi$axis.1, axis.2 = coi$axis.2,
                 series.sim = series.sim)
  return(invisible(output))
}
