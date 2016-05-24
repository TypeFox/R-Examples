wt <-
function(x, start = 1, dt = 1, dj = 1/20, 
         lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3),
         make.pval = T, method = "white.noise", params = NULL, 
         n.sim = 100, save.sim = F) {
 
                                 
  ###############################################################################
  ## Call function WaveletTransform
  ## Retrieve the wavelet transform, power, phases, amplitudes
  ###############################################################################
  
  # wavelet transform
  WT = WaveletTransform(x, dt = dt, dj = dj, 
                        lowerPeriod = lowerPeriod, upperPeriod = upperPeriod)
          
  Wave  = WT$Wave 
  Phase = WT$Phase
  Ampl  = WT$Ampl

  Power = WT$Power
  Power.avg = rowMeans(Power)
  
  Period = WT$Period
  Scale  = WT$Scale
  nr  = WT$nr
  nc  = WT$nc
  
  rm(WT)

  ###############################################################################
  ## Compute p values for significance check
  ###############################################################################
    
  Power.pval = NULL
  Power.avg.pval = NULL
  series.sim = NULL
  
  if (make.pval == T) {
       
      Power.pval = matrix(0, nrow = nr, ncol = nc)
      Power.avg.pval = rep(0, nr)
      
      if (save.sim == T) { series.sim = matrix(NA, nrow=nc, ncol=n.sim) }
       
      pb = txtProgressBar(min = 0, max = n.sim, style = 3) # create a progress bar
      for(ind.sim in 1:n.sim){
      
          x.sim = SurrogateData(x, method = method)
          
          if (save.sim == T) { series.sim[,ind.sim] = x.sim }
          
          WT.sim = WaveletTransform(x.sim, dt = dt, dj = dj, 
                                    lowerPeriod = lowerPeriod, upperPeriod = upperPeriod)
                                   
          Power.sim = WT.sim$Power
          Power.avg.sim = rowMeans(Power.sim)  
          
          rm(WT.sim)
          
          Power.pval[Power.sim >= Power] = Power.pval[Power.sim >= Power] + 1
          Power.avg.pval[Power.avg.sim >= Power.avg] = Power.avg.pval[Power.avg.sim >= Power.avg] + 1
          setTxtProgressBar(pb, ind.sim) # set progress bar
      }
      close(pb) # close progress bar

      # p-values
      
      Power.pval = Power.pval / n.sim 
      Power.avg.pval = Power.avg.pval / n.sim  
         
  }  
  
  ###############################################################################
  ## Compute the cone of influence COI
  ###############################################################################

  coi = COI(start = start, dt = dt, nc = nc, nr = nr, Period = Period)
  
  ###############################################################################
  ## Prepare the output
  ###############################################################################

  output = list(Wave = Wave, Phase = Phase, Ampl = Ampl,
                Power = Power, Power.avg = Power.avg,
                Power.pval = Power.pval, Power.avg.pval = Power.avg.pval, 
                Period = Period, Scale = Scale,     
                coi.1 = coi$x, coi.2 = coi$y,
                nc = nc, nr = nr,    
                axis.1 = coi$axis.1, axis.2 = coi$axis.2,
                series.sim = series.sim)
                
  return(invisible(output))
}
