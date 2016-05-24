WaveletCoherency <-
function(x, y, dt = 1, dj = 1/20, 
         lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3),
         window.type.t=1, window.type.s=1, window.size.t=5, window.size.s=1/4){
   
  #############################################################################                              
  ## compute wavelet transforms, retrieve results
  #############################################################################
  
  WT = WaveletTransform(x, dt = dt, dj = dj, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod)
  
  Period = WT$Period
  Scale  = WT$Scale
  nc = WT$nc
  nr = WT$nr  
          
  Wave.x  = WT$Wave 
  Phase.x = WT$Phase
  Ampl.x  = WT$Ampl
  Power.x = WT$Power
  
  WT = WaveletTransform(y, dt = dt, dj = dj, lowerPeriod = lowerPeriod, upperPeriod = upperPeriod)
  
  Wave.y  = WT$Wave 
  Phase.y = WT$Phase
  Ampl.y  = WT$Ampl
  Power.y = WT$Power
    
  #############################################################################                              
  ## compute cross-wavelet transform and cross-wavelet power
  #############################################################################
   
  Wave.xy  = (Wave.x * Conj(Wave.y)) / matrix(rep(Scale, nc), nrow=nr)
  Power.xy = Mod(Wave.xy)
  
  ############################################################################
  ## Windows to use for smoothing cross-wavelet transform and individual wavelet power
  ## Inspired by:
  ## Luis Aguiar-Conraria and Maria Joana Soares, "GWPackage"
  ############################################################################
  
  window.func = function(type=0, n=5) {
  
      if  ( (n!=floor(n)) ||  n<=0 )
            {  stop("Smoothing requires a positive and integer window size!") }
            
      # no smoothing
        if  (is.element(type, c(0,'none'))) {
             window = 1
        }         
            
      # Bartlett  window   (triangular window with L=n-1)       
        if  (is.element(type, c(1,'bar'))) {         
             if (n<=2) { stop("Bartlett window requires minimum size 3!") }
             window = 1 - abs((0:(n-1)) - (n-1)/2)/((n-1)/2)                          
        } 
        
      # Triangular window  (L=n)
        if  (is.element(type, c(2,'tri'))) {
             if (n==1) { stop("Triangular (non-Bartlett) window requires minimum size 2!") }  
             window = 1 - abs((0:(n-1)) - (n-1)/2)/(n/2)
        }  
          
      # Rectangular  window  (Boxcar or Dirichlet window)
        if  (is.element(type, c(3,'box'))) {
             window  = rep(1,n)
        }         
                                
      # Hanning  window  
        if  (is.element(type, c(4,'han'))) {
             if (n<=2) { stop("Hanning window requires minimum size 3!") }
             window  = 0.5 - 0.5*cos( 2*pi*( 0:(n-1) )/(n-1) )
                }
                
      # Hamming  window  
        if  (is.element(type, c(5,'ham'))) {
             if (n==1) { stop("Hamming window requires minimum size 2!") }
             window  = 0.53836 - (1-0.53836)*cos( 2*pi*( 0:(n-1) )/(n-1) )  
        }     
                
      # Blackmann  window  
        if  (is.element(type, c(6,'bla'))) {
             if (n==1) { stop("Blackman window requires minimum size 2!") }
             window  = 7938/18608 - (9240/18608)*cos( 2*pi*(0:(n-1))/(n-1) ) + (1430/18608)*cos( 4*pi*(0:(n-1))/(n-1) )
        }          

      # Normalization   
      window =  window/sum(window) 
      return(window)
    }



  #############################################################################
  ## computation of wavelet coherence
  #############################################################################

  # odd window sizes, given in terms of dj, dt resolution
  window.size.s = 2*floor(window.size.s/(2*dj))+1
  window.size.t = 2*floor(window.size.t/(2*dt))+1

  # window 2D matrix
  window2D = window.func(window.type.s, window.size.s) %*% t(window.func(window.type.t, window.size.t))

  # function to smooth a matrix according to a given smoothing window
  smooth2D = function(mat, window2D) {
      mat.pad = matrix(0, nrow = nrow(mat) + nrow(window2D)-1, ncol = ncol(mat) + ncol(window2D)-1)
      window2D.pad = mat.pad
      mat.pad[1:nrow(mat),1:ncol(mat)] = mat
      window2D.pad[1:nrow(window2D),1:ncol(window2D)] = window2D
      smooth.mat = fft(fft(mat.pad)* fft(window2D.pad), inverse = TRUE)/length(mat.pad)
      return( smooth.mat[ floor(nrow(window2D)/2) + (1:nrow(mat)), floor(ncol(window2D)/2) + (1:ncol(mat)) ] )
  }    

  # smooth individual wavelet powers and smooth the cross-wavelet transform
  sPower.x = Re(smooth2D(Power.x, window2D))
  sPower.y = Re(smooth2D(Power.y, window2D))
  sWave.xy = smooth2D(Wave.xy, window2D)
          
  # wavelet coherency   
  Coherency = sWave.xy / sqrt(sPower.x*sPower.y)
#   Coherency[which(sPower.x < 1e-4)] = 0
#   Coherency[which(sPower.y < 1e-4)] = 0
  Coherency[which(sPower.x*sPower.y == 0)] = 0
    
  # Wavelet coherence
  Coherence = Mod(sWave.xy)^2 / (sPower.x*sPower.y) 
#   Coherence[which(sPower.x < 1e-4)] = 0
#   Coherence[which(sPower.y < 1e-4)] = 0
  Coherence[which(sPower.x*sPower.y == 0)] = 0
  
  #############################################################################                              
  ## prepare the output
  #############################################################################
 
  output <- list(Wave.xy = Wave.xy, sWave.xy = sWave.xy,  
                 Power.xy = Power.xy, 
                 Coherency = Coherency, Coherence = Coherence, 
                 Wave.x = Wave.x, Wave.y = Wave.y, 
                 Phase.x = Phase.x, Phase.y = Phase.y,
                 Ampl.x = Ampl.x, Ampl.y = Ampl.y,
                 Power.x = Power.x, Power.y = Power.y,
                 sPower.x = sPower.x, sPower.y = sPower.y,
                 Period = WT$Period, Scale = WT$Scale, 
                 nc = nc, nr = nr)
                 
  return(invisible(output))
}
