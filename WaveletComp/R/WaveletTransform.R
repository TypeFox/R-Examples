WaveletTransform <-
function(x, dt = 1, dj = 1/20, 
         lowerPeriod = 2*dt, upperPeriod = floor(length(x)*dt/3)){
 
  ###############################################################################
  ## Provide parameters (which could be useful for other transforms as well)
  ###############################################################################
  
  # Original length and length of zero padding:
  series.length = length(x)
  pot2 = trunc(log2(series.length) + 0.5)
  pad.length = 2^(pot2+1)-series.length 

  # Define central angular frequency omega0 and fourier factor:
  omega0 = 6
#   fourier.factor   = (4*pi)/(omega0 + sqrt(2+omega0^2))
  fourier.factor = (2*pi)/omega0
  
  # Compute scales and periods:
  min.scale = lowerPeriod/fourier.factor             # Convert lowerPeriod to minimum scale 
  max.scale = upperPeriod/fourier.factor             # Convert upperPeriod to maximum scale 
  J = as.integer( log2(max.scale/min.scale) / dj)    # Index of maximum scale -1
    
  scales = min.scale * 2^((0:J)*dj)        # sequence of scales 
  scales.length = length(scales)           # J + 1
  periods = fourier.factor*scales          # sequence of periods

  # Computation of the angular frequencies
  N = series.length+pad.length
  omega.k = 1:floor(N/2)
  omega.k = omega.k * (2*pi)/(N*dt)                    # k <= N/2
  omega.k = c(0, omega.k, -omega.k[ floor((N-1)/2):1 ])

  ###############################################################################
  ## Define the Morlet wavelet transform function
  ###############################################################################
  
  morlet.wavelet.transform = function(x) {

    # Standardize x and pad with zeros
    x = (x-mean(x))/sd(x)     
    xpad = c(x, rep(0,pad.length)) 
    
    # Compute Fast Fourier Transform of xpad
    fft.xpad = fft(xpad)
    
    # Compute wavelet transform of x 
    # Prepare a complex matrix which accomodates the wavelet transform
    wave = matrix(0, nrow=scales.length, ncol=N)                 
    wave = wave + 1i*wave                   
    
    # Computation for each scale...
    # ... simultaneously for all time instances
    for (ind.scale in (1:scales.length)) {
        
         my.scale = scales[ind.scale]
         
         norm.factor = pi^(1/4) * sqrt(2*my.scale/dt)
         expnt       = -( (my.scale * omega.k - omega0)^2 / 2 ) * (omega.k > 0)
         daughter    = norm.factor * exp(expnt)
         daughter    = daughter * (omega.k > 0)
         
         wave[ind.scale,] = fft( fft.xpad * daughter, inverse=TRUE) / N
    }
       
    # Cut out the wavelet transform
    wave = wave[,1:series.length]
    return(wave)
  }        
  
  ###############################################################################
  ## Compute the wavelet transform, power, phases, amplitudes
  ###############################################################################
 
  Wave = morlet.wavelet.transform(x) 

  # Compute wavelet power
  Power = Mod(Wave)^2 / matrix(rep(scales, series.length), nrow=scales.length)
  
  # Phase  
  Phase = Arg(Wave)
  
  # Amplitude
  Ampl  = Mod(Wave) / matrix(rep(sqrt(scales), series.length), nrow=scales.length)
  
  ###############################################################################
  ## Prepare the output
  ###############################################################################

  output = list(Wave = Wave, 
                Phase = Phase, Ampl = Ampl,
                Period = periods, Scale = scales,
                Power = Power, 
                nc = series.length, nr = scales.length)
                 
  return(invisible(output))
}
