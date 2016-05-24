anthem<-function(){
   requireNamespace("tuneR")
   
   ####notes encodees comme suit: la3 = 0, si3 = 2, do4 = 3, etc. (un cran = un demi-ton) 
   keys = c(-2, 3, -2, 0, 2, -5, -5, 0, -2, -4, -2,
            -9, -9, -7, -7, -5, -4, -4, -2, 0, 2, 3, 5)
   freqs = 220*2^(keys/12) ##gamme temperee car faut pas deconner non plus...
   
   ####durees des notes encodees comme suit: ronde = 4, blanche = 2, noire = 1,... 
   durations = c(1/2, 1, 3/4, 1/4, 1, 1/2, 1/2, 1, 3/4, 1/4, 
                 1, 1/2, 1/2, 1, 1/2, 1/2, 1, 1/2, 1/2, 1, 1/2, 1/2, 1)
   
   BPM = 60
   samp.rate = 44100
   
   anthem = tuneR::bind(tuneR::sine(freq=freqs[1], duration = durations[1]*(60/BPM)*(samp.rate), samp.rate = samp.rate), 
                        tuneR::sine(freq=freqs[2], duration = durations[2]*(60/BPM)*(samp.rate), samp.rate = samp.rate))
   
   for (t in 3:length(keys))
   {
      anthem = tuneR::bind(anthem, tuneR::sine(freq=freqs[t], duration = durations[t]*(60/BPM)*(samp.rate), samp.rate = samp.rate))
   }
   
   correg()
   
   
   tuneR::writeWave(anthem, filename="russianthem.wav", extensible = FALSE)
   
   tuneR::play(anthem)
   return("CorReg")
}