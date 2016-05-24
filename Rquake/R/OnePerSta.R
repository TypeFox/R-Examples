OnePerSta<-function(twpx, phase="Y")
  {

####  given a pick dataframe
###  return a cleaned list with one pick per station
###   given a specific phase


    nona = which( is.na(twpx$tag) )
    
    if(length(nona)>0)
      {
        twpx = RSEIS::deleteWPX(twpx, nona)
      }
    
    
    A1T = Qrangedatetime(twpx)
    s1 = RSEIS::secdifL(A1T$min,  twpx)
    
    usta = unique(twpx$name)

    for(i in 1:length(usta))
      {
        jsta =  usta[i]
        ksta = which( twpx$name ==  jsta &  twpx$phase == phase)

        w1 =  which.min( s1[ksta]  )
        keep =  ksta[w1]

        
        twpx$onoff[ ksta[-w1] ] = -2
      }

    ww = which( twpx$onoff<0 )
    
    twpx = RSEIS::deleteWPX(twpx, ww)
    
    
    return(twpx)


    
  }
