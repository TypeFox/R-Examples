Y2Pphase<-function(twpx, phase)
  {
    ##  convert all named  phases to P
    ##  if there are P and other  picks for the same station, use P, discard other
    WY  = which(twpx$phase==phase)
    if(length(WY)<1)  return(twpx)
        WP  = which(twpx$phase=="P")
    if(length(WP)>0)
      {
       mpy =  match(  twpx$name[WY], twpx$name[WP]  )
       wrid = which(!is.na(mpy))
       if(length(wrid)>0)
         {
           twpx =   RSEIS::deleteWPX(twpx,wrid )
           WY  = which(twpx$phase==phase)
           twpx$phase[WY]  = "P"
         }
     }
    else
      {
        twpx$phase[WY]  = "P"
      }
    
    Awpx = twpx

    
    return(Awpx)
  }
