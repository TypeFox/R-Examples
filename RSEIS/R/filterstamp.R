filterstamp = function(fl=1/2, fh=10, type="BP")
  {
#####  normally fl fh and type are vectors of filters
    ###  however, if fl is a list, then the
    ###  information is taken from the elements of the list

    
    if(missing(fh)) { fh = fl }
    if(missing(type)) { type = "LP" }
    
    if(is.list(fl))
      {
        
        fh = fl$fh
        type= fl$type

        fl = fl$fl
      }


    
    n=max(c(length(fl), length(fh)))
    Notes = as.vector(1:(n))


    
    if(length(fl)<n)
      {
        fl = rep(fl[1], n)

      }
    if(length(fh)<n)
      {
        fh = rep(fh[1], n)

      }
    if(length(type)<n)
      {
        type = rep(type[1], n)

      }
    
    for(i in 1:n)
      {

        
        khigh = format.default(fh[i], digits=3)
        lhigh  = "Hz"
        if(fh[i]<1)
          {
            khigh = format.default(1/fh[i], digits=3)
            lhigh  = "s"
          }
        klow = format.default(fl[i], digits=3)
        llow  = "Hz"
        
        if(fl[i]<1)
          {
            klow = format.default(1/fl[i], digits=3)
            llow  = "s"
          }

        if(type[i]=="BP")
          {
            Notes[i] = paste(sep=' ', type[i], "FILTER",klow,llow , "to",khigh , lhigh)
          }
        if(type[i]=="HP")
          {
            Notes[i] = paste(sep=' ', type[i], "FILTER",klow,llow )
          }
        if(type[i]=="LP")
          {
            Notes[i] = paste(sep=' ', type[i], "FILTER",khigh , lhigh)
          }
        
        
      }


    return(Notes)
  }

