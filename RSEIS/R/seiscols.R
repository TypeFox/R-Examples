seiscols<-function(GH, acols="black", M="STNS" )
  {
    if(missing(acols))
      {

        acols = c("black", "darkmagenta", "forestgreen", "blueviolet", 
          "tan3", "lightseagreen", "deeppink", "cyan3", "bisque3", 
          "magenta1", "lightsalmon3", "darkcyan", "darkslateblue", 
          "chocolate4", "goldenrod4", "mediumseagreen")

      }
    if(missing(M)) { M = "STNS" }

    if(identical(M,"STNS") )
      {
    ustn = unique(GH$STNS)
    mcol = match(GH$STNS, ustn)

  }
    
    if(identical(M,"COMPS") )
      {
    ustn = unique(GH$COMPS)
    mcol = match(GH$COMPS, ustn)

  }


    zcol = mcol %% length(acols)
    zcol[zcol == 0 ] = length(acols)

    mcol = zcol
    
   
    

    return(acols[mcol])

  }
