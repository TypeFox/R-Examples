seisorder<-function(GH, ORD,  VNE= c("V", "N", "E") )
  {
###  given an RSEIS list used for swig 
###   get a convenient ordering

    ##  the order depends on the input ORD which has station information
###   and VNE which is the ordering of the components
    

    if(missing(VNE)) VNE= c("V", "N", "E")
    
    oro = order(ORD$dist)

    oname = ORD$name[oro]


    if(is.null(VNE))
      {
        gtags = GH$STNS
        otags = oname
      }
    else
      {

        alltags = NULL
        for(i in 1:length(VNE))
          {

            tag1 = paste(oname, VNE[i], sep="." )

            alltags = cbind(alltags,  tag1)

          }
        otags =  as.vector( t(alltags ))
        gtags = paste(GH$STNS,GH$COMPS , sep=".")
        
      }
    
    ma = match( otags, gtags )

    ma = ma[!is.na(ma)]
    return(ma)
    

  }

