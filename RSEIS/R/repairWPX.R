repairWPX<-function(wpx)
  {
###    returns 0 for failure, 1 for success
####  return an empty WPX
    zpx = cleanWPX()
    nam1 = names(zpx)
    nam2 = names(wpx)

    nn = as.vector( unlist( lapply(wpx, "length") ) )

    k = which(nam2=="sec")

    Npix  =nn[k]
    
###  check to make sure all the elements of a basic wpx list are there	
    m1 = match(nam1, nam2)

    ww = which(is.na(m1))
    
    for(i in 1:length(m1))	
      {
        if(!is.na(m1[i]))
          {
            zpx[[i]] = wpx[[m1[i] ]]
          }
        else
          {
            zpx[[i]] = rep(NA, times=Npix)
          }
      }

    #######  fix dates?

    
 nn = as.vector( unlist( lapply(zpx, "length") ) )
    Npix  = max(nn)
   wn =  which(nn==1)

    if(length(wn)>0)
      {
    for(i in 1:length(wn))
      {
        zed = zpx[[ wn[i]  ]]
        zpx[[ wn[i]  ]] = rep(zed[1], times=Npix)
      }
  }

    nn = as.vector( unlist( lapply(zpx, "length") ) )

    if(length(unique(nn))>1)
      {
        print("repairWPX: There are still problem with this WPX")
      }
    #######   fix names?


    
    ##   WPX = data.frame(WPX, stringsAsFactors = FALSE)
    invisible(zpx)
  }
