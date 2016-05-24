`selAPX` <-
function(APX, ista=NULL , icomp=c("V", "N", "E") )
  {
    ###############   select from an APX WPX or APIX list a subset
    if(missing(icomp)) { icomp=NULL  }

 
    MPX = data.frame(APX, stringsAsFactors = FALSE)

    #####aname = APX$name
    #####acomp = APX$comp

    if(is.null(icomp))
      {
        atag = APX$name
        itag = ista
        w1 = which(!is.na(match(atag, itag)))
        m = w1
      }
    else
      {
        
        M = RPMG::meshgrid(1:length(ista), 1:length(icomp) )
        itag =  paste(sep=".",ista[M$x], icomp[M$y])
        atag = paste(sep=".", APX$name, APX$comp)
        m =  which(!is.na(match(atag, itag)))
        
      }
    
    
  
    
    SWP =  MPX[  m ,  ]

    SWP = as.list(SWP)


    
    return(SWP)
  }

