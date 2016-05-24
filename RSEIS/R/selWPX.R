`selWPX` <-
  function(APX, ind=NULL,  ista=NULL , icomp=c("V", "N", "E") )
  {
###############   select from WPX
    if(missing(icomp)) { icomp=NULL  }

    
    MPX = data.frame(APX, stringsAsFactors = FALSE)

#####aname = APX$name
#####acomp = APX$comp



    if(!missing(ind))
      {

        getind =  ind


      }
    else
      {

        

        if(is.null(icomp))
          {
            atag = APX$name
            itag = ista
            w1 = which(!is.na(match(atag, itag)))
            getind = w1
          }
        else
          {

            getind  = which(ista %in% APX$name && icomp  %in% APX$comp )
            
          }

      }
    
    
    
    
    SWP =  MPX[ getind  ,  ]

    SWP = as.list(SWP)


    
    return(SWP)
  }

