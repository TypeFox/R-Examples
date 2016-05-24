SaveCSV <-
  function(nh, g )
  {

    ##   destdir = "."
    
    twpx= g$WPX

  ##   print("SaveCSV test")
  ##   print(twpx$sec)
  ##   print(data.frame(twpx))

    if( identical(legitWPX(twpx, quiet=FALSE),0)  ) {

      cat("Need to click on screen, register picks (GPIX, or YPIX)", sep="\n")
    }
    else
      {


        ##  print("SaveCSV it is going here")

        
        nona = which( is.na(twpx$tag) )
        
        if(length(nona)>0)
          {
            twpx = RSEIS::deleteWPX(twpx, nona)
          }
        if(length(twpx$tag)<1 )
          {
            
            g$action="donothing"
            invisible(list(global.vars=g))
          }
        
        RDATES = Qrangedatetime(twpx)
        
        fout1 = PCfiledatetime(RDATES$min, 0)
        
        fout2 = paste(fout1,"csv", sep="." )
        
        ##      fout3 = paste(destdir, fout2, sep="/")

        write.csv(twpx, file=fout2)
        g$LWPX = twpx
        g$WPX = RSEIS::cleanWPX()
        
      }
    g$zloc = list(x=NULL, y=NULL)
    g$action="donothing"
    invisible(list(global.vars=g))
  }
