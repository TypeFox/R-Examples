lastPIX<-function(nh, g )
  {

    if(is.null(g$LWPX))
    {
      g$zloc = list(x=NULL, y=NULL)
      g$action="replot"
      invisible(list(global.vars=g))
      
    }
    else
      {
    
    if( identical(legitWPX(g$LWPX, quiet=FALSE),0)  ) {
      
      cat("No saved Picks", sep="\n")
    }
    else
      {

        g$WPX = g$LWPX

      }

  }

    
    g$zloc = list(x=NULL, y=NULL)
    g$action="replot"
    invisible(list(global.vars=g))
    

  }
