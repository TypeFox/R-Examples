RQ <- function(nh, g, idev=3)
  {   ####  relocation button for RSEIS::swig
   
    ## 

    if(is.null(nh$sta))
      {
        print("no station file")
        invisible(list(global.vars=g))
      }

    if(is.null(nh$vel))
      {
        
        fuj1.vel = defaultVEL(1)
        ##  data(fuj1.vel)
        nh$vel = fuj1.vel
      }

    if(is.null(nh$pickfile))
      {
        print("RQ: no pickfile....converting")
        if(is.null(g$WPX) )invisible(list(global.vars=g))
        twpx = g$WPX
        nona = is.na(twpx$tag)
        
        if(any(nona))
          {
            
            twpx =  RSEIS::deleteWPX(twpx, which(nona))
          }

        uphase = unique( twpx$phase )
        
        if( any( "G" %in% uphase ) )
          {
            
            twpx = Y2Pphase(twpx, "G" )

          }

        if( any( "Y" %in% uphase ) )
          {
            
            twpx = Y2Pphase(twpx, "Y" )

          }

        

        
        A1T = Qrangedatetime(twpx)
        s1 = RSEIS::secdifL(A1T$min,  twpx)

        nh$pickfile =  INITpickfile(stas=nh$sta, src=NULL, WPX=twpx)

      }
    
    dev.set(dev.next() )
    
    eqsol = NLSlocate(nh, vel=nh$vel,  PLOT=TRUE )
    dev.set( g$MAINdev)

    g$zloc = list(x=NULL, y=NULL)
    g$action="donothing"
    invisible(list(global.vars=g))
    
  }
