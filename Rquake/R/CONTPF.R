CONTPF <- function(nh, g, idev=3)
  {
    


    ## print(nh$sta)

  ## print( nh$pickfile)

    
    if(length(nh$sta)<1)
      {
        print("CONTPF: No Station file")
        dev.set( g$MAINdev)
        g$zloc = list(x=NULL, y=NULL)
        g$action="donothing"
        invisible(list(global.vars=g))
      }
    if(length(nh$pickfile)<1)
      {
        print("CONTPF:  No Pickfile 1 ")
        dev.set( g$MAINdev)
        g$zloc = list(x=NULL, y=NULL)
        g$action="donothing"
        invisible(list(global.vars=g))
      }

    
    
      if( all(is.na( nh$sta)) )
      {
         print("no station file 2 ")
         invisible(list(global.vars=g))
      }
    
    if(is.null(nh$pickfile))
      {
        print("CONTPF: no pickfile....converting")
        if(length(g$WPX)<1)
          {
            print("CONTPF:  No WPX")
            dev.set( g$MAINdev)
            g$zloc = list(x=NULL, y=NULL)
            g$action="donothing"
            invisible(list(global.vars=g))
          }
        
         if(is.null(g$WPX) )invisible(list(global.vars=g))


        

       ## print("CONTPF:  trying")
        twpx = g$WPX
        
        
        
        nona = which( is.na(twpx$name) )
        
        if(length(nona)>0)
          {
            twpx = RSEIS::deleteWPX(twpx, nona)
          }

      ##   print("CONTPF:  dumping")
      ##  print(data.frame(twpx))

        
        
        A1T = Qrangedatetime(twpx)
        s1 = RSEIS::secdifL(A1T$min,  twpx)
        
            pickfile =  INITpickfile(stas=nh$sta, src=NULL, WPX=twpx)
        # print("CONTPF:  dumping")
      #  print(pickfile)
        
      }
    else
      {

        pickfile = nh$pickfile

      }
    
    
    
    
    dev.set( dev.next() )

    
    cproj = contPFarrivals(pickfile, nh$sta, proj=NULL, image=TRUE , phase="G" , add=FALSE)

    if(is.null(cproj))
      {
        cat(file="", "ERROR: No picks, try again", sep="\n")

      }
    
    dev.set( g$MAINdev)
    g$zloc = list(x=NULL, y=NULL)
    g$action="donothing"
    invisible(list(global.vars=g))
  }
