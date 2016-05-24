SavePF <-
function(nh, g )
  {
    
    
    if(is.null(nh$pickfile))
      {
        print("SavePF: saving raw pix")



        if( identical(legitWPX(g$WPX),0)  ) {
          
          cat("Need to click on screen, register picks (GPIX, or YPIX)", sep="\n")
        }
        else
          {
            
            twpx= g$WPX
        
            nona = which( is.na(twpx$tag) )
            
            if(length(nona)>0)
              {
                twpx = RSEIS::deleteWPX(twpx, nona)
              }
            
            PCsaveWPX(twpx, destdir="." )
          }
        
      }
    else
      {

        if(is.null( nh$sta))
          {
            print("no station file")
            invisible(list(global.vars=g))
          }
        
    
        if( all(is.na( nh$sta)) )
          {
            print("no station file")
            invisible(list(global.vars=g))
          }
    
        
        
        
        PFoutput(nh$pickfile, nh$sta,  sol=NULL, format=c(1,2) )
      }
    
    
    
    g$zloc = list(x=NULL, y=NULL)
    g$action="donothing"
    invisible(list(global.vars=g))
  }
