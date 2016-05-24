`setMarkup` <-
function(LABS=NULL, PROJ=NULL)
  {
    if(missing(LABS)) { LABS= NULL }
    if(missing(PROJ)) { PROJ=NULL }
    
    if(!is.null(LABS))
      {
        N = length(LABS)
        MM = list()
        for(i in 1:N)
          {
            lab = LABS[i]
            print(paste(sep=' ', i, lab) )
            print("click in figure pairs for arrows")
            L = locator(2, type="p", col='red')
            segments(L$x[1], L$y[1], L$x[2], L$y[2], col='red')
            
            apos = readline(prompt="Type in the position (0, 1,2,3,4)  ROT(T, F)  Arrow(T, F):")
            kin = unlist(strsplit(apos, " "))
            kin = kin[kin != "" ]
            pos = as.numeric(kin[1])
            rot = kin[2]
            arr = kin[3]
            angdeg = 180*atan2(L$y[2]-L$y[1],L$x[2]-L$x[1])/pi
            x1=L$x[1]
            y1=L$y[1]
            x2=L$x[2]
            y2=L$y[2]
            
            MM[[i]] = list(x1=x1, y1=y1, x2=x2, y2=y2, lab=lab, pos=pos,
                angdeg=angdeg, ROT=rot, ARR=arr, CEX=1)
            
          }
        
      }
    else
      {
        print("click in figure pairs for arrows")
        L = locator()
        N = length(L$x) / 2
        MM = list(1)
        
        for(i in 1:N)
          {
            j = (i-1)*2+1
            
            arrows(L$x[j],L$y[j] ,L$x[j+1],L$y[j+1] )
             angdeg = 180*atan2(L$y[j+1]-L$y[j],L$x[j+1]-L$x[j])/pi

          
                lab = readline(prompt="Type in the label:")
             
            print(paste(sep=' ', i, lab) )
            rot = FALSE
            apos = readline(prompt="Type in the position (0,1,2,3,4):")
            pos = as.numeric(apos)
            x1=L$x[j]
            y1=L$y[j]
            x2=L$x[j+1]
            y2=L$y[j+1]
            

             MM[[i]] = list(x1=x1, y1=y1, x2=x2, y2=y2, lab=lab, pos=pos,
                angdeg=angdeg, ROT=rot, ARR=TRUE, CEX=1)
           
            
          }
      }

    if(!is.null(PROJ))
      {

        for( i in 1:length(MM))
          {
            LL1 = XY.GLOB(MM[[i]]$x1, MM[[i]]$y1, PROJ)
            LL2 = XY.GLOB(MM[[i]]$x2, MM[[i]]$y2, PROJ)
            MM[[i]]$lat1= LL1$lat
            MM[[i]]$lon1= LL1$lon
            MM[[i]]$lat2= LL2$lat
            MM[[i]]$lon2= LL2$lon


          }



      }


    names(MM) = LABS
    
    
    return(MM)
  }

