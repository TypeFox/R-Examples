matsquiggle<-function(XMAT, dt1, dist=NULL, thick=1 , FLIP=FALSE, filcol='blue', tracecol="black", add=FALSE, PLOT=TRUE)
  {
#####  plot a seismic section as in exploration seismology
####  with varian - squiggle display

    ##  source("varsquig.R")

    d = dim(XMAT)
    x = seq(from=0, length=d[1], by=dt1)
    
    
    if(missing(dist))
      {
        y1 = seq(from=0, by=1, length=d[2])
       
      }
    else
      {

        y1 = dist


      }
    
   ###  print(y1)

    Ry  = range(y1)
    
   ###  fatness = (y1[2]-y1[1])/2

    if(missing(thick))
      {
        fatness = diff(Ry)/length(y1)
      }
    else
      {
        fatness = thick
      }

    if(is.null(fatness))
      {
        fatness = diff(Ry)/length(y1)

      }


    if(missing(add)) { add= FALSE }
    if(missing(PLOT)) { PLOT=TRUE }

    


    ###  if add = FALSE do not start a new plot

    if(!add){

      
    if(FLIP)
      {
        
       
        plot(range(y1), range(x), type='n', axes=FALSE, xlab="", ylab=""  )
      }
    else
      {
         
        plot(range(x),range(y1) , type='n', axes=FALSE, xlab="", ylab=""  )
      }

  }


    if(missing(filcol)) { filcol=rep('blue', length(y1))  }
    if(missing(tracecol)) { tracecol=rep('blue', length(y1)) } 
    
    
    
    if(length(filcol)==1) { filcol=rep(filcol , length(y1))  }
    if(length(tracecol)==1) { tracecol=rep(tracecol , length(y1))  }
    
    
    if(PLOT)
  {
    for(i in 1:length(y1))
      {
        if(FLIP)
          {
             L = list(y=range(x)   , x =c(y1[i]-fatness, y1[i]+fatness  ) )
          }
        else
          {
            L = list(x=range(x)   , y =c(y1[i]-fatness, y1[i]+fatness  ) )

          }
        
        varsquig(x, XMAT[,i], L=L , FLIP=FLIP, filcol=filcol[i], tracecol=tracecol[i], var=TRUE)
      }
  }


    invisible(list(thick=fatness, rx = range(x), ry=range(y1) ) )
  }
