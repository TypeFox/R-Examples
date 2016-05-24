`comp.env` <-
function(ex, Y, PLOT=TRUE, stamps=stamps)
  {
   #######  Takes in an common x predictor and
    ####  compares the envelopes of each column in the Y matrix
    ####  all the Y's must have the same length as ex, of course
    if(missing(PLOT)) { PLOT=TRUE }
    if(missing(stamps)) { stamps=NULL }

    
    nx = length(ex)
    ny = dim(Y)
    
    if(nx != ny[1] )
      {
        print("error, need nx=ny")
        return(0)
      }

    zw = list()

    for(i in 1:ny[2])
      {
        ave = mean(Y[,i])
        s1 = Y[,i]-ave
        
        w1 = wiggle.env(ex, s1)
        w1$ave = ave
        zw[[i]] = w1
        if(i==1) { trange = range(w1$y) } else {  trange = range(trange, range(w1$y))   }
        
      }

    if(PLOT==TRUE)
      {
        plot(range(ex), trange, type='n', xlab="s", ylab="Envelope")
        for(i in 1:length(zw))
          {
            lines(zw[[i]], col=i)

          }

        if(!is.null(stamps))
            {
              legend("topright", lty=1, col=1:ny[2], legend=stamps)
            }
   
        
        

      }
    invisible(zw)

  }

