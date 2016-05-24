`wlet.do` <-
function(why, dt, noctave=6, nvoice=20, flip=TRUE, ploty=TRUE, zscale=1, col=terrain.colors(100), STAMP=STAMP, units="", scaleloc=c(0.4,0.95) )
  {
    ### usage: wlet.do(x,  delta, noctave = 4, nvoice = 20,  flip=FALSE,  col=rainbow(100))

   

    if(missing(col)) { col=rainbow(100) }
    if(missing(noctave)) { noctave=6 }
    if(missing(nvoice)) { nvoice=20 }
    if(missing(flip)) { flip=TRUE }
    if(missing(ploty)) { ploty=TRUE }
    if(missing(zscale)) {  zscale=1  }
    if(missing(STAMP)) { STAMP=NULL }

   opar = par()
    ## par(mfrow=c(2,1))
    ## par(xaxs='i')

    ## plot.ts(ts(why, deltat=dt) )

    ####  this used to be a problem but I think they fixed
    ####  it.  If not, go back to jmlcwt
    kaha= cwt(why, noctave, nvoice=nvoice, w0=5, twoD=TRUE, plot=FALSE)
 ##   kaha= cwt(why, noctave, nvoice=nvoice, w0=5, twoD=TRUE, plot=TRUE)

###  get the scale for the y-axis
    

###    ii = sort(rep(c(1:noctave), times=nvoice))
###    jj = rep(c(0:(nvoice-1)), times=noctave)
###    sa = 2^(ii+jj/nvoice)
###   take the log
###    lsa = log2(sa)

    
    if(flip==TRUE)
      {
        baha =  mirror.matrix(Mod(kaha))
        
      }
    else
      {
        baha = Mod(kaha)
       
      }

    baha = list(img=baha, noctave=noctave , nvoice=nvoice, flip=flip)
    
    ##  wlet.plot(baha, why, dt, col=col, zscale=zscale)


    PE = plotwlet(baha, why, dt , zscale=zscale,  col=col,  ygrid=FALSE, STAMP=STAMP, units=units, scaleloc=scaleloc)

    
 
    
    invisible(list(baha=baha, PE=PE))
    
 }

