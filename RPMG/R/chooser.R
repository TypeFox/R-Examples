`chooser` <-
  function(opts=c(1, 2, 5, 10, 15, 20) , ncol=5, nsel=NA, newdev=TRUE, STAY=FALSE,
           cols="red", main="", newplot=TRUE, xlim=c(0,1), ylim=c(0,1) , just="CEN", ... )
  {
###  choosfilt()
    
    if(missing(opts))
      {
        opts = c(2, 5, 10, 15, 20)
      }


    if(is.null(opts))
      {
        opts =c(2, 5, 10, 15, 20)
        

      }

    if(missing(ncol)) {  ncol=5  }
    if(missing(nsel)) {  nsel = NA   }
    if(missing(newdev)) {  newdev=TRUE   }
    if(missing(STAY)) {  STAY=FALSE   }

    if(missing(xlim)) { xlim = c(0,1) }
    if(missing(ylim)) { ylim=c(0,1) }
    if(missing(just)) {  just="CEN"  }


    just = toupper(just)
    mjust = match(just, c("CEN", "LEFT", "RIGHT"))
    if(is.na(mjust) )
      {
        just="CEN"
      }

    if(missing(newplot)) { newplot=TRUE }
  
    if(is.na(nsel)) { nsel = length( opts) } 

    if(missing(main)) { main =  paste( sep=" ", "Choose by Clicking up to",nsel, "selections" ) }
    
    
    olddev = dev.cur()
    
    if(newdev) dev.new()

    
    N = length( opts)
    nrow = round((N/ncol)+.5)

    if(ncol==1)
      {

        nrow = N

      }

    
    dx = 1/ncol
    dy =  1/nrow

    lolab = as.character(opts) 
    lab = paste(sep='\n', lolab)

    if(length(lab)<(ncol*nrow))
      {
        lab = c(lab, rep(NA,  times = (ncol*nrow) -  length(lab))) 
      }

    if(missing(cols))
      {
        cols =pastel.colors(N, seed=1)

        cols[is.na(lab)] = NA

      }

    ##print(c(ncol, nrow))
    
    
    ##  B =  itoxyz(1:N, ncol, nrow, 1)

  if(newplot)
    {
         plot(xlim, ylim, type='n', axes=FALSE,xlab='', ylab='', main=main)
        ##  title( main=main )
    }

    
      if(FALSE)
        {

          xlim =c(.25, .75)
          ylim = c(.25, .75)

          rect(xlim[1], ylim[1], xlim[2], ylim[2]) 
          
        }

   

    dx = (xlim[2]-xlim[1]   ) /ncol
    dy =  (ylim[2]-ylim[1])   /nrow

 
    x = seq(from=xlim[1], by=dx, length=ncol)

    
    y = seq(from=ylim[1], length=nrow, by=dy)


    M =   meshgrid(x, y)

    
    
    rect(M$x, M$y, M$x+dx, M$y+dy, col=cols)

    if(just=="CEN") text(M$x+dx/2, M$y+dy/2, lab)      
    if(just=="LEFT") text(M$x, M$y, lab, adj=c(0,0) )      
    if(just=="RIGHT") text(M$x+dx, M$y+dy, lab, adj=c(1,0))      

    if(nsel==0)
      {
        if(olddev>1)dev.set(olddev)
        return(NULL)


      }

    z = locator(n=nsel, type='p', ...)

    if(length(z$x)<1)
      {
        
        if(STAY==FALSE) dev.off(dev.cur())
        return(NULL)
        
      }


    thex = z$x-xlim[1]
    they = z$y-ylim[1]

    
    ii = 1+floor(thex/dx)
    jj = 1+floor(they/dy)
    B =  jj+(ii-1)*(nrow)

    i = B

    GIVE =  opts[i]

####    print(c(i, ii, jj, GIVE))
####   print(  cbind(1:N, opts) )

    if(STAY==FALSE) dev.off(dev.cur())

    attr( GIVE,"params"  ) <- list(ind=i, dx=dx, dy=dy, nrow=nrow, ncol=ncol)
    
    return(GIVE)

    

  }

#######   choosdecim(opts=c(0, 1, 2, 5, 10, 15, 20, "None"))
#######   choosdecim(opts=floor( runif(10, 1,100) )   )
