BUTREPLOT<- function(opts , ncol=5, HOZ=TRUE, TOP=TRUE, cols="white", main="", xlim=c(0,1), ylim=c(0,1),  newplot=TRUE)
  {
##    source("/home/lees/NEWseis/BUTREPLOT.R")

    ##    BUTREPLOT(stdlab)
    
    ##    BUTREPLOT(stdlab, cols=rep(grey(.9), length(stdlab)) )

    ##    BUTREPLOT(stdlab, cols=rainbow( length(stdlab)) )

##    BUTREPLOT(STDLAB, HOZ=FALSE, TOP=TRUE);
    
##    BUTREPLOT(STDLAB, HOZ=TRUE, TOP=TRUE);

##    BUTREPLOT(STDLAB, HOZ=TRUE, TOP=FALSE);

##    BUTREPLOT(STDLAB, HOZ=FALSE, TOP=FALSE);


    
##    BUTREPLOT(STDLAB, newplot=TRUE); 


    if(missing(ncol)) {  ncol=5  }
    if(missing(HOZ)) {   HOZ=TRUE   }
    if(missing(TOP)) {    TOP=TRUE  }

    
    if(missing(main)) { main="" }
    if(missing(newplot)) { newplot=TRUE }
    if(missing(xlim)) { xlim = c(0,1) }
    if(missing(ylim)) { ylim=c(0,1) }

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
        cols =RPMG::pastel.colors(N, seed=1)

        cols[is.na(lab)] = NA
      }


    cols[is.na(lab)] = NA

    ##print(c(ncol, nrow))
    
    ##  B =  RPMG::itoxyz(1:N, ncol, nrow, 1)

    if(newplot)
      {
        plot(xlim, ylim, type='n', axes=FALSE,xlab='', ylab='', main=main)
      }

    dx = (xlim[2]-xlim[1]   ) /ncol
    dy =  (ylim[2]-ylim[1])   /nrow

    x = seq(from=xlim[1], by=dx, length=ncol)

    y = seq(from=ylim[1], length=nrow, by=dy)

    M =   RPMG::meshgrid(x, y)

  ######  ind = seq(from=1:length(lab))

      

    if(TOP==TRUE){
      M$y = M$y[rev(1:nrow),  ] 
    }
    
    if(HOZ==TRUE){
      M$x = t(M$x)
      M$y = t(M$y)
    }

   
    rect(M$x, M$y, M$x+dx, M$y+dy, col=cols)
    text(M$x+dx/2, M$y+dy/2, lab)
    
    invisible(list(M=M, dx=dx, dy=dy, rx=range(c(x, x+dx))  , ry=range(c(y, y+dy))  ))

  }

