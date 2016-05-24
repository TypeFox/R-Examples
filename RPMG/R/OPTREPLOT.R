OPTREPLOT<- function(opts , ncol=5, sel=1,  HOZ=TRUE, TOP=TRUE, cols="white", scol="black", bcol="white" , tcol="black", slwd=1, blwd=3, main="", xlim=c(0,1), ylim=c(0,1), cex=1,  mpct = 0.1 , newplot=TRUE)
  {
##    source("/home/lees/NEWseis/OPTREPLOT.R")

    ##    OPTREPLOT(stdlab)
    
    ##    OPTREPLOT(stdlab, cols=rep(grey(.9), length(stdlab)) )

    ##    OPTREPLOT(stdlab, cols=rainbow( length(stdlab)) )

##    OPTREPLOT(STDLAB, HOZ=FALSE, TOP=TRUE);
    
##    OPTREPLOT(STDLAB, HOZ=TRUE, TOP=TRUE);

##    OPTREPLOT(STDLAB, HOZ=TRUE, TOP=FALSE);

##    OPTREPLOT(STDLAB, HOZ=FALSE, TOP=FALSE);


    
##    OPTREPLOT(STDLAB, newplot=TRUE); 


    if(missing(ncol)) {  ncol=5  }
    if(missing(HOZ)) {   HOZ=TRUE   }
    if(missing(TOP)) {    TOP=TRUE  }
    if(missing(sel)) {    sel=1:length(opts)  }

    
    if(missing(main)) { main="" }
    if(missing(newplot)) { newplot=TRUE }
    if(missing(xlim)) { xlim = c(0,1) }
    if(missing(ylim)) { ylim=c(0,1) }
    if(missing(cex)) { cex=1  }
    if(missing(slwd)) { slwd=1  }
    if(missing(blwd)) { blwd=3  }

    if(missing(bcol))  { bcol = "white" }
    if(missing(tcol))  { tcol = "black" }

    N = length( opts)
    nrow = round((N/ncol)+.5)

    if(ncol==1)
      {

        nrow = N

      }

    if(missing(mpct)) mpct = 0.1
    
  
    
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


    cols[is.na(lab)] = NA

    ##print(c(ncol, nrow))
    
    ##  B =  itoxyz(1:N, ncol, nrow, 1)

    if(newplot)
      {
        plot(xlim, ylim, type='n', axes=FALSE,xlab='', ylab='', main=main)
      }

    dx = (xlim[2]-xlim[1]   ) /ncol
    dy =  (ylim[2]-ylim[1])   /nrow

    x = seq(from=xlim[1], by=dx, length=ncol)

    y = seq(from=ylim[1], length=nrow, by=dy)

    M =   meshgrid(x, y)

  ######  ind = seq(from=1:length(lab))

    mdx = dx*mpct
    mdy = dy*mpct


    if(TOP==TRUE){
      M$y = M$y[rev(1:nrow),  ] 
    }
    
    if(HOZ==TRUE){
      M$x = t(M$x)
      M$y = t(M$y)
    }

   
    rect(M$x+mdx, M$y+mdy, M$x+dx-mdx, M$y+dy-mdy, border=bcol, col=cols, lwd=slwd)
    text(M$x+dx/2, M$y+dy/2, lab, col=tcol,  cex=cex)

    if(length(sel)>0)
      {
        rect(M$x[sel]+mdx, M$y[sel]+mdy, M$x[sel]+dx-mdx, M$y[sel]+dy-mdy, border=scol, lwd=blwd )
      }


    
    invisible(list(M=M, dx=dx, dy=dy, rx=range(c(x, x+dx))  , ry=range(c(y, y+dy))  ))

  }

