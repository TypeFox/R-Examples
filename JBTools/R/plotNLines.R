plotNLines <- function(
    ##title<< Plot many lines in one plot
    x.data = matrix(1:dim(y.data)[2], ncol = dim(y.data)[2], nrow =dim(y.data)[1], byrow=TRUE)                     
                                ##<< numeric matrix/data frame: x-values with one series of
                                ##   values per row
    ,y.data                     ##<< numeric matrix/data frame: y-values with one series of
                                ##   values per row
    ,option=c('normal','diff.scales','stacked')[1] ##<< character: which type of plot to use (see details)
    ,n.lines.max=30             ##<< integer: only for 'stacked plots': how many lines to draw per panel
    ,grid=FALSE                 ##<< logical: only for 'stacked plots': whether to draw a primitive grid
    ,scale=1                    ##<< numeric: only for 'stacked plots': scale factor to scale the y scale of the stacked plots
    ,plot.scale=TRUE            ##<< TRUE: only for 'stacked plots': whether to add a small scale showing the y axis scale
    ,function.add=function()(1) ##<< function: only for 'stacked plots': function to call after plotting the individual panels
    ,bgc='white'                ##<< color: color of the plot background
    ,pch.axlink=1               ##<< integer: only for 'stacked plots': pch value for the symbol that links each line to its x-axis
    ,type='l'                   ##<< standard plotting parameter
    ,colors=c()                 ##<< colors to use for the different series
    ,xlim=c()                   ##<< standard plotting parameter
    ,ylim=c()                   ##<< standard plotting parameter
    ,xlab=c()                   ##<< standard plotting parameter
    ,ylab=c()                   ##<< standard plotting parameter
    ,labels=c()                 ##<< standard plotting parameter
    ,lty=1                      ##<< standard plotting parameter
    ,title                      ##<< standard plotting parameter
    ,yaxt='s'                   ##<< standard plotting parameter
    ,...                        ##<< further arguments passed to the plot() calls
)
##description<< plotNLines function uses different techniques to visualize many line plots on one display.
##details<< Many parameters are identical to standard plotting parameters (see ?par,?plot)
##and are not explained here.
##The function offers three options:
##'normal':
## plots all plots in one coordinate system colored according to colors
##'diff scales':
## plots all plots in the same region but uses different y axis scales. The visibility may
## be limited to ~7 plots.
##'stacked':
##plots many plots in one region, all shifted vertically a bit to increase visibility. This allows
##for the easy comparison of many similar plots (e.g. time series) but reduces the details that are
##visible.
##\if{html}{\out{<img src="../doc/plotNLines_demo.png" alt="image ..plotNLines_demo should be here"/>}}\ifelse{latex}{}{}
{
  ## Preparation
  
  n.plots <- dim(y.data)[1]

  if (is.null(dim(x.data))) {
    if (length(x.data) == dim(y.data)[2]) {
      x.data <- matrix(x.data, ncol = length(x.data), nrow = dim(y.data)[1], byrow = TRUE)
    } else {
      stop('If x.data is a vector, it has to be of the same length as dim(y.data)[2]!')
    }
  }
  
  
  if ( dim(x.data)[2] != dim(y.data)[2]) 
    stop('Error: dimensions of x and y data differ!\n')   
  
  if (length(option)==3)
    option=option[1]
  
  if (length(lty)==1)
    lty=rep(lty,times=n.plots)
  
  if (length(colors)==0) {
    if (!(option=='diff.scales')) {
      colors = rep(1,times=dim(x.data)[1])
    } else {
      colors=1:n.plots
    }
  } else if (length(colors)==1) {
    colors=rep(colors,n.plots)
  } else if (!(length(colors)<2) & !(length(colors)==n.plots)) {
    stop('Colvector has to be 1 or have the same length as dim(x.values)[1]!')
  }
  
  ## Scaled: All plots in one with different y axes
  if (option=='diff.scales') {
    if (n.plots>5)
      {
        cat('Error: Too many differnt rows in data to plot scales for each variable!')
        return()
      }
    
    old.margins <- par()$mar
    old.mgp     <- par()$mgp
    old.ps      <- par()$ps
    new.margins <- old.margins
    new.margins[4] <- 1 + (n.plots %/%2)*1.5
    new.margins[2] <- 1 + ((n.plots + 1)%/%2)*1.5
    
    if(n.plots<3) {
      ylab <- rownames(y.data)
    } else {ylab=''}
    
    par(mar=new.margins,tck=0.02,mgp=old.mgp/c(1,2,1))
    plotBG(bgc)
    plot(x.data[1,], y.data[1,], type = type, xlab = '', ylab = '', col = colors[1], lty = lty[1], yaxt = yaxt, ...)
    mtext(ylab[1], side = 2, col = colors[1], line = 1,cex = 1.2)
    mtext(xlab,side = 1,line = old.mgp[1], cex = 1.2)
    
    xlab=''
    for(i in 2:n.plots) {
      function.add()
      par(new=TRUE)
      side.axis   <- 4-((i%%2)*2)
      offset.axis <- ((i-1)%/%2)*1.5
      if (i==2 && length(ylab)>1) {
        ylab <- ylab[2]
      } else {ylab <- ''}
      plot(x.data[i,],y.data[i,],type=type,axes=FALSE,xlab='',ylab='',col=colors[i],lty=lty[i],
           yaxt=yaxt,...)
      if (yaxt=='s')
        {
          axis(side=side.axis,col=colors[i],line=offset.axis,col.axis=colors[i],cex.lab=1)
          mtext(ylab,side=side.axis,col=colors[i],line=offset.axis+1,cex=1.2)
        }
    }
  }
  
  
  
  ## Normal: All plot in one with different colors
  if (option=='normal') {
    if (length(xlim)==0)
      xlim=c(min(x.data[!is.na(y.data)],na.rm=TRUE),max(x.data[!is.na(y.data)],
                                          na.rm=TRUE))
    if (length(ylim)==0)
      ylim=c(min(y.data,na.rm=TRUE),max(y.data,na.rm=TRUE))
                                        #       plotBG(bgc)
    plot(x.data[1,],y.data[1,],col=colors[1],type=type,
         ,xlim=xlim,ylim=ylim,lty=lty[1],yaxt=yaxt,xlab='',ylab='',...)
    
    for (i in 2:dim(x.data)[1])
      points(x.data[i,],y.data[i,],col=colors[i],type=type,lty=lty[i],...)
    mtext(xlab,side=1,outer=TRUE,line=par()$mgp[1])
    mtext(ylab,side=2,outer=TRUE,las=0)
  }
  
  ## Stacked: All plots sperately, each slightly shifted vertically
  if (option =='stacked') {
    n.points      <- dim(x.data)[2]
    n.lines       <- dim(x.data)[1]
    nr.plots      <- ceiling(n.lines/n.lines.max)
    n.lines.pp    <- floor(n.lines/nr.plots)
    plots.open    <- par()$mfrow[1]*par()$mfrow[2]
    if(length(labels)==0)
      labels <- 1:n.lines
    if (plots.open > 1 && nr.plots >plots.open)
      {
        cat('Error: too many lines to plot. Open more plot areas!\n')
        return()
      }
    if (plots.open==1 && nr.plots > 1)
      par(mfrow=c(1,nr.plots),oma=c(1.5,1.5,0.1,3),mar=c(1.5,1.7,0,0),tcl=0.2,las=1,mgp=c(3,0.05,0))
    
    sd.all        <- sd(as.vector(t(y.data)),na.rm=TRUE)
    mean.all      <- mean(y.data,na.rm=TRUE)
    y.data.trans  <- ((y.data - mean.all)/sd.all)*scale

    offset <- seq(1,n.lines.pp,by=1)
    
    if (length(xlim)==0) {
      x.lim <- c(min(x.data[!is.na(y.data)],na.rm=TRUE),max(x.data[!is.na(y.data)],na.rm=TRUE))
      points.plot <- matrix(TRUE,ncol=n.points,nrow=n.lines)
    } else {
      x.lim <- xlim
      points.plot <- x.data > x.lim[1] & x.data < x.lim[2]
    }
    line.matrix   <- 1
    for (h in 1:nr.plots) {
      if (bgc!='white' )
        plotBG(bgc)
      plot(x.data[line.matrix,points.plot[line.matrix,]],
           y.data.trans[line.matrix,points.plot[line.matrix,]]+offset[1],xlim=x.lim,
           ylim=c(1-(2*scale),n.lines.pp+2*scale),type=type,col=colors[line.matrix],yaxt='n',
           ylab='',xlab='',lty=lty[1],yaxt=yaxt, ...)
      pos.zero <- which.min(abs(y.data.trans[line.matrix,points.plot[line.matrix,]]))
      if (!(length(pos.zero))==0)
        points(x.data[line.matrix,which(points.plot[line.matrix,])[pos.zero]],offset[1],
               col=colors[line.matrix],cex=1.5,pch=pch.axlink)
      if (grid)
        points(x=userCoords(c(0,1))$x,y=rep(offset[1],2),lty=2,type='l')
      line.matrix=line.matrix+1
      for (i in 2:n.lines.pp) {
        points(x.data[line.matrix,points.plot[line.matrix,]],
               y.data.trans[line.matrix,points.plot[line.matrix,]]+offset[i],type=type,
               col=colors[line.matrix],lty=lty[i],...)
        pos.zero <- which.min(abs(y.data.trans[line.matrix,points.plot[line.matrix,]]))
        if (!(length(pos.zero))==0)
          points(x.data[line.matrix,which(points.plot[line.matrix,])[pos.zero]],
                 offset[i],col=colors[line.matrix],cex=1.5,pch=pch.axlink)
        line.matrix <- line.matrix + 1
      }
      if (grid)
        for (l in 1:length(offset))
          points(x=userCoords(c(0,1))$x,y=rep(offset[l],2),lty=2,type='l')
      if (yaxt!='n')
        axis(2,at=1: n.lines.pp,labels=labels[((h-1)*n.lines.pp+1):(h*n.lines.pp)])
      function.add()
    }
    mtext(xlab,side=1,outer=TRUE,line=par()$mgp[1])
    mtext(ylab,side=2,outer=TRUE,las=0)
    if (plot.scale) {
      labels.scale = pretty(range(y.data,na.rm=TRUE),n=3)
      at           = 1+((labels.scale/sd.all))*scale
      axis(side=4,at=at,labels=labels.scale,line=.5,mgp=c(3,0.5,0))
    }
  }
  
}
