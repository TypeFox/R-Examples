MixturePlot = function(x=NULL,y=NULL,z=NULL,w=NULL,des=NULL,
                       res=400,lims=c(rep(0,6)),color.palette = heat.colors,
                       constrts=FALSE,contrs=TRUE,n.breaks=10,levels=NULL,
                       cols=FALSE, despts=TRUE, mod=NA,x3lab="Fraction X3",
                       x2lab="Fraction X2", x1lab="Fraction X1", 
                       corner.labs = NULL,
                       colorkey=list(dx=0.04,x0=0.95,y0=0.45,y1=0.90,add=TRUE,mode="all"),
                       pseudo=FALSE)
{
  
  ##############################################
  #Argument list
  # des design matrix with known points (data frame)
  # x, y, z locations for known points
  # w values at x,y,z locations
  # res number of color blocks between 0 and 1 of x
  # lims vector of lower and upper constraints for x1,x2,x3
  # constrts if TRUE constraints found in lims will be added to the graph
  # contrs if TRUE contour lines will be added to the graph
  # n.breaks number of breaks between levels (used for contours if 'levels' not specified)
  # levels takes a list of contour levels (e.g. levels=c(1,3,5,10) will draw contours at those heights)
  # cols if TRUE regions between contour lines will be colored
  # despts if TRUE plots the design points in data frame des or vectors x, y, z
  # color.palette is the color palette to use.
  # mod is an indicator for the model 1=linear, 2=quadratic, or 4=special cubic
  # x3lab label for the x3 axis
  # x2lab label for the x2 axis
  # x1lab label for the x1 axis
  # corner.labs labels for x3, x2 and x1 vertices
  # colorkey  list with locations of the color key
  # psuedo if TRUE uses psuedo components to zoom in on constrained region.  Will create the smallest
  #  equilateral triangle that still contains the whole constrained region.
  
  
  if((constrts | pseudo) & sum(lims == rep(0,6)) == 6 ){
    stop("Component limits must be specified with the 'lims' argument to plot constraints or to use pseudo components")
  }
  # check for valid mod
  if(mod < 1 | mod == 3| mod > 4) 
    stop("mod must be one of the following: 1 = linear, 2 = quadratic, or 4 = Special Cubic")
  
  if (is.null(des)){
    if (is.null(x)) 
      stop("There must be a data frame containing the design, or vectors of known points")
  }  else {
    x<-des$x3
    y<-des$x2
    z<-des$x1
    w<-des$y
  }
  
  base = high = NULL 
  #Depends on the following libraries
  #library(lattice)
  #library(grid)
  ######################
  ## Creation of Grid ##
  ######################
  trian <- expand.grid(base=seq(0,1,l=res), high=seq(0,sin(pi/3), l=res))#87
  trian <- subset(trian, (base*sin(pi/3)*2)>high)
  trian <- subset(trian, ((1-base)*sin(pi/3)*2)>high)
  new2 <- data.frame(z=trian$high*2/sqrt(3))
  new2$x <- trian$base-trian$high/sqrt(3)
  new2$y <- 1 - new2$x-new2$z
  
  
  if(pseudo){
    l.bnds <- lims[seq(1,5,by=2)]
    sum.bnds <- sum(l.bnds)
    new2$x <- l.bnds[3]+(1-sum.bnds)*new2$x
    new2$y <- l.bnds[2]+(1-sum.bnds)*new2$y
    new2$z <- l.bnds[1]+(1-sum.bnds)*new2$z
    x.pseudo <- (x-l.bnds[3])/(1-sum.bnds)
    y.pseudo <- (y-l.bnds[2])/(1-sum.bnds)
    z.pseudo <- (z-l.bnds[1])/(1-sum.bnds)
  }
  
  if(is.na(mod)==FALSE){
    ## Fit the data to a model
    if (mod==1) {
      ## This is the Scheffe Linear model
      fm1 = lm(w~x+y+z-1)
    } 
    if (mod==2) {
      ## This is the Scheffe Quadratic model
      fm1 = lm(w~x+y+z+x*y+x*z+y*z-1)
    } 
    if (mod==4) {
      ## This is the Scheffe Special Cubic Model
      fm1 = lm(w~x+y+z+x*y+x*z+y*z+x*y*z-1)
    }
  }
  ## Create a new dataset using the model
  if (is.na(mod)==TRUE){
      stop("There must be a model specified")
    }
    else{	
    trian$w = predict(fm1, newdata=data.frame(new2))
  }	
  
  ## Function for laying out barycentric coordinates
  grade.trellis <- function(from=0.2, to=0.8, step=0.2, col=1, lty=3, lwd=.5){
    if (constrts) {
      #Constraints on x1
      f1<-lims[1]
      t1<-lims[2]
      s1=t1-f1
      x1 <- seq(f1, t1, s1)
      x2 <- x1/2
      y2 <- x1*sqrt(3)/2
      x3 <- (1-x1)*0.5+x1
      y3 <- sqrt(3)/2-x1*sqrt(3)/2
      panel.segments(x2, y2, 1-x2, y2, col=col, lty=2, lwd=2.0)
      #Constraints on x2 (note backwards f2-1-upper, t2=1-lower
      f2<-1-lims[4]
      t2<-1-lims[3]
      s2<-t2-f2
      x1 <- seq(f2, t2, s2)
      x2 <- x1/2
      y2 <- x1*sqrt(3)/2
      x3 <- (1-x1)*0.5+x1
      y3 <- sqrt(3)/2-x1*sqrt(3)/2
      panel.segments(x1, 0, x2, y2, col=col, lty=2, lwd=2.0)
      #Constraints on x3
      f3<-lims[5]
      t3<-lims[6]
      s3<-t3-f3
      x1 <- seq(f3, t3, s3)
      x2 <- x1/2
      y2 <- x1*sqrt(3)/2
      x3 <- (1-x1)*0.5+x1
      y3 <- sqrt(3)/2-x1*sqrt(3)/2
      panel.segments(x1, 0, x3, y3, col=col, lty=2, lwd=2.0)
    }
    #Grid lines
    x1 <- seq(from, to, step)
    x2 <- x1/2
    y2 <- x1*sqrt(3)/2
    x3 <- (1-x1)*0.5+x1
    y3 <- sqrt(3)/2-x1*sqrt(3)/2
    panel.segments(x1, 0, x2, y2, col="darkgrey", lty=lty, lwd=lwd)
    panel.text(x1, 0, label=x1, pos=1)
    panel.segments(x1, 0, x3, y3, col="darkgrey", lty=lty, lwd=lwd)
    panel.text(x2, y2, label=rev(x1), pos=2)
    panel.segments(x2, y2, 1-x2, y2, col="darkgrey", lty=lty, lwd=lwd)
    panel.text(x3, y3, label=rev(x1), pos=4)
    
  }
  
  grade.trellis.pseudo <- function(from=0.2, to=0.8, step=0.2, col=1, lty=3, lwd=.5){
    #Constraints on x1
    x1 <- (lims[2]-l.bnds[1])/(1-sum.bnds)
    x2 <- x1/2
    y2 <- x1*sqrt(3)/2
    x3 <- (1-x1)*0.5+x1
    y3 <- sqrt(3)/2-x1*sqrt(3)/2
    panel.segments(x2, y2, 1-x2, y2, col=col, lty=2, lwd=2.0)
    #Constraints on x2 (note backwards f2-1-upper, t2=1-lower
    x1 <- 1-(lims[4]-l.bnds[2])/(1-sum.bnds)
    x2 <- x1/2
    y2 <- x1*sqrt(3)/2
    x3 <- (1-x1)*0.5+x1
    y3 <- sqrt(3)/2-x1*sqrt(3)/2
    panel.segments(x1, 0, x2, y2, col=col, lty=2, lwd=2.0)
    #Constraints on x3
    f3<-lims[5]
    t3<-lims[6]
    s3<-t3-f3
    x1 <- (lims[6]-l.bnds[3])/(1-sum.bnds)
    x2 <- x1/2
    y2 <- x1*sqrt(3)/2
    x3 <- (1-x1)*0.5+x1
    y3 <- sqrt(3)/2-x1*sqrt(3)/2
    panel.segments(x1, 0, x3, y3, col=col, lty=2, lwd=2.0)
    #Grid lines
    x1 <- seq(from, to, step)
    x2 <- x1/2
    y2 <- x1*sqrt(3)/2
    x3 <- (1-x1)*0.5+x1
    y3 <- sqrt(3)/2-x1*sqrt(3)/2
    labx1 <- l.bnds[3]+(1-sum.bnds)*x1
    labx2 <- l.bnds[2]+(1-sum.bnds)*x1
    labx3 <- l.bnds[1]+(1-sum.bnds)*x1
    panel.segments(x1, 0, x2, y2, col="darkgrey", lty=lty, lwd=lwd)
    panel.text(x1, 0, label=labx1, pos=1)
    panel.segments(x1, 0, x3, y3, col="darkgrey", lty=lty, lwd=lwd)
    panel.text(x2, y2, label=rev(labx2), pos=2)
    panel.segments(x2, y2, 1-x2, y2, col="darkgrey", lty=lty, lwd=lwd)
    panel.text(x3, y3, label=rev(labx3), pos=4)
    
  }
  
  ## Perform the actual plotting
  if(is.null(levels)){
    p <- levelplot(w~base*high, trian, aspect="iso", xlim=c(-0.1,1.1), ylim=c(-0.1,0.96),
                   x3lab=NULL, x2lab=NULL, contour=contrs, cuts=n.breaks, labels=TRUE, pretty=TRUE, region=cols, 
                   col.regions = color.palette(n=n.breaks+1), cex.lab=1.3, 
                   par.settings=list(axis.line=list(col=NA), axis.text=list(col=NA)),
                   panel=function(..., at=pretty(trian$w,n=11), contour=TRUE, labels=pretty(trian$w,n=11)){
                     panel.levelplot(..., at=pretty(trian$w,n=11), contour=TRUE, labels=pretty(trian$w,n=11),
                                     lty=2, lwd=0.5, col=1)}
    )
  }else{
    p <- levelplot(w~base*high, trian, aspect="iso", at=levels, xlim=c(-0.1,1.1), ylim=c(-0.1,0.96),
                   x3lab=NULL, x2lab=NULL, contour=contrs, labels=TRUE, pretty=TRUE, region=cols, 
                   col.regions = color.palette(n=n.breaks+1), cex.lab=1.3, 
                   par.settings=list(axis.line=list(col=NA), axis.text=list(col=NA)),
                   panel=function(..., at=pretty(trian$w,n=11), contour=TRUE, labels=pretty(trian$w,n=11)){
                     panel.levelplot(..., at=pretty(trian$w,n=11), contour=TRUE, labels=pretty(trian$w,n=11),
                                     lty=2, lwd=0.5, col=1)}
    )
  }
  ##labels and legend
  grid.newpage()
  pushViewport(viewport(xscale = p$x.limits, yscale = p$y.limits))
  do.call(panel.levelplot, trellis.panelArgs(p, 1))          
  ## update the trellis panel
  #trellis.focus("panel", 1, 1, highlight=TRUE)
  panel.segments(c(0,0,0.5), c(0,0,sqrt(3)/2), c(1,1/2,1), c(0,sqrt(3)/2,0),lwd=2)
  if(pseudo){
    grade.trellis.pseudo()	
  }else{
    grade.trellis()
  }
  panel.text(0, 0, label=corner.labs[2], pos=2)
  panel.text(1/2, sqrt(3)/2, label=corner.labs[3], pos=3)
  # This point is x1 vertex
  #panel.points(1/2,sqrt(3)/2,col="black",cex=1.4,pch=19)
  panel.text(1, 0, label=corner.labs[1], pos=4)
  panel.text(.5,-.075,x3lab)
  panel.text(.18,.5,x2lab,srt=60)
  panel.text(.82,.5,x1lab, srt=-60)
  
  if (despts) {
    if(pseudo){
      # Plot the design points
      # using the transformation x=z*sqrt(3)/4+x+.065*z, y=z*sqrt(3)/2
      xpts<-(z.pseudo*sqrt(3)/4)+.065*z.pseudo+x.pseudo
      ypts<-z.pseudo*sqrt(3)/2
      panel.points(xpts,ypts,pch=19,cex=1.4,col="black")
    }else{
      xpts<-(z*sqrt(3)/4)+.065*z+x
      ypts<-z*sqrt(3)/2
      panel.points(xpts,ypts,pch=19,cex=1.4,col="black")
    }
  }
  
  ck.x=colorkey$x0
  ck.y.b=colorkey$y0 #.45
  ck.y.t=colorkey$y1 #.90
  ck.y = seq(ck.y.b,ck.y.t,len=n.breaks+2)
  d.x = colorkey$dx
  d.y = diff(ck.y[1:2])
}
#Function over
