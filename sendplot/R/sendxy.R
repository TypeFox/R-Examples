#  sendxy.R
#
# wrapper to sendPlot.R for a single, interactive scatterplot
#



sendxy <- function(plot.call,
                   x, y, 
                   xy.lbls = NA, x.lbls = NA,y.lbls=NA,
                   x.links=NA, y.links=NA,
                   xy.links=NA,asLinks=NA,
                   xlim = NA,
                   ylim = NA,
                   mai=NA,
                   mai.prc=FALSE,
                   plt.extras=NA,
                   bound.pt=FALSE,source.plot=NA,
                   paint=FALSE,img.prog = NA,
                   resize="800x1100",
                   ps.paper="letter",ps.width=8,ps.height=11,
                   fname.root="test",dir="./",header="v2",
                   up.left=c(205,131),low.right=c(633,883),
                   spot.radius=5, automap=FALSE, automap.method="mode"
                   ){

  cat("NOTE:  sendxy function is deprecated\n      Please see xy.send \n\n\n")
  

  # check plot call length -- this is wrapper for single plot
  if(length(plot.call)>1){
    cat("NOTICE: you have chosen a wrapper for a single plot.\n The first plot call in plot.calls will be used \n Additional plotting arguments, i.e. points, lines, abline, axes, should be placed in plt.extras \n ")
    plot.calls = plot.call[1]
  }
  plot.calls = plot.call
  
    
  # single plot - make matrix of ones for layout
  mat = matrix(rep(1, 170), ncol=10, nrow=17)
  
  # if margins are specified use default
  if(is.na(mai[1])){    
    mai.mat = NA
  }else{
    mai.mat = matrix(mai, ncol=4)
  }
  
  # if xlim or ylim are not specified
  # take the range of values
  if(is.na(xlim[1])) xlim = range(x, na.rm=TRUE)
  if(is.na(ylim[1])) ylim = range(y, na.rm=TRUE)
  
  # additional plot calls (i.e. axes, points, lines, abline)
  # put into correct format
  if(length(plt.extras)>1){
    plt = list()
    plt.idx = 1
    for(i in 1:length(plt.extras)){
      plt[plt.idx] = plt.extras[i]
      plt.idx = plt.idx + 1
    }
    plt.extras = list()
    plt.extras$plot1 = plt
  }
  # run sendplot using default or set arguments
  sendplot(mat = mat, plot.calls = plot.calls, mai.mat = mai.mat ,mai.prc=mai.prc, type="scatterplot",x = x, y = y, xlim = xlim, ylim = ylim, plt.extras=plt.extras, xy.lbls = xy.lbls,x.lbls = x.lbls,y.lbls=y.lbls,x.links=x.links,y.links=y.links,xy.links=xy.links,asLinks=asLinks,bound.pt=bound.pt,resize=resize, ps.paper=ps.paper,ps.width=ps.width,ps.height=ps.height,fname.root=fname.root,dir=dir, paint=paint,source.plot=source.plot, img.prog=img.prog,up.left=up.left,low.right=low.right,spot.radius=spot.radius, header=header, automap=automap, automap.method=automap.method)

  
}# end sendxy
