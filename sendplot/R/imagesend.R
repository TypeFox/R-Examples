#
# wrapper for single heatmap/image
#



imagesend <- function(plot.call, 
                    x.pos,
                    y.pos,
                    xy.type,
                    plot.extras = NA,
                    mai.mat=NA, mai.prc=FALSE,
                    xy.labels=NA,
                    image.size="800x1100",
                    spot.radius = 5,
                    fname.root="Splot",
                    dir="./",
                    window.size = "800x1100", # in px
                    ... # arguments in makeImap other than x.pos, y.pos, spot.radius,xy.type,  xy.labels, dir, fname.root, returnVl
                    ){


    # check plot call length -- this is wrapper for single plot
  if(length(plot.call)>1){
    cat("NOTICE: you have chosen a wrapper for a single plot.\n The first plot call in plot.calls will be used \n Additional plotting arguments, i.e. points, lines, abline, axes, should be placed in plt.extras \n ")
    plot.calls = plot.call[1]
  }
  plot.calls = plot.call
  
  # single plot - make matrix of ones for layout
  mat = matrix(rep(1, 170), ncol=10, nrow=17)
  
  Splot = initSplot(mat=mat, plot.calls= plot.calls, Iflag=TRUE, figTypes="image", plot.extras=plot.extras,source.plot="png", mai.mat=mai.mat,mai.prc=mai.prc,  image.size=image.size, returnVl=TRUE, saveFlag=FALSE)

  Splot = makeImap(Splot, figure=1, xy.type=xy.type, x.pos=x.pos, y.pos=y.pos, spot.radius=spot.radius, xy.labels=xy.labels, fname.root=fname.root, dir=dir, returnVl=TRUE, ...)

  Splot = makeSplot(Splot, fname.root=fname.root, dir=dir, window.size=window.size, returnObj=TRUE)
  

}
