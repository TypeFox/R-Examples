plot.genarch <- function(x, min.vi = .5, main=NULL, cex.axis=1, 
                         cex.names=1, cex.main=1, maxval=NULL, 
                         minval=NULL, viridis=F, ...){
  data <- x
  # lets store the graphic paratmeters so we leave people 
  # unscathed for their next plots
  old.par <- par()
  old.par <- old.par[c(1:12,14:18,20,24:53,55:72)]
  
  
  
  # lets reduce the genarch to just what we need
  data <- data[-1]
  
  # if the user is getting rid of too much data lets stop and warn them
  while(sum(as.numeric(data[[3]][,2]) >= min.vi) < 2){
    min.vi <- min.vi - .05
    #stop("select a lower min.vi to include more CGEs", call.=F)
  }
  
  # reformat and keep only the data to be plotted
  data <-cbind(t(data[[1]])[-1,],
               as.numeric(data[[3]][,2]))[as.numeric(data[[3]][,2]) >= min.vi,]
  colnames(data)[3] <- "vi"
  
  # get some extra room
  par(mar=c(2, 2, 2, 6))
  
  # make colors for barplot
  foo.colors <- heat.colors(100)[100:1]
  # add viridis
  if(viridis == T){
    foo.colors <- viridis::viridis(100)[100:1]
  }
  
  # set up graph max and min y axis if user does not provide
  if(is.null(maxval)){
    maxval <- max(as.numeric(data[, 1]) + as.numeric(data[, 2]))
    if(maxval < 0) maxval <- 0
  }
  if(is.null(minval)){
    minval <- min(as.numeric(data[, 1]) - as.numeric(data[, 2]))
    if(minval > 0) minval <- 0
  }

  # give the graph a title if the user fails to
  if(is.null(main)) main <- "Model Weighted Averages and Unconditional SE"
  
  mp <- barplot(as.numeric(data[, 1]), 
                    names.arg=row.names(data), 
                    col=foo.colors[round(100 * as.numeric(data[,3]))],
                    ylim = c(minval - .4 * abs(minval), 
                             maxval + .4 * abs(maxval)),
                    main = main, 
                    cex.axis=cex.axis, cex.names=cex.names, cex.main=cex.main)
      segments(mp, 
               as.numeric(data[, 1]) - as.numeric(data[, 2]), 
               mp, 
               as.numeric(data[, 1]) + as.numeric(data[, 2]), 
               lwd=2)
      # Now plot the horizontal bounds for the error bars
      # 1. The lower bar
      segments(mp - 0.1, as.numeric(data[, 1]) - as.numeric(data[, 2]), 
               mp + 0.1, as.numeric(data[, 1]) - as.numeric(data[, 2]), 
               lwd = 2)
      # 2. The upper bar
      segments(mp - 0.1, as.numeric(data[, 1]) + as.numeric(data[, 2]), 
               mp + 0.1, as.numeric(data[, 1]) + as.numeric(data[, 2]), 
               lwd = 2)
      # add a legend for the variable importance
      locs <- par("usr")  
      color.legend(locs[2] ,  #xl
                   locs[4]- ((locs[4]-locs[3]) * 0.5),  #yb
                   locs[2] + (locs[2]*.05),    #xr
                   locs[4],        #yt
                   legend = c("0.00", "0.25", "0.50", "0.75","1"),
                   rect.col = foo.colors,
                   cex = 0.6,
                   align="rb",
                   gradient = "y")
      par(xpd = TRUE)
      text(x = (((locs[2]) + (locs[2] + (locs[2]*.05)))/2), 
           y = locs[4],
           labels = "Variable\nImportance", cex = 0.5, pos = 3)
      # lets reset peoples graphics paramters so they make sinces again
      par(old.par)
      
}