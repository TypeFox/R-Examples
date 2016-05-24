# File plot_2parOF.R
# Part of the hydroPSO R package, http://www.rforge.net/hydroPSO/ ; 
#                                 http://cran.r-project.org/web/packages/hydroPSO
# Copyright 2010-2012 Mauricio Zambrano-Bigiarini & Rodrigo Rojas
# Distributed under GPL 2 or later

################################################################################
#                               'plot_2parOF'                                  #
################################################################################
# Author : Mauricio Zambrano Bigarini                                          #
# Started: Nov 30th, 2010                                                      #      
# Updates: 17-Jan-2010 ; 25-Jan-2011                                           #
#          15-Feb-2012 ; 08-Mar-2012 ; 23-Mar-2012 ; 20-Nov-2012               #
################################################################################
# Purpose: For two user-defined parameters, it plots the values of the         #
#          objective funtion in a two dimensional box, where the boundaries    # 
#          of each parameter are used as axis.                                 #
################################################################################

# params : matrix or data.frame with the parameter sets used during calibration
# p1.name: character with the name of the 1st parameter to be plotted
# p2.name: character with the name of the 2nd parameter to be plotted
# type   : character, indicating the type of plot. Valid values are: 
#          -) "sp"       : spatial plot
#          -) "scatter3d": 3d scatterogram
# auto.key : logical, indicating if the legend has to be drawn or not
# key.space: character,indicating the position of the legend with respect to the plot
#
plot_2parOF <- function(params, 
                        gofs,
                        p1.name, 
                        p2.name, 
                        type="sp", 
                        MinMax=c("min", "max"),
                        gof.name="GoF", 
                        main=paste(gof.name, "Surface"),
                        GOFcuts,
                        colorRamp= colorRampPalette(c("darkred", "red", "orange", "yellow", "green", "darkgreen", "cyan")),
                        points.cex=0.7, 
                        alpha=0.65,
                        axis.rot=c(0, 0),
                        auto.key=TRUE, 
                        key.space= "right"
                        ) {

    # Checking 'params'
    if (missing(params)) 
      stop("Missing argument: 'params' must be provided !!" )

    # Number of parameter sets
    n <- nrow(params)
    
    # Checking 'gofs'
    if (missing(gofs)) {
      stop("Missing argument: 'gofs' must be provided !!" )
    } else if (length(gofs) != n)
        stop("Invalid argument: 'length(gofs) != nrow(params)' (", length(gofs), "!=", n, ") !!" )    
        
    # Setting 'MinMax' 
    MinMax <- match.arg(MinMax)    

    # If the user provided 'p1.name', it checks that the field 'p1.name' exists in 'params'
    if (!missing(p1.name)) {
      if ( !(p1.name %in% colnames(params)) )
      stop("Invalid argument: The field '", p1.name, "' doesn't exist in 'params'")
    } # IF end

    # If the user provided 'p2.name', it checks that the field 'p2.name' exists in 'params'
    if (!missing(p2.name)) {
      if ( !(p2.name %in% colnames(params)) )
      stop("Invalid argument: The field '", p2.name, "' doesn't exist in 'params'")
    } # IF end

    # Checking the value of 'type'
    if (is.na(match(type, c("sp", "scatter3d") ) ) ) {
     stop( "Invalid argument: 'type' must be in c('sp', 'scatter3d')" )
    } else if (type=="scatter3d") {  
            if ( is.na( match("scatter3d", installed.packages()[,"Package"] ) ) ) {
               warning("Package 'scatter3d' is not installed =>  type='sp'")
               type <- "sp"
            } # IF end 
         } # IF end

    # Selecting only the 2 parameters provided by the user + the objective function
    p              <- data.frame(params[ ,c(p1.name, p2.name)], GoF= gofs)
    colnames(p)[3] <- gof.name
    
    # Ordering parameter sets to have the best ones over the rest
    #ifelse(MinMax=="max", decreasing<-FALSE, decreasing<-TRUE)
    decreasing <- TRUE
    if (MinMax=="max") decreasing <- FALSE
    p <- p[order(p[, gof.name], decreasing = decreasing), ]
    
    if (type=="sp") {
      colnames(p) <- c("x", "y", gof.name)
      
      # If the user didn't provide 'GOFcuts', the 5 quantiles are used
      if (missing(GOFcuts)) 
         GOFcuts <- unique(fivenum(as.numeric(p[, gof.name]), na.rm=TRUE))
      
      sp::coordinates(p) <- ~ x+y
      
      # Reduced margins among figures. From: https://stat.ethz.ch/pipermail/r-help/2007-January/123556.html
      theme.novpadding <- list(
                              layout.heights= list(
                                         top.padding = 2, main.key.padding = 0,
                                 	 key.axis.padding= 0, axis.xlab.padding= 0,
                                 	 xlab.key.padding= 0, key.sub.padding= 0,
                                 	 bottom.padding= 0
                                 	 ),
                                layout.widths=list(
                                         left.padding= 3, key.ylab.padding= 0,
                         	         ylab.axis.padding= 0, axis.key.padding=0,
                         	         right.padding= 0
                         	         )
                                )
      sp::spplot(p, scales=list(draw=TRUE, cex=0.75, tick.number=4, x=list(rot=axis.rot[1]), 
                 y= list(rot=axis.rot[2]) ), cuts=GOFcuts, col.regions=colorRamp(13), 
                 aspect="fill", auto.key=auto.key, key.space= key.space, 
                 xlab=list(p1.name, font=2), ylab=list(p2.name, font=2), 
                 main=main, cex=points.cex, alpha=alpha, 
                 par.settings = theme.novpadding)
  
  } else if (type=="scatter3d") {

        x <- params[1:n,p1.name]
        y <- params[1:n,p2.name]
        z <- params[1:n,gof.name]

        scatterplot3d::scatterplot3d(x, y, z, highlight.3d=TRUE, col.axis="grey", 
                      type="p", col.grid="lightblue", pch=20, main=main,
                      xlab=p1.name, ylab=p2.name, zlab=gof.name)
    } # IF end

} # 'plot_2parOF' END
