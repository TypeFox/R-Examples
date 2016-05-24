level_shardsplot <- function(
    object, par.names,rows = 1:NCOL(object$data),
    centers = rep(NA, length(par.names)),
    class.labels = NA, 
    revert.colors = rep(FALSE, length(par.names)),
    log.classes = rep(FALSE, length(par.names)),
    centeredcolors = colorRamp(c("red", "white", "blue")),
    mfrow = c(2,2), plot.type = c("eight", "four", "points", "n"),
    expand = 1, stck = TRUE, grd = FALSE, standardize = FALSE, 
    label = FALSE, plot = TRUE, vertices = TRUE, 
    classcolors = "topo", wghts = 0, 
    xlab = "Dimension 1", ylab = "Dimension 2", xaxs = "i", 
    yaxs = "i", ...)
{

    require("som")
    if(class(object) != "som") stop("Object must be of class som")
    classes <- 100
    oldpar <- par(mfrow = mfrow)
    
    for(i in rows){
        values <- object$code[,i]
        val.class <- rep(classes, length(values))
        class.breaks <- seq(from = min(values), to = max(values), length.out = classes+1)
        if(!is.na(centers[i])){
            if(log.classes[i]){
                    col.values <- log(values)
                    centers[i] <- log(centers[i])
            } else {
                    col.values <- values
            }
            myrgb <- centeredcolors(rerange(seq(from = min(col.values), to = max(col.values), length.out = 100), 
                center = centers[i]))
    
            theClasscolors <- rgb(red = myrgb[,1], green = myrgb[,2], blue = myrgb[,3], maxColorValue = 255)
            which.center <- 100 * (centers[i] - min(col.values)) / (max(col.values) - min(col.values))
            if(which.center < 1 || which.center > 100) 
                which.center <- NA
            legend.colors <- theClasscolors[c(1, which.center, classes)]
        } else {
            if(is.character(classcolors) && length(classcolors) == 1){
                theClasscolors <- switch(classcolors,
                    "rainbow" = rainbow(classes),
                    "topo"    = topo.colors(classes),
                    "gray"    = gray((1:classes) / classes),
                    stop("argument classcolors only supports 'rainbow', 'topo', and 'gray'")
                )
                if(revert.colors[i]){
                    theClasscolors <- rev(theClasscolors)
                }
            } else stop("argument classcolors must be 'rainbow', 'topo', and 'gray'")
        legend.colors <- theClasscolors[c(1, classes / 2, classes)]
        }
    
        shardsplot(object = object, plot.type = plot.type, expand = expand,
            stck = stck, grd = grd, standardize = standardize, label = label,
            plot = plot, vertices = vertices, classcolors = theClasscolors,
            wghts = wghts, xlab = xlab, ylab = ylab, xaxs = xaxs, yaxs = yaxs,
            plot.data.column = i, log.classes = log.classes[i],
            revert.colors = revert.colors[i], xaxt = "n", yaxt = "n",...)
    
        if(any(!is.na(class.labels))){
            opar <- par(xpd = NA)
            legend(((object$ydim-1)/2)+1, -0.2, pt.bg = legend.colors, 
                xjust = 0.5, yjust = 1, legend = class.labels[i,], 
                pch = 21, horiz = TRUE, col = "black", pt.cex = 1.3)
            par(opar)    
        }
        mtext(3,text = par.names[i], line = 1, adj = 0.5, cex = 1.5)
    }
    par(oldpar)
}
