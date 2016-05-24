plot.dwt.multiple <- function( x, levels = NULL, ylim = NULL, draw.dashed.lines = TRUE, draw.level.labels = TRUE, col = c("red","blue"), ...)
{
    inner.plot.level <- function(dwt.object, w.range, v.range, col = c("red", "blue")) {
        prev.plotcolor <- col[2]
        index <- 0
   
        for(i in w.range) {
            plotted.level <- FALSE

            if(prev.plotcolor == col[2] && plotted.level == FALSE) {
                prev.plotcolor <- col[1]

                lines((index+1):(index + length(dwt.object@W[[i]][,1])), as.vector(dwt.object@W[[i]][,1]), col = prev.plotcolor)

                index <- index + length(dwt.object@W[[i]][,1])
                plotted.level <- TRUE
            }

            if(prev.plotcolor == col[1] && plotted.level == FALSE) {
                prev.plotcolor <- col[2]

                lines((index+1):(index + length(dwt.object@W[[i]][,1])), as.vector(dwt.object@W[[i]][,1]), col = prev.plotcolor)

                index <- index + length(dwt.object@W[[i]][,1])
                plotted.level <- TRUE
            }
        }
     
        for(i in v.range) {
            plotted.level <- FALSE
            if(prev.plotcolor == col[2] && plotted.level == FALSE) {
                prev.plotcolor <- col[1]

                lines((index+1):(index + length(dwt.object@V[[i]][,1])), as.vector(dwt.object@V[[i]][,1]), col = prev.plotcolor)

                index <- index + length(dwt.object@V[[i]][,1]) 
                plotted.level <- TRUE
            }

            if(prev.plotcolor == col[1] && plotted.level == FALSE) {
                prev.plotcolor <- col[2]

                lines((index+1):(index + length(dwt.object@V[[i]][,1])), as.vector(dwt.object@V[[i]][,1]), col = prev.plotcolor)
 
                index <- index + length(dwt.object@V[[i]][,1])
                plotted.level <- TRUE

            }
        }    
    }

    inner.plot <- function(levels, range, filter.name = NULL, type = "h", level.labels = NULL, ylim = NULL, xlim = NULL, ylab = NULL, xlab = "", draw.dashed.lines = TRUE, draw.xaxis.label = FALSE)
    {

        draw.dash <- function(range) 
        {
            range.shift <- 0
            for(i in 1:length(range)) {  
                abline(v=range.shift+range[i], lty = 2)       
                range.shift <- range.shift + range[i]
            }
        } 

        draw.level.labels <- function(level.labels, range) 
        {
             range.shift <- 0
            if(length(level.labels) > 3) {
                level.labels <- c(level.labels[1:3], "...")
            }
             for(i in 1:length(level.labels)) {
                 mtext(level.labels[[i]], side = 3, at = range.shift + range[[i]]/2, line = 1)
                 range.shift <- range.shift + range[[i]]
             }
        }

        if(is.null(ylab)) {
            ylab = filter.name
        }
        plot(1:length(levels), y = levels, type = type, ylab = ylab, axes = FALSE, frame.plot = TRUE, 
        ylim = ylim, xaxs = "i", xlim = c(0,length(levels)))     
        axis(2)

        if(draw.xaxis.label) {
            axis(1, labels = TRUE)
        }
        else {
            axis(1, labels = FALSE)
        }

        if(is.null(xlab) == FALSE) {
            mtext(xlab, side = 1, line = 2)    
        }

        if(draw.dashed.lines) {
            par(new = FALSE)
            draw.dash(range)
        }

        if(is.null(level.labels) == FALSE) {
           par(new=FALSE)
           draw.level.labels(level.labels, range)              
        }

    }
  
    concatenate.levels <- function( x, w.range, v.range)
    {
        list.levels <- list()

        for(i in 1:length( x)) {
            dwt.object <-  x[[i]]
            concatenated.levels <- NULL
            for(j in w.range)
            {
                concatenated.levels <- c(concatenated.levels, as.vector(dwt.object@W[[j]][,1]))
            }

            for(j in v.range)
            {
                concatenated.levels <- c(concatenated.levels, as.vector(dwt.object@W[[j]][,1])) 
            }
            list.one.level <- list(concatenated.levels)
            list.levels <- c(list.levels, list.one.level)
        }

        list.levels    
    }    

    range.levels <- function(dwt.object, w.range, v.range)
    {
        range.levels <- NULL
        for(j in w.range)
        {
            range.levels <- c(range.levels, length(dwt.object@W[[j]][,1]))   
        }
  
        for(j in v.range)
        {
            range.levels <- c(range.levels, length(dwt.object@W[[j]][,1]))   
        }

        range.levels    
    }    

    if(class( x) =="dwt") {
         x <- list( x)
    }
    if (class( x) != "list") { 
        stop("Invalid argument: ' x' must be of class 'dwt', 'list'")
    }
    else {
        for(i in 1:length( x)) {
            if(class( x[[i]]) != "dwt") {
                stop("Invalid argument: elements of ' x' must be of class 'dwt'")
            }
        }
    }   

    if(is.null(levels)) {
        minlevel <-  x[[1]]@level
        for(i in 1:length( x)) {
            minlevel <- min(minlevel,  x[[i]]@level)
        }
        
        w.range <- 1:minlevel
        v.range <- max(w.range)  
    }
    if(class(levels) == "numeric") {
        if(length(levels) == 1) {
        w.range <- 1:levels
        v.range <- max(w.range)  
        }
        else {
        w.range <- levels
        v.range <- max(w.range)
        }
    }
    if(class(levels) == "integer") {
        if(length(levels) == 1) {
        w.range <- 1:levels
        v.range <- max(w.range)  
        }
        else {
        w.range <- levels
        v.range <- max(w.range)
        }
    }
    if(class(levels) == "list") {
        w.range <- levels[[1]]
        v.range <- levels[[2]]
    }

    par(omi = c(.6, 0, .4, 0))
    par(mfrow = c(length( x),1))

    concatenated.levels <- concatenate.levels( x, w.range, v.range)
    range <- range.levels( x[[1]], w.range, v.range)

    if(is.null(ylim)) {

        minimum <- min(concatenated.levels[[1]])
        maximum <- max(concatenated.levels[[1]])
        if(length(concatenated.levels) > 1) {
            for(i in 2:length(concatenated.levels)) {
                if(minimum > min(concatenated.levels[[i]])) {
                    minimum <- min(concatenated.levels[[i]])
                }	
            }
       
            for(i in 2:length(concatenated.levels)) {
                if(maximum < max(concatenated.levels[[i]])) {
                    maximum <- max(concatenated.levels[[i]])
                }	
            }
        }
        ylim <- c(floor(minimum), ceiling(maximum))    
    }

    if(draw.level.labels) {
        level.labels <- NULL
        for (j in 1:length(w.range)) {
            label <- substitute(W[level], list(level = w.range[j]))
            level.labels <- c(level.labels, label)
        }

        for (i in 1:length(v.range)) {
            label <- substitute(V[level], list(level = v.range[i]))
            level.labels <- c(level.labels, label)
        }
    }
    else {
        level.labels <- NULL
    }

    for(i in 1:length( x)) {

        levels <- concatenated.levels[[i]]  
        if(i == 1 && length( x) == 1) {
             par(mai = c(0, .8, 0, .4))
             inner.plot(levels, range,  x[[i]]@filter@wt.name, type = "n", level.labels = level.labels, ylim = ylim, draw.dashed.lines = draw.dashed.lines, draw.xaxis.label= TRUE, xlab = "n")
             inner.plot.level( x[[i]], w.range, v.range)
        }
     
        if(i == 1 && length( x) != 1) {
             par(mai = c(.1, .8, 0, .4))
             inner.plot(levels, range,  x[[i]]@filter@wt.name, type = "n", level.labels = level.labels, ylim = ylim, draw.dashed.lines = draw.dashed.lines, draw.xaxis.label = FALSE)
             inner.plot.level( x[[i]], w.range, v.range, col = col)
        }
        if(i != 1 && i != length( x)) {
             par(mai = c(.1, .8, 0, .4))
             inner.plot(levels, range,  x[[i]]@filter@wt.name, type = "n", ylim = ylim, draw.dashed.lines = draw.dashed.lines, draw.xaxis.label = FALSE)        
             inner.plot.level( x[[i]], w.range, v.range, col = col)
        }
        if(i == length( x) && i != 1) {
             par(mai = c(.1, .8, 0, .4))
             inner.plot(levels, range,  x[[i]]@filter@wt.name, type = "n", ylim = ylim, draw.dashed.lines = draw.dashed.lines, draw.xaxis.label = TRUE, xlab = "n")
             inner.plot.level( x[[i]], w.range, v.range, col = col)
        }
    }    
}

