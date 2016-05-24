#' Plot global data
#' 
#' Plot a quick world map with reasonable coloring.
#'
#' @param x A \code{\link{cmip5data}} object
#' @param dates numeric. Which date value(s) should we plot?
#' @param splitPacific logical. Try to split image in the Pacific?
#' @param capMinMax logical. Cap data min and max by quantile? This may produce better coloring.
#' @param verbose logical. Print info as we go?
#' @return A ggplot object.
#' @details Uses \code{ggplot2::geom_raster}.
#' @examples
#' d <- cmip5data(1970:1975)   # sample data
#' worldPlot(d)
#' @export
worldPlot <- function(x, dates=unique(x$val$time), splitPacific=TRUE, capMinMax=TRUE, verbose=FALSE) {
 
    # Sanity checks
    stopifnot(class(x)=="cmip5data")
    stopifnot(is.numeric(dates)) 
    stopifnot(is.logical(capMinMax) & length(capMinMax)==1)
    stopifnot(is.logical(verbose) & length(verbose)==1)
    length(dates) <- min(length(dates), 16)   # can't see anything smaller...
    stopifnot(require(ggplot2))
    
    # Preliminaries
    lon <- x$lon
    lat <- x$lat
    val <- x$val
    
    # Filter for Z and time
    if(length(unique(val$Z)) > 1) {
        # Suppress stupid NOTEs from R CMD CHECK
        Z <- NULL        
        val <- dplyr::filter(val, Z == min(val$Z))  # only use the first lev/depth        
    }
    
    if(length(unique(val$time)) > 1) {
        # Suppress stupid NOTEs from R CMD CHECK
        time <- NULL        
        val <- dplyr::filter(val, time %in% dates)       
    }
    
    if(nrow(val) == 0) {
        warning("No data found after date filtering!")
        return()
    }
    
    # Split at the Pacific ocean
    if(splitPacific) {
        half.numlon <- ceiling(length(lon) * 0.527) # yields a division at 190 for a 360 grid   
        shiftIndex <- c(half.numlon:length(lon), 1:(half.numlon-1))        
        if(verbose) cat("Splitting Pacific at longitude position", half.numlon, "(", lon[half.numlon], ")\n")
        val$lon <- lon[shiftIndex[match(val$lon, lon)]]
    }
    
    if(capMinMax) {
        if(verbose) cat('Setting quantile-derived min/max\n')
        quant <- quantile(val$value, na.rm=TRUE)
        minNum <- quant[[3]]-10*(quant[[3]]-quant[[2]])
        maxNum <- quant[[3]]+10*(quant[[4]]-quant[[3]])
        minNum <- max(minNum, quant[[1]])
        maxNum <- min(maxNum, quant[[5]])
        val$value[val$value < minNum] <- minNum
        val$value[val$value > maxNum] <- maxNum        
    }
    
    # Plot
    
    # Suppress stupid NOTEs from R CMD CHECK
    value <- NULL        
    nrowcol <- sqrt(length(unique(val$time)))
    
    p <- ggplot2::ggplot(val, ggplot2::aes(lon, lat))
    p <- p + ggplot2::geom_raster(ggplot2::aes(fill=value))
    p <- p + ggplot2::scale_x_continuous(expand=c(0,0))
    p <- p + ggplot2::scale_y_continuous(expand=c(0,0))
    p <- p + ggplot2::scale_fill_gradientn(colours=rainbow(4))
    p <- p + ggplot2::facet_wrap(~time, nrow=nrowcol, ncol=nrowcol)
    p + ggplot2::ggtitle(paste0(x$model, " ", x$experiment, " ", x$variable, " (", x$valUnit, ")"))
} # worldPlot2
