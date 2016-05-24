#' @rdname LSC-utils
#' @description
#' \code{\link{plot_LSC_2plus1D}} plots LSC for a (2+1)D field.
#' @keywords hplot
#' @param z an object of class \code{"LSC"}
#' @param type a \code{"temporal"} or a \code{"spatial"} summary plot of LSC
#' @param time.frames a vector of length \eqn{\leq 6} to indicate what frames
#' should be displayed (only for \code{type = "temporal"}). If \code{NULL} 
#' (default) then it chooses them automatically based on valleys and peaks in the 
#' spatial average LSC.
#' @param zlim minimum and maximum z values for which colors should be 
#' plotted, defaulting to the range of the finite values of \code{z}.
#' @param lsc.unit character string (default: \code{"bits"}) to write 
#' next to the color legend
#' @param col colors: either a string decribing a pallette from the 
#' \code{RColorBrewer} package (see also \url{http://colorbrewer2.org/}), or a list of colors 
#' (see \code{\link[graphics]{image}} for suggestions).
#' @param data (optional) original data to compare to LSC (relevant only for \code{type = "spatial"})
#' @param heights passed to \code{\link[graphics]{layout}} for dividing the plotting
#' region vertically. If \code{data = NULL} a vector of length 2; otherwise a
#' vector of length 3.
#' @export
#' @examples
#' \dontrun{
#' data(contCA00)
#' 
#' temp_lsc = states2LSC(states = contCA00$predictive_states - min(contCA00$predictive_states) + 1)
#' temp_lsc_3D = array(temp_lsc, dim = c(25, 20, 40))
#' class(temp_lsc_3D) = c("LSC", "LSC_2plus1D")
#' plot_LSC_2plus1D(temp_lsc_3D, type = "temporal")
#' plot_LSC_2plus1D(temp_lsc_3D, type = "spatial")
#' }
#'
#' 
# 

plot_LSC_2plus1D <- function(z, type = "temporal", 
                             time.frames = NULL, 
                             zlim = NULL, 
                             heights = NULL, 
                             lsc.unit = "bits", 
                             data = NULL, col = NULL){
  
  LSC <- z

  if (!is.null(data)) {
    cc <- (dim(data)[1] - dim(LSC)[1] ) / (dim(data)[3] - dim(LSC)[3] )  
    h.p <- (dim(data)[3] - dim(LSC)[3] )  / cc
  }
  
  if (is.null(col)) {
    col <- colorRampPalette(brewer.pal(9, name="YlOrRd"))(100)
  } 
  if (length(col) == 1 & is.character(col)){
    col <- colorRampPalette(brewer.pal(9, name=col))(100)
  }

  LSC.colors <- col
  orig.colors <- tim.colors(100)
  
  if (type == "temporal") {
    LSC.temporal.avg <- exp(apply(log(LSC), 3, mean))
    LSC.temporal.median <- apply(LSC, 3, median) 
    LSC.temporal.mean.trimmed <- apply(LSC, 3, mean, trim = 0.1) 
    LSC.temporal.sd <- apply(LSC, 3, function(mat) sd(c(mat))) 
    
    
    TT <- length(LSC.temporal.avg)
    #LSC.temporal.avg.smooth <- loess(LSC.temporal.avg~c(1:TT), span = 0.2)$fit
    LSC.temporal.avg.smooth <- gam(LSC.temporal.avg~s(c(1:TT)) )$fitted.values
    LSC.temporal.median.smooth <- gam(LSC.temporal.median~s(c(1:TT)) )$fitted.values
    LSC.temporal.sd.smooth <- gam(LSC.temporal.sd~s(c(1:TT)) )$fitted.values
    
    if (is.null(time.frames)) {
      if (is.null(heights)) {
        if (is.null(data)) {
          heights <- c(3, 4) 
        } else {
          heights <- c(2, 2, 4)
        }
      }
      time.frames.tmp <- c(1, which.min(LSC.temporal.avg), which.max(LSC.temporal.avg), TT)
      for (split in seq_len(floor(2 * log2(TT)))) {
        window.size <- floor(TT/split)
        temp.LSC <- LSC.temporal.avg.smooth[(split-1)*window.size + seq_len(window.size)]
        time.frames.tmp <- c(time.frames.tmp, which.max(temp.LSC), which.min(temp.LSC))
      }
      time.frames.tmp <- sort(unique(time.frames.tmp), decreasing = FALSE)
      
      time.frames <- time.frames.tmp[1]
      for (tt in 2:length(time.frames.tmp)){
        if (time.frames[length(time.frames)] + 2*log2(TT) < time.frames.tmp[tt]) { 
          time.frames <- c(time.frames, time.frames.tmp[tt])
        }
      }
    }
    
    if (length(time.frames) > 5) {
      time.frames <- time.frames[sort(sample(1:length(time.frames), 5), decreasing = FALSE)]
    }
    num.frames <- length(time.frames)
    par(cex.lab = 2, cex.axis = 2, lwd = 1)
    if (!is.null(data)) {
      layout(matrix(c(seq_len(2*(num.frames + 1)), rep(2*(num.frames + 1) + 1, num.frames + 1)), 
             byrow = TRUE, nrow = 3), heights = heights)#, widths = c( 5 rep(4, num.frames), 3 ) )
    } else {
      layout(matrix(c(1, 2:(num.frames + 1), num.frames + 2, rep((num.frames + 2) + 1, num.frames)), 
                    byrow = TRUE, nrow = 2), heights = heights, widths = c(2, rep(4, num.frames) ) )
    }
    
    LR.margin <- 0.2
    TB.margin <- 0.2

    if (!is.null(data)){
      # data plots
      par(mar = c(TB.margin, LR.margin, TB.margin + 1, 5 - LR.margin - 1), cex.lab = 2, 
          cex.axis = 2, lwd = 1)
      image2(data[, ,tt+h.p], col = orig.colors, legend = "only", main = "", 
             zlim = range(data[, ,time.frames]))
      
      par(mar = c(TB.margin, LR.margin, TB.margin, LR.margin), cex.lab = 2, cex.axis = 2, lwd = 1)
      for (tt in time.frames){
        image2(data[, ,tt+h.p], col = orig.colors, legend = FALSE, main = "", 
               zlim = range(data[, ,time.frames]))
      }
    }
    # LSC plots
    
    if (is.null(zlim)) {
      zlim <- range(LSC[, ,time.frames])
    }
    par(mar = c(TB.margin + 2, LR.margin + 3, TB.margin + 1, LR.margin * 2), 
        cex.lab = 2, cex.axis = 2, lwd = 1)
    make_legend(data = 0, col = LSC.colors, zlim = zlim, side = 2, col.label = lsc.unit)
    
    par(mar = c(TB.margin ,LR.margin,TB.margin,LR.margin),cex.lab = 2, cex.axis = 2, lwd = 1)
    for (tt in time.frames) {
      image2(LSC[, ,tt], col = LSC.colors, legend = FALSE, main = "", zlim = zlim)
      mtext(paste(tt), side = 1, line = TB.margin+1, cex=1.25)
    }
    par(mar = rep(0, 4), cex.lab = 2, cex.axis = 2, lwd = 2)
    plot.new()
    
    par(mar = c(4, LR.margin, 3, LR.margin), cex.lab = 2, cex.axis = 2, lwd = 2)
    
    plot(LSC.temporal.avg, ylab = "", xlab = "", axes = FALSE, type="l")
    lines(LSC.temporal.avg.smooth, lwd = 2, col = "red")
    par(new=TRUE)
    axis(2)
    box()
    axis(1, at  = time.frames, col = "blue", lty=2)
    grid()
    axis(4)
    mtext("Time", side = 1, cex=1.5, line = 3)
    mtext("bits", side = 2, cex=1.5, line = 3)
    
    #plot(LSC.temporal.sd, axes=FALSE, type = "l", col = "red", lwd=1)
    #lines(LSC.temporal.sd.smooth, lwd = 2, col = "red")
    #lines(LSC.temporal.mean.trimmed, lty = 2, lwd=1)
    #lines(LSC.temporal.median, lty = 2, lwd = 1)
    #lines(LSC.temporal.median.smooth, lty = 2, lwd = 1, col = "red")
    
    points(time.frames, LSC.temporal.avg[time.frames], pch = 15, col = "blue", cex = 2)
    abline(v = time.frames, lty = 2, col = "blue")
    #usr <- par("usr")
    #text(time.frames+2, rep(usr[3] + 0.05*(usr[4]-usr[3]), num.frames), paste(time.frames), cex = 1.5)
    #par(op)
    
    invisible(LSC.temporal.avg)
  }
  
  if (type == "spatial") {  
    LSC.temporal.avg <- apply(LSC, 3, mean, na.rm = TRUE)
    TT <- length(LSC.temporal.avg)
    LSC.temporal.avg.smooth <- gam(LSC.temporal.avg~s(c(1:TT)))$fitted.values
    
    #plot(LSC.temporal.avg)
    #lines(LSC.temporal.avg.smooth, col = 2)
    
    weights <- LSC.temporal.avg.smooth / sum(LSC.temporal.avg.smooth)
    weights <- LSC.temporal.avg / sum(LSC.temporal.avg)
    
    
    LSC.spatial <- list()
    LSC.spatial[["median"]] <- apply(LSC, 1:2, median, na.rm = TRUE) 
    LSC.spatial[["weighted mean"]] <- apply(LSC, 1:2, weighted.mean, w = weights, na.rm = TRUE) 
    LSC.spatial[["mean"]] <- apply(LSC, 1:2, mean, na.rm = TRUE) 
    LSC.spatial[["sd"]] <- apply(LSC, 1:2, sd, na.rm = TRUE)
    
    color.pal <- list()
    color.pal[["mean"]] <- LSC.colors
    color.pal[["weighted mean"]] <- color.pal[["mean"]]
    color.pal[["median"]] <- color.pal[["mean"]]
    color.pal[["sd"]] <- colorRampPalette(brewer.pal(9, name="PuBu"))(100)
    
    zlim.spatial <- list()
    mean.range <- range(LSC.spatial[["mean"]], LSC.spatial[["weighted mean"]], LSC.spatial[["median"]])
    
    zlim.spatial[["mean"]] <- zlim.spatial[["median"]] <- zlim.spatial[["weighted mean"]] <- mean.range
    for (name in c("sd")) {
      zlim.spatial[[name]] <- range(LSC.spatial[[name]])
    }
    
    if (!is.null(zlim)) {
      zlim.spatial <- zlim
    }
    if (is.null(heights)) {
      heights <- c(10, 3) 
    }

    num.summaries <- length(LSC.spatial)
    #  layout(matrix(c(rep(c(1,2,3), each = 2), 3 + 1:(2*3)), byrow = TRUE, nrow = 2),
    #         widths = rep(c(5,2), 3), heights = c(1,3) )
    
    par(mar = c(0, 2, 0.5, 1))
    layout(matrix(seq_len(num.summaries*2), byrow = FALSE, ncol = num.summaries),
           heights = heights )
    
    
    niter <- 0
    for (name in names(LSC.spatial)) {
      niter <- niter + 1
      par(cex.lab = 1.5, cex.axis = 1.5, lwd = 1)
      #par(mar = c(0,1,0.5,1))
      if (!is.null(zlim)) {
        temp.pdf <- density( LSC.spatial[[name]], from = zlim.spatial[[name]][1], 
                             to = zlim.spatial[[name]][2] )
      } else { 
        temp.pdf <- density( LSC.spatial[[name]] )
        
      }
      
      #plot(temp.pdf, main = "", xlab = "", 
      #     xlim = zlim.spatial[[name]], ylab = "pdf", axes = FALSE)
      
      #abline(v = range(LSC.spatial[[name]]), col = "gray")
      #box()
      
      #image2(LSC.spatial[[name]], col =color.pal[[name]], legend = "only", zlim = zlim.spatial[[name]])
      
      par(mar = c(0.5, 1, 3, 1))
      image2(LSC.spatial[[name]], col = color.pal[[name]], legend = FALSE, zlim = zlim.spatial[[name]])
      mtext(paste0("(", letters[niter], ") ", name), side = 3, line = 1, cex = 1.25)
      
      par(mar = c(4, 1, 0, 1))
      make_legend(data=LSC.spatial[[name]], col = color.pal[[name]], side = 1, 
                  zlim = zlim.spatial[[name]], 
                  col.ticks = pretty(temp.pdf$x, n = 6), cex.axis = 1.75)
      temp.pdf$y <- temp.pdf$y / max(temp.pdf$y) * 0.9
      lines(temp.pdf, lwd = 2)
      mtext("bits", side = 1, line = 3, cex = 1.25)
    }
  }
}
