bwplots <-
function(x, by, log = FALSE, wend = 0.05, notch = TRUE, xlab = "",
	ylab = deparse(substitute(x)), ylim = NULL, main = "", 
	label = NULL, plot.order = NULL, xpos = NA, width, space = 0.25, 
	las = 1, cex.axis = 1, adj = 0.5, add = FALSE, ssll = 1, colr = 8, 
	pch = 3, ...)
{
     # Original script to plot box-and-whisker plots shared by Doug Nychka on
     # S-News, April 28, 1992; script modified to generate IDEAS style box-
     # and-whisker plots, Garrett, June 15, 1995.  Additional modifications made
     # in September 1999, on other ocassions, and in August 2006 during
     # conversion to R.
     #
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variable, 'x', data must be run through ltdl.fix.df to convert any <dl
     # -ve values to positive half that value, and set zero2na = TRUE if it is
     # required, to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     if(is.matrix(x)) data <- data.frame(x)
     if(data.class(x) == "numeric")
         data <- list(x, ...)
     if(is.list(x))
         data <- x
     if(!missing(by))
         data <- cat2list(unlist(data), by)
     # 
     # at this point the data should be in list form regardless of how 
     # the pieces were originally passed
     #
     if(wend > 0.25) wend <- 0.05
     quant <- c(wend, 0.25, 0.5, 0.75, 1 - wend)
     if(log)
         logy <- "y"
     else logy <- ""
     cols <- length(data)
     range.data <- range(as.numeric(unlist(data)), na.rm = TRUE)
     if(missing(label)) {
         if(is.null(names(data)))
          label <- format(1:cols)
         else label <- names(data)
     }
     if(is.null(plot.order))
         plot.order <- 1:cols
     else {
         cat(" label ", label, "\n", "plot.order ", plot.order, "\n")
         labels <- label
         for(i in 1:cols)
          label[i] <- labels[plot.order[i]]
         cat(" label ", label, "\n")
     }
     if(is.na(xpos)) {
         xpos <- 1:cols
     }
     if(missing(width)) {
         width <- min(diff(sort(xpos))) * space
         if(cols == 1)
          width <- space
     }
     if(length(width) == 1)
         width <- rep(width, cols)
     if(!add) {
         if(is.null(ylim))
             plot(range(c(xpos - (0.5 * width)/space, xpos + (0.5 * width)/space)), 
                 range.data, log = logy, type = "n", xaxt = "n", ylab = ylab, 
                 xlab = xlab, main = main, ...)
         else 
             plot(range(c(xpos - (0.5 * width)/space, xpos + (0.5 * width)/space)),
                 range.data, log = logy, type = "n", xaxt = "n", ylab = ylab, 
                 xlab = xlab, ylim = ylim, main = main, ...)
     }
     for(i in 1:cols) {
         temp <- data[[plot.order[i]]]
         temp <- temp[!is.na(temp)]
         sssz <- length(temp)
         if(sssz >= ssll) {
          bb <- quantile(temp, quant)
          mid <- xpos[i]
          low <- mid - width[i] * 0.5
          hih <- mid + width[i] * 0.5
          if(sssz > 5) {
               x <- c(mid, mid, NA, mid, mid)
               y <- c(bb[1], bb[2], NA, bb[4], bb[5])
               lines(x, y)
               if(notch) {
                    v <- sort(temp)
                    j <- qbinom(0.025, sssz, 0.5)
                    if (j < 1) j <- 1
                    lci <- v[j]
                    uci <- v[sssz - j + 1]
                    x <- c(hih, low, low, mid, low, low, hih, hih, mid, hih, hih)
                    y <- c(bb[2], bb[2], lci, bb[3], uci, bb[4], bb[4], uci, bb[
                         3], lci, bb[2])
                    polygon(x, y, col = colr)
                    lines(x, y)
               }
               else {
                    x <- c(hih, low, low, hih, hih)
                    y <- c(bb[2], bb[2], bb[4], bb[4], bb[2])
                    polygon(x, y, col = colr)
                    lines(x, y)
                    lines(c(hih, low), c(bb[3], bb[3]))
               }
          }
          else {
               x <- c(hih, low, low, hih, hih)
               y <- c(bb[2], bb[2], bb[4], bb[4], bb[2])
               polygon(x, y, col = colr)
               lines(x, y)
               lines(c(hih, low), c(bb[3], bb[3]))
          }
          if(wend >= 9.9999999999999995e-008) {
               if(wend <0.25) {
                   points(mid, min(temp), pch = pch)
                   points(mid, max(temp), pch = pch)
               }
          }
         }
     }
     #    if((length(label) > 7) & missing(label))
     #         las <- 2
     axis(1, xpos, label, tick = FALSE, las = las, cex.axis = cex.axis, adj = adj)
     invisible()
}
