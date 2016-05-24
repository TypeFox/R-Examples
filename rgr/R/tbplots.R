tbplots <-
function(x, by, log = FALSE, logx = FALSE, notch = TRUE, xlab = "", 
	ylab = deparse(substitute(x)), ylim = NULL, main = "", 
	label = NULL, plot.order = NULL, xpos = NA, width, space = 0.25, 
	las = 1, cex = 1, adj = 0.5, add = FALSE, ssll = 1, colr = 8, 
	...)
{
     # Original script to plot box-and-whisker plots shared by Doug Nychka on 
     # S-News, April 28, 1992; script modified to generate an IDEAS style box-
     # and-whisker plots, Garrett, June 15, 1995.  Additional modifications
     # made in September 1999, and to plot Tukey boxplots in October 2003,
     # and in August 2006 during conversion to R.
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
     if(!add) frame()
     quant <- c(0, 0.25, 0.5, 0.75, 1)
     if(logx) log <- TRUE
     if(log)
         logplot <- "y"
     else logplot <- ""
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
                 range.data, log = logplot, type = "n", xaxt = "n", ylab = ylab, 
                 xlab = xlab, main = main, ...)
         else 
             plot(range(c(xpos - (0.5 * width)/space, xpos + (0.5 * width)/space)),
                 range.data, log = logplot, type = "n", xaxt = "n", ylab = ylab, 
                 xlab = xlab, ylim = ylim, main = main, ...)
     }
     xmax <- 10^12
     for(i in 1:cols) {
         temp <- data[[plot.order[i]]]
         temp <- na.omit(temp)
         sssz <- length(temp)
         if(sssz >= ssll) {
          bb <- quantile(temp, quant)
          if(log&&logx) {
               bb.save <- bb
               bb <- log10(bb)
          }
          hw <- bb[4] - bb[2]
          lwl <- bb[2] - 1.5 * hw
          lol <- bb[2] - 3 * hw
          uwl <- bb[4] + 1.5 * hw
          uol <- bb[4] + 3 * hw
          if(log&&logx) {
               lwl <- 10^lwl
               lol <- 10^lol
               uwl <- 10^uwl
               uol <- 10^uol
               bb <- bb.save
          }
          mid <- xpos[i]
          lwend <- min(temp[temp > lwl])
          uwend <- max(temp[temp < uwl])
          for(j in 1:sssz) {
               if((temp[j] < lwend) || (temp[j] > uwend)) {
                    pch <- 3
                    pchcex <- 1
                    if((temp[j] <= lol) || (temp[j] >= uol)) {
                         pch <- 1
                         pchcex <- 1.8
                    }
                    points(mid, temp[j], pch = pch, cex = pchcex)
               }
          }
          low <- mid - width[i] * 0.5
          hih <- mid + width[i] * 0.5
          if(sssz > 5) {
               x <- c(mid, mid, NA, mid, mid)
               y <- c(lwend, bb[2], NA, bb[4], uwend)
               lines(x, y)
               if(notch) {
                    v <- sort(temp)
                    j <- qbinom(0.025, sssz, 0.5)
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
               lines(x, y)
               lines(c(hih, low), c(bb[3], bb[3]))
          }
         }
     }
     #    if((las = 1) & (length(label) > 7) & missing(label))
     #         las <- 2
     axis(1, xpos, label, tick = FALSE, las = las, cex.axis = cex, adj = adj)
     invisible()
}
