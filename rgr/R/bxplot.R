bxplot <-
function(xx, xlab = deparse(substitute(xx)), log = FALSE, ifbw = FALSE, wend = 0.05, 
     xlim = NULL, main = "", ifn = TRUE, colr = 8, cex = 1, ...)
{
     # Function to plot a single horizontal box plot; the default is a Tukey boxplot, 
     # ifbw = T generates an IDEAS style box-and-whisker plot, and wend defines the
     # the end of the whisker - default is 5th and 95th %ile.  For the special case
     # of wend = 0, the whiskers extend to the data minima and maxima and no plusses
     # are plotted at the extremes; when wend = 0.25 only the box s plotted.  Setting
     # log = T generates log-scaled plots and a log transformation of the data for the
     # Tukey boxplot outlier calculations.  The box is infilled with a grey, colr = 8,
     # if no colour is required, set colr = 0.  Settng xlim may result in outliers not
     # being plotted as the x-axis is shortened; however, the statistics used to define
     # the boxplot, or box-and-whisker plot, are still based on the total data set.
     # The data set size is printed in the lower right corner, as is the number of
     # omitted observations if setting xlim has resulted in truncation.  To plot just
     # part of a data set create a subset first, or use the x[x<some.value] construct
     # in the call.
     #
     # NOTE: Prior to using this function the data frame/matrix containing the
     # variable, 'x', data must be run through ltdl.fix.df to convert any <dl
     # -ve values to positive half that value, and set zero2na = TRUE if it is
     # required, to convert any zero values or other numeric codes representing 
     # blanks to NAs.
     #
     temp.x <- remove.na(xx)
     x <- temp.x$x[1:temp.x$n]
     nx <- temp.x$n
     rangex <- range(x)
     if(log) {
         logplot <- "x"
         if((!is.null(xlim)) && (xlim[1] <= 0))
             xlim[1] <- rangex[1]
     }
     else logplot <- ""
     if(is.null(xlim)) {
         plot(x = rangex, y = c(0, 0), xlab = xlab, ylab = "", log = logplot,
             ylim = c(-1, 1), yaxt = "n", type = "n", main = main, ...)
         limits <- par("usr")
         if(log)
             limits[1:2] <- 10^limits[1:2]
         nxx <- 0
     }
     else {
         plot(x = rangex, y = c(0, 0), xlab = xlab, ylab = "", log = logplot,
             xlim = xlim, ylim = c(-1, 1), yaxt = "n", type = "n", main = main, ...)
         limits <- par("usr")
         if(log)
             limits[1:2] <- 10^limits[1:2]
         dropped <- x[(x > xlim[2]) | (x < xlim[1])]
         nxx <- length(dropped)
     }
     if(ifbw) {
         if(wend <= 9.9999999999999995e-008) {
          lowend <- rangex[1]
          hihend <- rangex[2]
          q <- as.vector(quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
          q1 <- q[1]
          q2 <- q[2]
          q3 <- q[3]
         }
         else {
          if(wend > 0.25)
               wend <- 0.05
          q <- as.vector(quantile(x, probs = c(wend, 0.25, 0.5, 0.75, 1 - wend), 
               na.rm = TRUE))
          lowend <- q[1]
          q1 <- q[2]
          q2 <- q[3]
          q3 <- q[4]
          hihend <- q[5]
          if(wend < 0.25)
               points(rangex, c(0, 0), pch = 3)
         }
     }
     else {
         q <- as.vector(quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
         xcut <- numeric(4)
         if(log) {
          lq1 <- log10(q[1])
          lq3 <- log10(q[3])
          hw <- lq3 - lq1
          xcut[1] <- 10^(lq1 - 3 * hw)
          xcut[2] <- 10^(lq1 - 1.5 * hw)
          xcut[3] <- 10^(lq3 + 1.5 * hw)
          xcut[4] <- 10^(lq3 + 3 * hw)
         }
         else {
          hw <- q[3] - q[1]
          xcut[1] <- q[1] - 3 * hw
          xcut[2] <- q[1] - 1.5 * hw
          xcut[3] <- q[3] + 1.5 * hw
          xcut[4] <- q[3] + 3 * hw
         }
         lowend <- min(x[x > xcut[2]])
         q1 <- q[1]
         q2 <- q[2]
         q3 <- q[3]
         hihend <- max(x[x < xcut[3]])
         for(i in 1:nx) {
             if((x[i] < lowend) || (x[i] > hihend)) {
                 pch <- 3
                 sym.cex <- 1
                 if((x[i] <= xcut[1]) || (x[i] >= xcut[4])) {
                      pch <- 1
                      sym.cex <- 1.6
                 }
                 if(is.null(xlim) || ((x[i] <= xlim[2]) & (x[i] >= xlim[1]))) {
                    points(x[i], 0, pch = pch, cex = sym.cex)
                 }
             }
         }
     }
     polygon(c(q1, q1, q3, q3, q1), c(0.4, -0.4, -0.4, 0.4, 0.4), col = colr)
     ypos <- c(0, 0, 0.4, 0.4, 0, 0, 0, -0.4, -0.4, 0.4, -0.4, -0.4, 0)
     xpos <- c(lowend, q1, q1, q3, q3, hihend, q3, q3, q2, q2, q2, q1, q1)
     lines(xpos, ypos)
     if(ifn) {
         if(log)
             xpos <- 10^(log10(limits[2]) - (log10(limits[2]) - log10(limits[1])) * 
             0.05)
         else xpos <- limits[2] - (limits[2] - limits[1]) * 0.05
         ypos <- limits[3] + (limits[4] - limits[3]) * 0.15
         text(xpos, ypos, labels = paste("N =", nx), adj = 1, cex = cex)
         if(nxx != 0) {
             ypos <- limits[3] + (limits[4] - limits[3]) * 0.06
             text(xpos, ypos, labels = paste(nxx, "points omitted"), adj = 1,
                  cex = cex * 0.8)
         }
     }
     invisible()
}
