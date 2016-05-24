#' Produces a Figure of the correlation coefficients. 
#' 
#' Produce a figure with 3 panels, each panel is for a different climate variable. An example of this figure in included in "On the influence of tree size on the climate - growth relationship of New Zealand kauri (Agathis australis): insights from annual, monthly and daily growth patterns. J Wunder, AM Fowler, ER Cook, M Pirie, SPJ McCloskey. Trees 27 (4), 937-948"
#' @param corr.1  the output from using function 'correlation.function' for the 1st climate variable comparing near and far-pith. Depending on the order of the inputs into the function 'correlation.function' will determine the color of the resultant box-plot (1st item is black, and 2nd item is gray) 
#' @param corr.2 same as corr.1 but for the second climate variable
#' @param corr.3 same as corr.1 but for the third climate variable
#' @param corr.1.full the output from using function 'correlation.function' for the full dataset. correlations for the full dataset are shown as a grey dashed line. 
#' @param corr.2.full same as corr.1.full but for the second climate variable
#' @param corr.3.full same as corr.1.full but for the third climate variable
#' @param col.names.season  col.names.season<- list("SON_2", "DJF_2", "MAM_2", "JJA_2", "SON_1", "DJF_1", "MAM_1", "JJA_1", "SON", "DJF", "MAM", "JJA")
#' @return This returns a figure. 
#' @examples 
#' \dontrun{ Figure.correlation.barplot(corr.temp, corr.prec, corr.SOI, 
#'    corr.temp.full, corr.prec.full, corr.SOI.full, col.names.season) }
#' @export

Figure.correlation.barplot <- function(corr.1, corr.2, corr.3, corr.1.full, corr.2.full, corr.3.full, col.names.season) {
  
  # This sub-function draws the pox plot. 
  # 'summary' is the input vector containing c(min, lower quartile, median, upper quartile, max) for correlation functions generated from the bootstrapped chronologies.
  # the other inputs control the colour and position of the box-plot 
  #plot.box <- function(summary, at, col.i, offset){
  #  	xx <- as.numeric(substr(summary[,at], 9, 18))
  #		arrows(at+offset, xx[1], at+offset, xx[2], length = 0, col = col.i, lwd =2)
  #		arrows(at+offset, xx[5], at+offset, xx[6], length = 0, col = col.i, lwd =2)
  #		rect(at-0.15+offset, xx[2], at+0.15+offset, xx[5], border = col.i)
  #		arrows(at-0.15+offset, xx[3], at+0.15+offset, xx[3], length = 0, lwd = 2, 
  #col = col.i)
  #}
  
  plot.bar <- function(corrx, at, col.i, offset){
    xx <- corrx[at]
    rect(at-0.15+offset, 0, at+0.15+offset, xx, border = col.i, col = col.i)
  }
  
  
  # plots a star relating the the strength of evidence of a difference between the correlation functions constructed from the two chronologies. 
  # corrx: is a list outputed from 'correlation.function'.
  plot.sign <- function(corrx){
    prb <- pt(abs(corrx$t.meanequal), df=99 , lower.tail=FALSE)*2  # calculates the p-value, $t.meansequal is the test statistic for H0: mean correlation from subset x = mean correlation from subset y.
    for ( i in 4:12){    # plots points based on level of significance
      if(prb[i] <= 0.001){points(i-0.05, -0.4, pch = '*', col = "red", cex = 2)
                          points(i+0.05, -0.4, pch = '*', col = "red", cex =2)
      }else{
        if(prb[i] <= 0.05 ){points(i, -0.4, pch = '*', col = "red", cex = 2)
                            #      }else{
                            #  if(prb[i] <= 0.1) {points(i, -0.3, pch = '*', col = "red")}
        }}}}
  # ******* The following was added******  
  plot.sign.full <- function(corrx){
    #	plot(NULL, xlim = c(1,12), ylim = c(0,1), yaxt = "n", xaxt = "n", xlab="", ylab="", bty = "n")
    col.i <- c(gray(0), gray(0.5))
    offset.i <- c(-0.20, 0.20)
    for(ss in 1:2){
      # $t.mean is the test statistics for H0: mean correlation for subset = 0
      prb <- pt(abs(corrx$t.mean[ss,]), df=99 , lower.tail=FALSE)*2
      for ( i in 4:12){
        if( corrx$t.mean[ss,i] >= 0) {
          if(prb[i] <= 0.001){points(i + offset.i[ss] -0.07, +0.4, pch = '+', col = col.i[ss], cex = 1.5)
                              points(i+ offset.i[ss] +0.07, +0.4, pch = '+', col = col.i[ss], cex = 1.5)
          }else{
            if(prb[i] <= 0.05 ){points(i+offset.i[ss], +0.4, pch = '+', col = col.i[ss], cex =1.5)
            }} }
        if( corrx$t.mean[ss,i] < 0) {
          if(prb[i] <= 0.001){points(i + offset.i[ss] -0.075, +0.4, pch = '-', cex = 1.9, col = col.i[ss])
                              points(i+ offset.i[ss] +0.075, +0.4, pch = '-', col = col.i[ss], cex = 1.9)
          }else{
            if(prb[i] <= 0.05 ){points(i+offset.i[ss], +0.4, pch = '-', col = col.i[ss], cex=1.9)
                                
                                #          }else{
                                #	    if(prb[i] <= 0.1) {points(i+offset.i[ss], +0.5, pch = '*', col = col.i[ss])}
            }}}
      }} }
  # ******* End of change*******
  
  # This function using the above two sub-functions and performs the ploting. 
  # x: 'summary' for the first subset
  # y: 'summary' for the second subset
  # full: 'summary for the full dataset
  # corrx: used to construct the strength of evidemce of a difference between x and y
  plot.fun.in <- function(corrx, col.names.season){
    plot(NULL, xlim=c(4.8, 12.2), ylim=c(-0.5, 0.5), yaxs = 'i', xaxt = "n", xlab="", ylab="", las = 1, cex.lab = 1, cex.axis = 1.5)
    rect(8.5, -0.5, 11.5, 0.5, density = NULL, angle = 45,
         col = gray(0.9), border = gray(0.9), lty = NULL, lwd = par("lwd"),
         xpd = NULL)
    rect(0, -0.5, 12.5, 0.5, density = NULL, angle = 0,
         col = NULL, border = 1, lty = NULL, lwd = par("lwd"),
         xpd = NULL)
    axis(1,  at = 5:12,  las=2, labels = rep("", length(col.names.season[5:12])) )
    mtext("Correlation coefficients", side=2, line=4, cex= 1)
    #mtext("Season", side=1, line=3, cex = 0.75)
    abline(h=0)
    abline(h=0.5)
    #lines(full,  col= gray(0.75), lwd = 2, lty = 2)
    
    for( i in 4:12){
      plot.bar(corrx$corr.site.1, i, col.i = gray(0), offset = -0.15)
    }
    for( i in 4:12){
      plot.bar(corrx$corr.site.2, i, col.i = gray(0.5), offset = 0.15)
    }
    
    plot.sign(corrx)
    plot.sign.full(corrx)
    #axis(1,  at = 5:12,  las=1, labels =unlist( col.names.season)[5:12], cex = 0.75)
    #axis(1,  at = 5:12,  las=1, labels =rep(expression(paste("SON",scriptstyle(t-1))),8), cex = 0.75)	
    axis(1,  at = 5:12,  las=1, labels =c(expression(paste("SON ",scriptstyle(t-1))), expression(paste("DJF ",scriptstyle(t-1))), 
                                          expression(paste("MAM ",scriptstyle(t-1))), expression(paste("JJA ",scriptstyle(t-1))),
                                          expression(paste("SON ",scriptstyle(t))), expression(paste("DJF ",scriptstyle(t))),
                                          expression(paste("MAM ",scriptstyle(t))), expression(paste("JJA ",scriptstyle(t))))  , cex.lab= 4, cex.axis = 1.5)	
  }
  layout( matrix(seq(1,3), 3, 1), heights = c(4,4,4))   # ****
  # The follows uses the above functions to construct the plots, and produced the three pannels for each of the climate varaibles 
  # ***** par(mfrow = c(3,1), mar= c(4, 4, 2, 1)+0.1)   # figure layout
  par(mar= c(4, 5.5, 3, 1)+0.1)   #****
  plot.fun.in(corr.1, col.names.season)
  
  mtext("a", 3, adj = -0.7*(par("mar")[2] * par("csi"))/(dev.size()[1] - sum(par("mar")[c(2, 4)]) * par("csi")), line = 0.75, cex=1.5)
  #mtext("a", 3, adj = +11.5*(par("mar")[2] * par("csi"))/(dev.size()[1] - sum(par("mar")[c(2, 4)]) * par("csi")), line = 0.75, cex=1.5)
  
  plot.fun.in(corr.2, col.names.season)
  mtext("b", 3, adj = -0.7*(par("mar")[2] * par("csi"))/(dev.size()[1] - sum(par("mar")[c(2, 4)]) * par("csi")), line = 0.75, cex=1.5)
  #mtext("b", 3, adj = 11.5*(par("mar")[2] * par("csi"))/(dev.size()[1] - sum(par("mar")[c(2, 4)]) * par("csi")), line = 0.75, cex=1.5)
  
  plot.fun.in(corr.3, col.names.season)
  mtext("c", 3, adj = -0.7*(par("mar")[2] * par("csi"))/(dev.size()[1] - sum(par("mar")[c(2, 4)]) * par("csi")), line = 0.75, cex=1.5)
  #mtext("c", 3, adj = 11.5*(par("mar")[2] * par("csi"))/(dev.size()[1] - sum(par("mar")[c(2, 4)]) * par("csi")), line = 0.75, cex=1.5)
  mtext("Season", side=1, line=3, cex = 1)
}