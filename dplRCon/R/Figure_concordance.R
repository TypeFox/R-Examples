#' Figure: Plots the Concordances 
#'
#' Plots the concordance indices. This is used to generate the Figure in Concordance: A measure of similarity between matrices of time series with applications in dendroclimatology".
#' @param con.indices This is the output generated from the function overall.concordance.period. 
#' @param x.lim This is a vector with the start and end year for the required plot. 
#' @author Maryann Pirie
#' @return A figure with the concordance indices and their confidence intervals, a smoother fitted through the indices and the confidence interval for the mean concordance. 
#' @examples 
#' \dontrun{
#' #Subset "near-pith" is the material within 0 -20cm from the estimated pith
#' spline200.sub0.20.n   <- TruncSeriesPithoffset( ring.raw, ring.stand, dbh.po.nc, c(1,200))
# Subset "far-pith" is the material further than 20cm from the estimated pith
#' spline200.sub20.2000.n  <- TruncSeriesPithoffset( ring.raw, ring.stand, dbh.po.nc, c(200,200000))

#series.bootstapped
#' boot.0.20   <-  series.bootstrap( spline200.sub0.20.n$sub.series.stand, stat, 999,
#'     names.stat, aver.by.tree = FALSE)
#' boot.20.2000   <- series.bootstrap(spline200.sub20.2000.n$sub.series.stand, stat, 999, 
#'    names.stat, aver.by.tree = FALSE)
#'
#' overall.precision.HUP <- overall.concordance.period(spline200.sub20.2000.n$sub.series.stand , 
#'    spline200.sub0.20.n$sub.series.stand, 
#'    boot.20.2000$boot.series.mean,  boot.0.20$boot.series.mean ,1 , concordance.indices, 
#'    c(1880,1999), trim.alpha=0.005, concordance.beta=0.5)
#' 
#' figure.function.concordance (overall.precision.HUP, x.lim=c(1880,2000))}
#' @export
#' 
figure.function.concordance <- function(con.indices, x.lim ){
  #  plot.new()
  #	par( mar = c(5, 4, 1, 2)) 
  plot(con.indices$pre.in.period[,3], 
       type = 'n', 
       main = "",      
       lwd = 2, 
       xlab = "Year",
       ylab = "Concordance indices",
       yaxs = 'i',
       ylim = c(0, 1),
       xlim= x.lim  , 
       las=1
  )
  time <- seq(tsp(con.indices$ci.p)[1], length = length(con.indices$ci.p))
  points(ts(con.indices$pre.in.period[,3], start = con.indices$period[1]), pch = '+', cex = 0.75)
  ci <- NULL
  for( i in 1:length(time)){
    arrows(time[i], con.indices$ci.p[[i]][1], time[i] , con.indices$ci.p[[i]][2], angle=90, length = 0 )
  }
  lines(smooth.spline(na.contiguous(con.indices$pre.in.period[,3]), spar=0.3), col=1, lwd = 3)     
  abline(h=0.8, lty = 2)
  
  mtext(paste("Overall Concordance = (", round(con.indices$ci.con[1], 4), ", ", round(con.indices$ci.con[2], 4), ")" ), 
        side = 1, line = 3, adj = 0, cex = 0.8)
  
}
