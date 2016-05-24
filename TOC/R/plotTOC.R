.plotTOC <- function(object, labelThres=FALSE, modelLeg="Model", digits=3, nticks=5, digitsL=1, posL = NULL, offsetL = 0.5, ...){

  old.opt <- options()
  options(digits=digits)
  old.par <- par(no.readonly = TRUE)
  par(oma = c(0, 0, 0, 4))
  par(mgp = c(1.5, 1, 0))
  
population <- object@population
prevalence <- object@prevalence/population
units <- object@units

tocd <- object@table
if((!is.null(tocd$HitsP) & !is.null(tocd$"Hits+FalseAlarmsP"))==TRUE){
tocd$Hits <- tocd$HitsP
tocd$"Hits+FalseAlarms" <- tocd$"Hits+FalseAlarmsP"
}

graphics::plot(c(0, population*(1-prevalence), population), c(0, 0, prevalence * population), type="l", lty="dashed", 
     xlab=paste0("Hits+False Alarms (", units, ")"), ylab=paste0("Hits (", units, ")"), 
     lwd=2, col=rgb(128,100,162, maxColorValue=255), bty="n", xaxt="n", yaxt="n", xlim=c(0, population), 
     ylim=c(0, prevalence * population), asp=1/prevalence, ...)

xlabels <- c(0, format((1:nticks)*population/nticks, digits))
ylabels <- c(0, format((1:nticks)*prevalence * population/nticks, digits))

axis(1, pos = 0, labels=xlabels, at=xlabels, xaxp = c(0, population, nticks), cex.axis=0.9, ...)
axis(2, pos = 0, labels=ylabels, at=ylabels, yaxp = c(0, prevalence * population, nticks), cex.axis=0.9, ...)


# maximum
lines(c(0, prevalence * population, population), c(0, prevalence * population, prevalence * population), 
      lty="dotdash", lwd=2, col=rgb(79,129,189, maxColorValue=255)) 

# hits+misses
lines(c(0, population), rep(prevalence*population, 2), lwd=3, col=rgb(146,208,80, maxColorValue=255))

# uniform
lines(c(0, population), c(0, prevalence*population), lty="dotted", lwd=2, col=rgb(0,0,255, maxColorValue=255))

#lines(tocd$"Hits+FalseAlarms", tocd$maximum, lty="dotdash", lwd=2, col=rgb(79,129,189, maxColorValue=255))

# model
lines(tocd$"Hits+FalseAlarms", tocd$Hits, lwd=2, col=rgb(255,0,0, maxColorValue=255))
points(tocd$"Hits+FalseAlarms", tocd$Hits, pch=17, col=rgb(255,0,0, maxColorValue=255))
if(labelThres == TRUE) text(tocd$"Hits+FalseAlarms", tocd$Hits, round(as.numeric(tocd$Threshold), digitsL), pos = posL, offset = offsetL, ...)

par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
graphics::plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")

legend("right", c("Hits+Misses", "Maximum", modelLeg, "Uniform", "Minimum"), 
       col = c(rgb(146,208,80, maxColorValue=255), rgb(79,129,189, maxColorValue=255), rgb(255,0,0, maxColorValue=255), rgb(0,0,255, maxColorValue=255), rgb(128,100,162, maxColorValue=255)), 
       lty = c(1, 4, 1, 3, 2), pch = c(NA, NA, 17, NA, NA),
       merge = TRUE, bty="n", lwd=c(3, 2, 2, 2, 2))
par(old.par)
options(old.opt)

}
