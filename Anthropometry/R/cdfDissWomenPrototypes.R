cdfDissWomenPrototypes <- function(min_med,min_med_UNE,main,xlab,ylab,leg,cexLeg,...){
 plot(sort(min_med), (1:length(min_med))/length(min_med), type = "s", main = main, xlab = xlab, ylab = ylab, 
      xlim = c(0,1), ylim = c(0,1), yaxt = "n",...)
 axis(2,at = seq(0,1,0.1),labels = seq(0,1,0.1)) 
 lines(sort(min_med_UNE), (1:length(min_med_UNE))/length(min_med_UNE), type="s", lty = "dashed",...) 
 legend("bottomright", legend = leg, lty = c("solid", "dashed"), cex = cexLeg)
}
