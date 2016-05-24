# script for producing Bayes factor figures

# to call: bayfacfig(2, 5, c(-5.4534, -5.3955, -5.235, -5.19948, -5.321), 4)

bayfacfig <- function(indnr, modelnr, BF, markmod)
{  
  model_nr = c(1:modelnr)
  
  dev.set(1)
  postscript("bferrbars.eps", horizontal=FALSE, width=5, height=5,
             onefile=FALSE, paper="special", family="ComputerModern")
  
  hE <- errbar(model_nr, BF[model_nr], BF[model_nr]+rep(0,length(BF[model_nr])), 
               BF[model_nr]-rep(0,length(BF[model_nr])), lty=3, type='b',
               lwd=1, pch=19, col="blue", xlab="Model(number of terms)", ylab="Log Bayes Factor")
  
  # highlighling selected model
  points(markmod, BF[markmod], pch=19, col="red")
  
  legend('topleft', legend = c("Selected Model"), pch=19, col = "red", bty="n")

  dev.off()  
}  