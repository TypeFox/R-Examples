SD2legend <-
function(){
  
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("behavior","+2sd","mean","-2sd"), col = c("red","black", "green","black"),lwd = 3,ncol=4,bty ="n")
}
