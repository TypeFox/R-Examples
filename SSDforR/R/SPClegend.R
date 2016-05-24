SPClegend <-
function(){
  
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("behavior","Uband","mean","Lband"), col = c("red","blue", "green","orange"), lwd = 1,ncol=4,bty ="n")
}
