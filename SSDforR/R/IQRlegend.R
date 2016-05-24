IQRlegend <-
function(){
par(mar=c(1, 1, 1, 1))
plot.new()
legend("center", c("behavior","p75","median","p25"), col = c("red","blue", "green","orange"), lwd = 1,ncol=4,bty ="n")
}
