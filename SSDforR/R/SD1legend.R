SD1legend <-
function(){

par(mar=c(1, 1, 1, 1))
plot.new()
legend("center", c("behavior","+1sd","mean","-1sd"), col = c("red","black", "green","black"), lwd = 1,ncol=4,bty ="n")
}
