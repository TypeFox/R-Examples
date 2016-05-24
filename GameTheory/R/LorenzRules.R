LorenzRules <-
function(x){
	DATA<-x[1:nrow(x)-1,]
	par(family = "Times")
	plot(Lc(DATA[,2]),family="Times",col="#63A5DB",axis=NULL,lwd=1.5,lty=4)
	lines(Lc(DATA[,3]),col="#005B9A",lty=4,lwd=1)
	lines(Lc(DATA[,4]),col="#734A75",lty=4,lwd=1)
	lines(Lc(DATA[,5]),col="#B80606",lwd=1,lty=2)
	lines(Lc(DATA[,6]),col="#E38030",lwd=1,lty=5)
	legend(0, 1,
	 c("Proportional", "CEA", "CEL","Talmud","RA"), 
	 col = c("#63A5DB", "#005B9A", "#734A75","#B80606","#E38030"),
       text.col = "black",
       lty = c(4, 4, 4,2,5), 
       cex=0.85,
       merge = TRUE, 
       bg = "gray95")
}
