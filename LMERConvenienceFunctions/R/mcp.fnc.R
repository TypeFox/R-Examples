mcp.fnc <-
function(model,trim=2.5,col="red"){
	data<-model@frame
	par(mfrow=c(2,2))
	data$rstand = as.vector(scale(resid(model))) 
	plot(density(data$rstand),main="") 
	qqnorm(data$rstand,pch = ".",main ="") 
	qqline(data$rstand,col=col)
	plot(data$rstand~fitted(model),pch=".")
	abline(h = c(-trim, trim),col=col)
	### removed this until can figure out how to calculate dffits 
	### for a merMod. If you know how, please let the package 
	### maintainer know :-)
	#dffits = abs(resid(model, "dffits")) 
	#plot(dffits, type="h")
	par(mfrow=c(1,1))
}
