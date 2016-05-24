f.debug.pvalues <- function(.effects.CI, .pvalues.new, .pvalues.old, .arg.effects, .names.haplo){
##
##
## Debug for beregning av p-verdier ved Johnsons fordeling (pakke SuppDists)
##
#
## Output skal i eget directory...
dir.create("debug.haplin.pvalues", showWarnings=F)
#
# require(SuppDists) # Denne skal skiftes ut senere...
#
pdf(paste("debug.haplin.pvalues/",paste(.names.haplo,sep="",collapse="-"),".pdf",sep=""))
on.exit(dev.off(),add=T)
par(mfrow=c(3,3))
#
if(length(.arg.effects!=0)){
	.temp.pvalues <- sapply(colnames(.arg.effects),function(x){
		## Johnson parameters
		.param <- JohnsonFit(.arg.effects[,x], moment="quant")
		#
		## Plotter den empriske kumulative mot Johnson 
		.emp.cum <- ecdf(.arg.effects[,x])
		plot(.emp.cum, main = x)
		.arg.Johnson <- seq(max(min(.arg.effects[,x])*0.9,.Machine$double.xmin),min(max(.arg.effects[,x])*1.1,.Machine$double.xmax),by=0.01)
		.pJohnson <- pJohnson(.arg.Johnson,.param)
		lines(.arg.Johnson,.pJohnson, col="red")
	})
}
#
## Estimat, konfidensintervaller, nye og gamle p-verdier skrives til fil, sammen med differansen
write.table(cbind(.effects.CI, p.value.new = as.numeric(.pvalues.new), p.value.old = as.numeric(.pvalues.old), diff.pvalues = (as.numeric(.pvalues.new)-as.numeric(.pvalues.old))), paste("debug.haplin.pvalues/",paste(.names.haplo,sep="",collapse="-"),".txt",sep=""))
#
## Return empty
return(invisible())
#
}