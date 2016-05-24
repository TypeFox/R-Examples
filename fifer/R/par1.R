##' Change default par parameters
##'
##' @title Change default par parameters
##' @export
##' @author Dustin Fife
par1 = function(){
	if (.Platform$OS.type=="unix"){
		par(mar=c(3.25,3.25,1,1), mgp=c(2, .5, 0), tck=-.01, cex.axis=.8, family="Times", cex.lab=1.5, cex.main=1.5, font.main=1)
	} else {
		par(mar=c(3.25,3.25,1,1), mgp=c(2, .5, 0), tck=-.01, cex.axis=.8, cex.lab=1.5, cex.main=1.5, font.main=1)
	}
}


##' Change default par parameters
##'
##' @title Change default par parameters
##' @export
##' @author Dustin Fife
par2 = function(){
	par1()
	if (.Platform$OS.type=="unix"){	
		par(family="Helvetica", font.lab=2, cex.axis=1)
	} else {
		par(font.lab=2, cex.axis=1)
	}
	par(mar=c(3.25,3.25,1,1))
	par(cex.lab =1.5)	
}



