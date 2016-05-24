.onLoad <- function(lib, pkg) {
	invisible(.C("FPUcheck",PACKAGE="PearsonDS"))
}

#Removed (because of NOTE in R-2.15.0, resulting slow down should be negligible)
#.onAttach <- function(lib, pkg) {
# 	if (suppressWarnings(require(gsl,quietly=TRUE))) {
#    assign(".hasGSL",TRUE,pos=paste("package:",pkg,sep=""))
#  } else {
#    assign(".hasGSL",FALSE,pos=paste("package:",pkg,sep=""))
#  }   
#}
