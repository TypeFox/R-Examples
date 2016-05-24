#' Creates a pdf file with the design of the Paper Helicopter
#' 
#' The pdf file contains a template with lines and indications to build the 
#' paper helicopter described in many SixSigma publications.
#' 
#' The pdf file must be printed in A4 paper, without adjusting size to paper.
#' 
#' @return 
#' No value is returned. A pdf file is saved in the working directory
#' 
#' @references 
#' George Box.
#' Teaching engineers experimental design with a paper helicopter.
#' \emph{Quality Engineering}, 4(3):453--459, 1992.
#' 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.
#' 
#' @note 
#' See the \code{vignette("HelicopterInstructions")} to see assembling instructions.
#' 
#' @author EL Cano
#' 
#' @examples 
#' ss.heli()
#' vignette("HelicopterInstructions")
#' @export
ss.heli <- function(){
	pdf(file="helicopter.pdf", width=6, height=10)
	#dev.new(width=8, height=15)
	#grid.newpage()
	grid::grid.rect(gp=grid::gpar(fill="white"))
	grid::grid.text("Six Sigma with R | Paper Helicopter template",
			y=unit(1,"npc")-unit(0.5,"cm"), just="top")
	vpframe<-grid::viewport(width=unit(12,"cm"), height=unit(22.5,"cm"))
	grid::pushViewport(vpframe)
	grid::grid.rect(gp=grid::gpar(lwd=2))
	grid::grid.text("cut",x=unit(-2,"mm"), just=c("left","bottom"),rot=90)
	grid::grid.rect(y=1, height=unit(9.5,"cm"), just="top")
	grid::grid.text(expression(fold %up%""), x=unit(3,"cm"), y=unit(13.3,"cm"))
	grid::grid.text(expression("fold"%down%""), x=unit(9,"cm"), y=unit(13.3,"cm"))
	grid::grid.lines(x=c(0.5,0.5), y=unit(c(13,22.5),"cm"), gp=grid::gpar(lty=2))
	grid::pushViewport(grid::viewport(y=unit(13,"cm"), width=unit(2,"cm") , 
					height=unit(2,"cm"),angle=45))
	grid::grid.rect( gp=grid::gpar(lty=3))
	grid::grid.text("tape?", y=unit(1,"npc")-unit(2,"mm"), 
			gp=grid::gpar(cex=0.7, col="gray"), just="top")
	grid::popViewport()
	grid::grid.text("cut",x=unit(6.2,"cm"),y=unit(15,"cm"), rot=90, 
			just=c("left", "top"))
	grid::grid.rect(x=0, y=0, width=unit(4,"cm"), height=unit(9.5,"cm"),
			just=c("left","bottom"), gp=grid::gpar(lty=2))
	grid::grid.text(expression("fold"%down%""%down%""), 
			x=unit(3.5,"cm"), y=unit(1,"cm"), rot=90)
	grid::grid.text("cut", x=unit(2,"cm"), y=unit(9,"cm"))
	grid::grid.lines(x=unit(c(4,4),"cm"), y=unit(c(0,9.5),"cm"))
	grid::grid.rect(x=1, y=0, width=unit(4,"cm"), height=unit(9.5,"cm"),
			just=c("right","bottom"), gp=grid::gpar(lty=2))
	grid::grid.text(expression("fold"%up%""%up%""), 
			x=unit(8.5,"cm"), y=unit(1,"cm"), rot=90)
	grid::grid.text("cut", x=unit(10,"cm"), y=unit(9,"cm"))
	grid::grid.lines(x=unit(c(8,8),"cm"), y=unit(c(0,9.5),"cm"))
	grid::grid.rect(x=unit(4,"cm"), y=unit(2,"cm"), width=unit(1.5,"cm"), 
			height=unit(6,"cm"), gp=grid::gpar(lty=3), just="bottom")
	grid::grid.text("tape?", x=unit(42,"mm"), y=unit(5,"cm"), 
			gp=grid::gpar(cex=0.7, col="gray"), just="top", rot=90)
	grid::grid.rect(x=unit(8,"cm"), y=unit(2,"cm"), width=unit(1.5,"cm"), 
			height=unit(6,"cm"), gp=grid::gpar(lty=3), just="bottom")
	grid::grid.text("tape?", x=unit(76,"mm"), y=unit(5,"cm"), 
			gp=grid::gpar(cex=0.7, col="gray"), just="top", rot=90)
	
	grid::grid.xspline(x=unit(c(5.5,5.6,6.6,6.5),"cm"), 
			y=unit(c(0,2,2,0),"cm"),
			shape=c(0,-1,1,0),
			gp=grid::gpar(lty=3))
	
	grid::grid.text("clip?", x=0.5, y=unit(2,"mm"), gp=grid::gpar(cex=0.7, 
					col="gray"), just="bottom")
	#Body length limits
	grid::grid.lines(x=c(0,1), y=unit(c(3,3),"cm"), 
			gp=grid::gpar(lty=4, col="gray"))
	grid::grid.text("min\n(6.5cm)", x=unit(12.1,"cm"), y=unit(3,"cm"), 
			gp=grid::gpar(cex=0.7), 
			just="left")
	grid::grid.lines(x=c(0,1), y=unit(c(1.5,1.5),"cm"), 
			gp=grid::gpar(lty=4, col="gray"))
	grid::grid.text("std\n(8cm)", x=unit(12.1,"cm"), y=unit(1.5,"cm"), 
			gp=grid::gpar(cex=0.7), just="left")
	grid::grid.text("max\n(9.5cm)", x=unit(12.1,"cm"), y=unit(0,"cm"), 
			gp=grid::gpar(cex=0.7), just="left")
	grid::grid.text(expression(""%<-%" body length "%->%""), 
			x=unit(12.1,"cm"), y=unit(6,"cm"), rot=90, just="top")
	
	#Body width limits
	grid::grid.text(expression(""%<-%" body width "%->%""), 
			y=unit(-0.1,"cm"), just="top")
	grid::grid.text("min\n(4cm)", x=unit(4,"cm"), 
			y=unit(-0.1,"cm"), just="top", gp=grid::gpar(cex=0.7))
	grid::grid.text("min\n(4cm)", x=unit(8,"cm"), 
			y=unit(-0.1,"cm"), just="top", gp=grid::gpar(cex=0.7))
	grid::grid.text("max\n(6cm)", x=unit(3,"cm"), 
			y=unit(-0.1,"cm"), just="top", gp=grid::gpar(cex=0.7))
	grid::grid.text("max\n(6cm)", x=unit(9,"cm"), 
			y=unit(-0.1,"cm"), just="top", gp=grid::gpar(cex=0.7))
	grid::grid.lines(x=unit(c(3,3),"cm"), y=unit(c(0,9.5),"cm"), 
			gp=grid::gpar(lty=4, col="gray"))
	grid::grid.lines(x=unit(c(9,9),"cm"), y=unit(c(0,9.5),"cm"), 
			gp=grid::gpar(lty=4, col="gray"))
	grid::grid.lines(x=unit(c(3.5,3.5),"cm"), y=unit(c(0,9.5),"cm"), 
			gp=grid::gpar(lty=4, col="gray"))
	grid::grid.lines(x=unit(c(8.5,8.5),"cm"), y=unit(c(0,9.5),"cm"), 
			gp=grid::gpar(lty=4, col="gray"))
	
	#Wings length limits
	grid::grid.lines(x=c(0,1), y=unit(c(19.5,19.5),"cm"), 
			gp=grid::gpar(lty=4, col="gray"))
	grid::grid.text("min\n(6.5cm)", x=unit(12.1,"cm"), y=unit(19.5,"cm"), 
			gp=grid::gpar(cex=0.7), 
			just="left")
	grid::grid.lines(x=c(0,1), y=unit(c(21,21),"cm"), 
			gp=grid::gpar(lty=4, col="gray"))
	grid::grid.text("std\n(8cm)", x=unit(12.1,"cm"), y=unit(21,"cm"), 
			gp=grid::gpar(cex=0.7), just="left")
	grid::grid.text("max\n(9.5cm)", x=unit(12.1,"cm"), y=1, 
			gp=grid::gpar(cex=0.7), just="left")
	grid::grid.text(expression(""%<-%" wings length "%->%""), 
			x=unit(12.1,"cm"), y=unit(16,"cm"), rot=90, just="top")
	dev.off()
	
#	PostScriptTrace("helides.ps")
#	design <- readPicture("helides.ps.xml")
#	grid.picture(design)
	
}


