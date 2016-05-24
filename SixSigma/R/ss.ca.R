# Capability Analysis Functions
# 
# Author: Emilio Lopez
###############################################################################
if(getRversion() >= '2.15.1') utils::globalVariables(c("..density..", "value"))


#' Main calculations regarding The Voice of the Process in SixSigma: Yield, FTY, RTY,
#' DPMO
#' 
#' Computes the Yield, First Time Yield, Rolled Throughput Yield and Defects
#' per Million Opportunities of a process.
#' 
#' The three arguments must have the same length.
#' 
#' @param defects A vector with the number of defects in each product/batch, ...
#' @param rework A vector with the number of items/parts reworked
#' @param opportunities A numeric value with the size or length of the product/batch
#' @return 
#'   \item{Yield }{Number of good stuff / Total items}
#'   \item{FTY }{(Total - scrap - rework) / Total }
#'   \item{RTY }{prod(FTY)}
#'   \item{DPMO}{Defects per Million Opportunities}
#' @references 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.
#' 
#' Gygi C, DeCarlo N, Williams B (2005) \emph{Six sigma for dummies}. --For dummies,
#'   Wiley Pub.
#' @author Emilio L. Cano
#' @export
#' @examples 
#' ss.ca.yield(c(3,5,12),c(1,2,4),1915)
ss.ca.yield <- function(defects = 0, rework = 0, opportunities = 1){
	Yield <- (opportunities - sum(defects)) / opportunities
	FTY <- (opportunities - sum(defects) - sum(rework)) / opportunities
	RTY <- prod((opportunities - (defects + rework)) / opportunities)
	DPU <- sum(defects)
	DPMO <- (DPU / opportunities) * 10^6 
	ss.ca.yield <- (list(Yield = Yield, FTY = FTY, 
						RTY = RTY, DPU = DPU, DPMO = DPMO))
	as.data.frame(ss.ca.yield)
} 



#' Capability Indices
#' 
#' Compute the Capability Indices of a process, Z (Sigma Score), \eqn{C_p} 
#' and \eqn{C_{pk}}.
#' 
#' @usage 
#' ss.ca.cp(x, LSL = NA, USL = NA, LT = FALSE, f.na.rm = TRUE, 
#'   ci = FALSE, alpha = 0.05)
#' @usage 
#' ss.ca.cpk(x, LSL = NA, USL = NA, LT = FALSE, f.na.rm = TRUE, 
#'   ci = FALSE, alpha = 0.05)
#' @usage 
#' ss.ca.z(x, LSL = NA, USL = NA, LT = FALSE, f.na.rm = TRUE)
#' 
#' @aliases ss.ca.z ss.ca.cp ss.ca.cpk
#' 
#' @param x A vector with the data of the process performance
#' @param LSL Lower Specification Limit
#' @param USL Upper Specification Limit
#' @param LT Long Term data (TRUE/FALSE). Not used for the moment
#' @param f.na.rm Remove NA data (TRUE/FALSE)
#' @param ci If TRUE computes a Confidence Interval
#' @param alpha Type I error (\eqn{\alpha}) for the Confidence Interval 
#' @return A numeric value for the index, or a vector with the limits 
#' of the Confidence Interval
#' 
#' @references 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.
#' 
#' Montgomery, DC (2008) \emph{Introduction to Statistical Quality Control}
#'   (Sixth Edition). New York: Wiley&Sons\cr
#' @author EL Cano
#' @seealso \code{\link{ss.study.ca}}
#' @keywords cp cpk capability
#' @examples 
#' ss.ca.cp(ss.data.ca$Volume,740, 760)
#' ss.ca.cpk(ss.data.ca$Volume,740, 760)
#' ss.ca.z(ss.data.ca$Volume,740,760)
#' @export ss.ca.z ss.ca.cp ss.ca.cpk 
ss.ca.z <- function(x, LSL = NA, USL = NA, 
		LT = FALSE, f.na.rm = TRUE){
	if (is.na(LSL) & is.na(USL)) {
		stop ("No specification limits provided")
	}
	zz.m <- mean(x, na.rm = f.na.rm)
	zz.s <- sd(x, na.rm = f.na.rm)
	zul <- (USL - zz.m) / zz.s
	zll <- (zz.m - LSL) / zz.s
	
	if (is.na(zul)){
		z <- zll
	}
		else if (is.na(zll)){
			z <- zul
		}
		else {
			z <- min(zul, zll)
		}
		if (LT == FALSE){
			z <- z - 1.5
		} 
		return(as.vector(z))	
}


ss.ca.cp <- function(x, LSL = NA, USL = NA, 
		LT = FALSE, f.na.rm = TRUE, 
		ci = FALSE, alpha = 0.05){
	if (is.na(LSL) & is.na(USL)) {
		stop("No specification limits provided")
	}
	if (!is.numeric(x)){
		stop("Incorrect vector data")
	}
	cp.m <- mean(x, na.rm = f.na.rm)
	cp.s <- sd(x, na.rm = f.na.rm)
	cp.l <- (cp.m - LSL) / (3 * cp.s)
	cp.u <- (USL - cp.m) / (3 * cp.s)
	cp <- (USL - LSL) / (6 * cp.s)
	if (is.na(cp)){
		cp <- max(cp.l, cp.u, na.rm = TRUE)
	}
	if (ci == FALSE){
		return(as.numeric(cp))
	}
	else{
		return(c(
				as.numeric(cp)*sqrt((qchisq (alpha/2,length(x)-1,
											lower.tail=TRUE)/(length(x)-1))),				
				as.numeric(cp)*sqrt((qchisq (alpha/2,length(x)-1,
											lower.tail=FALSE)/(length(x)-1)))
		))							
	}
}

ss.ca.cpk <- function(x, LSL = NA, USL = NA, 
		LT = FALSE, f.na.rm = TRUE, 
		ci = FALSE, alpha = 0.05 ){
	if (is.na(LSL) & is.na(USL)) {
		stop("No specification limits provided")
	}
	if (!is.numeric(x)){
		stop("Incorrect vector data")
	}
	ss.n <- length(x[!is.na(x)])
	cpk.m <- mean(x, na.rm = f.na.rm)
	cpk.s <- sd(x, na.rm = f.na.rm)
	cpk.ul <- (USL - cpk.m) / (3 * cpk.s)
	cpk.ll <- (cpk.m - LSL) / (3 * cpk.s)
	cpk <- min(cpk.ul, cpk.ll, na.rm = TRUE)
	if (ci == FALSE){
		return(as.numeric(cpk))
	}
	else{
		return(c(
		cpk*(1-(qnorm(1-(alpha/2))*sqrt((1/(9*ss.n*cpk^2))+(1/(2*(ss.n-1)))))),
		cpk*(1+(qnorm(1-(alpha/2))*sqrt((1/(9*ss.n*cpk^2))+(1/(2*(ss.n-1))))))
		))
	}
		

}	

###############################################################################
#' Graphs and figures for a Capability Study
#' 
#' Plots a Histogram with density lines about the data of a process. Check normality
#' with qqplot and normality tests. Shows the Specification Limits and the 
#' Capability Indices.
#' 
#' @usage 
#' ss.study.ca(xST, xLT = NA, LSL = NA, USL = NA, Target = NA, 
#'   alpha = 0.05, f.na.rm = TRUE, f.main = "Six Sigma Capability Analysis Study", 
#'   f.sub = "")
#' @param xST Short Term process performance data 
#' @param xLT Long Term process performance data 
#' @param LSL Lower Specification Limit of the process
#' @param USL Upper Specification Limit of the process
#' @param Target Target of the process
#' @param alpha Type I error for the Confidence Interval
#' @param f.na.rm If TRUE NA data will be removed
#' @param f.main Main Title for the graphic output
#' @param f.sub Subtitle for the graphic output
#' @return Figures and plot for Capability Analysis
#' 
#' @references 
#' Cano, Emilio L., Moguerza, Javier M. and Redchuk, Andres. 2012.
#' \emph{Six Sigma with {R}. Statistical Engineering for Process
#'   Improvement}, Use R!, vol. 36. Springer, New York.
#'   \url{http://www.springer.com/statistics/book/978-1-4614-3651-5}.
#'   
#' Montgomery, DC (2008) \emph{Introduction to Statistical Quality Control}
#'   (Sixth Edition). New York: Wiley&Sons
#' 
#' @seealso \code{\link{ss.ca.cp}}
#' @author EL Cano
#' @export
#' @examples 
#' 	ss.study.ca(ss.data.ca$Volume, rnorm(40, 753, 3), 
#' 		LSL = 740, USL = 760, T = 750, alpha = 0.05, 
#'  			f.sub = "Winery Project")
ss.study.ca<-function (xST, xLT = NA, LSL = NA, USL = NA, 
		Target = NA, alpha = 0.05, 
		f.na.rm = TRUE,
		f.main = "Six Sigma Capability Analysis Study", 
		f.sub = ""){
	if (is.na(Target)){
		stop("Target is needed")
	}
	if (is.na(LSL) & is.na(USL)){
		stop("No specification limits provided")
	}
	#Facts
	mST = mean(xST, na.rm = f.na.rm)
	sST = sd(xST, na.rm = f.na.rm)
	nST = length(xST[!is.na(xST)])
	nLT = length(xLT[!is.na(xLT)])
	zST = ss.ca.z(xST, LSL, USL)
	cpST = ss.ca.cp(xST, LSL, USL)
	cpiST = ss.ca.cp(xST, LSL, USL, ci = TRUE, alpha = alpha)
	cpkST = ss.ca.cpk(xST, LSL, USL)
	cpkiST = ss.ca.cpk(xST, LSL, USL, ci = TRUE, alpha = alpha)
	DPMO <- (1 - pnorm(zST - 1.5)) * 10^6
	if (is.numeric(xLT)){
		mLT = mean(xLT, na.rm = f.na.rm)
		sLT = sd(xLT,na.rm = f.na.rm)
		cpLT = ss.ca.cp(xLT, LSL, USL, LT = TRUE)	
		cpiLT = ss.ca.cp(xLT, LSL, USL, LT = TRUE, ci = TRUE, alpha = alpha)
		cpkLT = ss.ca.cpk(xLT, LSL, USL, LT = TRUE)
		cpkiLT = ss.ca.cpk(xLT, LSL, USL, LT = TRUE, ci = TRUE, alpha = alpha)
		zLT = ss.ca.z(xLT, LSL, USL, LT = TRUE)
		DPMO <- (1 - pnorm(zLT)) * 10^6
	}
	else{
		mLT=NA
		sLT=NA
		cpLT=NA	
		cpiLT=NA
		cpkLT=NA
		cpkiLT=NA
		zLT<-zST-1.5
	}

######
	.ss.prepCanvas(f.main, f.sub)
#grid::grid.rect()##########
	vp.plots<-grid::viewport(name="plots",
			layout=grid::grid.layout(2,2,c(0.6,0.4),c(0.6,0.4)))
	grid::pushViewport(vp.plots)

	vp.hist <- grid::viewport(name="hist", layout.pos.row=1, layout.pos.col=1)
	grid::pushViewport(vp.hist)
#grid::grid.rect()##########
	grid::grid.text("Histogram & Density", y=1, just=c("center", "top") )

##############	

binwST <- diff(range(xST))/ sqrt(nST)
ggdata <- reshape2::melt(xST)
qqp <- ggplot(ggdata, aes(x=value))
hist <- qqp + geom_histogram(aes(y = ..density..), 
				binwidth = binwST,
				fill = "steelblue", 
				stat = "bin")
xST_density <- density(xST, bw = binwST)
if (!is.na(LSL)){
	hist <- hist +
		annotate(geom = "text", 
				x = LSL, 
				y = max(xST_density$y), 
				label = "LSL", 
				hjust = -0.1, 
				size = 5) 
} 
hist <- hist +	annotate(geom = "text",
				x = Target, 
				y = max(xST_density$y), 
				label = "Target",
				hjust = -0.1,
				size = 5)
if (!is.na(USL)){
	hist <- hist + 
		annotate(geom = "text",
				x = USL, 
				y = max(xST_density$y), 
				label = "USL",
				hjust = 1.1, 
				size = 5) 
}
	hist <- hist + xlab(NULL) + 
		ylab(NULL) + 
		theme(axis.text.y = element_blank())
if (!is.na(LSL)){
		hist <- hist + geom_vline(xintercept = LSL,
				linetype = 2,
				size = 1) 
	}
if (!is.na(USL)){
	hist <- hist + geom_vline(xintercept = USL,
			linetype = 2,
			size = 1) 
}
	hist <- hist + geom_vline(xintercept = Target,
				linetype = 3, 
				size = 1) +
		stat_density(geom="path", 
				position="identity", 
				size = 1) +
		stat_function( 
				fun = dnorm, 
				args = with(ggdata,	c(mean(value), sd(value))),
				linetype = 2, 
				size = 1
		) 

if (is.numeric(xLT)){
	binwLT <- diff(range(xLT))/ sqrt(nLT)
	ggdataLT <- reshape2::melt(xLT)
	hist <- hist + 
		stat_density(geom="path",
				data = ggdataLT, 
				position = "identity") + 
		stat_function(
				fun = dnorm, 
				args = with(ggdataLT, 
						c(mean = mean(value), sd = sd(value))),
				linetype=2
		)
} 

	print(hist, newpage=FALSE)
	
	grid::popViewport()
	vp.norm<-grid::viewport(name="normal",layout.pos.row=2, layout.pos.col=1,
			layout=grid::grid.layout(2,2,c(0.6,0.4),c(0.1, 0.9)))
	grid::pushViewport(vp.norm)
	grid::grid.text("Check Normality", y=1,just=c("center","top"))
#grid::grid.rect()##########
	vp.qq<-grid::viewport(name="qqp",layout.pos.row=2,layout.pos.col=1, 
			height=unit(0.5,"npc"))
	grid::pushViewport(vp.qq)
#grid::grid.rect()##########

	qqp <- qplot(sample = xST) + 
			xlab(NULL) + ylab(NULL) +
			theme(axis.text.x = element_blank(), 
              axis.text.y = element_blank()) 
	print(qqp,newpage=FALSE)
	grid::popViewport()
	vp.testn<-grid::viewport(name="testn",layout.pos.row=2, layout.pos.col=2)
	grid::pushViewport(vp.testn)
	ss.ts <- shapiro.test(xST)
	ss.tl <- nortest::lillie.test(xST)
	if (min(ss.ts$p.value, ss.tl$pvalue) < alpha){
		warning("Normality test/s failed")
	} 
	grid::grid.text("Shapiro-Wilk Test", y=.9,just=c("center","top"), 
			gp=grid::gpar(cex=.8))
	grid::grid.text(paste("p-value: ",format(ss.ts$p.value,digits=4)),
			gp=grid::gpar(cex=.8), y=.8)
	grid::grid.text("Lilliefors (K-S) Test", gp=grid::gpar(cex=.8))
	grid::grid.text(paste("p-value: ", format(ss.tl$p.value,digits=4)),
			gp=grid::gpar(cex=.8),y=.4)
	grid::popViewport()
	grid::grid.text("Normality accepted when p-value > 0.05",y=0.02, 
			just=c("center","bottom"), gp=grid::gpar(cex=.8))
	grid::popViewport()
	vpNumbers<-grid::viewport(name="numbers", 
			layout.pos.row=c(1,2), layout.pos.col=2,
			layout=grid::grid.layout(4,1))
	grid::pushViewport(vpNumbers)
	
grid::grid.rect(gp=grid::gpar(col="#BBBBBB",lwd=2))##########
	vpLegend<-grid::viewport(name="legend", layout.pos.row=1)
	grid::pushViewport(vpLegend)

grid::grid.rect(gp=grid::gpar(col="#BBBBBB",lwd=2))##########
	grid::grid.text(expression(bold("Density Lines Legend")), 
			y=0.95, just=c("center","top"))
	grid::grid.lines(x=c(0.05,0.3), y=c(0.7,0.7), gp=grid::gpar(lty=1, lwd=3))
	grid::grid.text("Density ST", x=0.35, y=0.7,just=c("left","center"),
			gp=grid::gpar(cex=0.8))
	
	grid::grid.lines(x=c(0.05,0.3), y=c(0.55,0.55), gp=grid::gpar(lty=2, lwd=3))
	grid::grid.text("Theoretical Dens. ST", x=0.35, y=0.55,just=c("left","center"), 
			gp=grid::gpar(cex=0.8))

if (is.numeric(xLT)){	
	grid::grid.lines(x=c(0.05,0.3), y=c(0.40,0.40), gp=grid::gpar(lty=1, lwd=1))
	grid::grid.text("Density LT", x=0.35, y=0.40,just=c("left","center"), 
			gp=grid::gpar(cex=0.8))
	
	grid::grid.lines(x=c(0.05,0.3), y=c(0.25,0.25), gp=grid::gpar(lty=2, lwd=1))
	grid::grid.text("Theoretical Density LT", x=0.35, y=0.25,just=c("left","center"), 
			gp=grid::gpar(cex=0.8))
}	
	grid::popViewport()
	vpSpec<-grid::viewport(name="spec", layout.pos.row=2)
	grid::pushViewport(vpSpec)
#grid::grid.rect()#############
	grid::grid.text(expression(bold("Specifications")), y=.95, just=c("center","top"))
	grid::grid.text(expression(bold("LSL: ")), 
			y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(LSL, y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("Target: ")), 
			y=unit(.95,"npc")-unit(2.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(Target, y=unit(.95,"npc")-unit(2.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("USL: ")), 
			y=unit(.95,"npc")-unit(3.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(USL, y=unit(.95,"npc")-unit(3.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::popViewport()
	vpProcess<-grid::viewport(name="proc", layout.pos.row=3,
			layout=grid::grid.layout(1,2))
	grid::pushViewport(vpProcess)
#grid::grid.rect()##############
	grid::grid.lines(x=c(0,1),y=c(1,1), gp=grid::gpar(col="#BBBBBB", lwd=3))
	grid::grid.text(expression(bold("Process")), y=.95, just=c("center","top"))
	vpSTp<-grid::viewport(layout.pos.col=1)
	grid::pushViewport(vpSTp)
	grid::grid.text("Short Term",x=0.05, y=.95, just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("Mean: ")), y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.4f",mST), y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("SD: ")), y=unit(.95,"npc")-unit(2.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.4f",sST), y=unit(.95,"npc")-unit(2.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("n: ")), y=unit(.95,"npc")-unit(3.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(nST, y=unit(.95,"npc")-unit(3.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold(Z[s]*": ")), y=unit(.95,"npc")-unit(4.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.2f",zST), y=unit(.95,"npc")-unit(4.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	
	grid::popViewport()
	vpLTp<-grid::viewport(layout.pos.col=2)
	grid::pushViewport(vpLTp)
	grid::grid.text("Long Term",x=.95, y=.95, just=c("right","top"), gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("Mean: ")), y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.4f",mLT), y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("SD: ")), y=unit(.95,"npc")-unit(2.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.4f",sLT), y=unit(.95,"npc")-unit(2.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("n: ")), y=unit(.95,"npc")-unit(3.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(nLT, y=unit(.95,"npc")-unit(3.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold(Z[s]*": ")), y=unit(.95,"npc")-unit(4.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.2f",zLT), y=unit(.95,"npc")-unit(4.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("DPMO: ")), y=unit(.95,"npc")-unit(5.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(round(DPMO,1), y=unit(.95,"npc")-unit(5.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::popViewport()
	grid::popViewport()
	vpIndices<-grid::viewport(name="ind", layout.pos.row=4,
			layout=grid::grid.layout(1,2))
	grid::pushViewport(vpIndices)
#grid::grid.rect()###############
grid::grid.lines(x=c(0,1), y=c(1,1), gp=grid::gpar(col="#BBBBBB",lwd=2))
	grid::grid.text(expression(bold("Indices")),y=.95,just=c("center","top"))
	vpSTi<-grid::viewport(layout.pos.col=1)
	grid::pushViewport(vpSTi)
	grid::grid.text("Short Term",x=0.05, y=.95, just=c("left","top"),gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold(C[p]*": ")), y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.4f",cpST), y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("CI: ")), y=unit(.95,"npc")-unit(3,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.7))
	grid::grid.text(paste("[",paste(sprintf("%.1f",cpiST[1]),sep=""),
					",",sprintf("%.1f",cpiST[2]),"]",sep=""), 
					y=unit(.95,"npc")-unit(3,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.7))
	grid::grid.text(expression(bold(C[pk]*": ")), y=unit(.95,"npc")-unit(4.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.4f",cpkST), y=unit(.95,"npc")-unit(4.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("CI: ")), y=unit(.95,"npc")-unit(6.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.7))
	grid::grid.text(paste("[",paste(sprintf("%.1f",cpkiST[1]),sep=""),
					",",sprintf("%.1f",cpkiST[2]),"]",sep=""), 
			y=unit(.95,"npc")-unit(6.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.7))

	grid::popViewport()
	vpLTi<-grid::viewport(layout.pos.col=2)
	grid::pushViewport(vpLTi)
	grid::grid.text("Long Term",x=.95, y=.95, just=c("right","top"), gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold(P[p]*": ")), y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.4f",cpLT), y=unit(.95,"npc")-unit(1.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("CI: ")), y=unit(.95,"npc")-unit(3,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.7))
	grid::grid.text(paste("[",paste(sprintf("%.1f",cpiLT[1]),sep=""),
					",",sprintf("%.1f",cpiLT[2]),"]",sep=""), 
			y=unit(.95,"npc")-unit(3,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.7))
	grid::grid.text(expression(bold(P[pk]*": ")), y=unit(.95,"npc")-unit(4.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(sprintf("%.4f",cpkLT), y=unit(.95,"npc")-unit(4.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.8))
	grid::grid.text(expression(bold("CI: ")), y=unit(.95,"npc")-unit(6.5,"lines"), 
			just=c("right","top"),
			gp=grid::gpar(cex=.7))
	grid::grid.text(paste("[",paste(sprintf("%.1f",cpkiLT[1]),sep=""),
					",",sprintf("%.1f",cpkiLT[2]),"]",sep=""), 
			y=unit(.95,"npc")-unit(6.5,"lines"), 
			just=c("left","top"),
			gp=grid::gpar(cex=.7))
	grid::popViewport()
	grid::popViewport()
	
}
#trellis.par.set(standard.theme(color=FALSE))
#ss.study.ca(x,LSL=740, USL=760, T=750, f.sub="Winery Project")
