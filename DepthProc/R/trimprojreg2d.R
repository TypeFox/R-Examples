#' @name trimProjReg2d
#' @title trimProjReg2d
#' @export
#' 
#' @description Computes projection trimmed regression in 2 dimensions.
#'
#' @param x Independent variable
#' @param y Dependent variable
#' @param alpha Percentage of trimmed observations
#'
#'  @author Zygmunt Zawadzki from Cracow University of Economics.
#'
#'
#' @examples
#' 
#'  #EXAMPLE 1
#'  data(pension)
#'  plot(pension)
#'  abline(lm(Reserves~Income,data = pension), lty = 3, lwd = 2) #lm
#'  abline(trimProjReg2d(pension[,1],pension[,2]), lwd = 2) #trimprojreg2d
#'  legend("bottomright", c("OLS","TrimLS"), lty = 1:2)
#'  
#'  #EXAMPLE 2
#'  data(under5.mort)
#'  data(inf.mort)
#'  data(maesles.imm)
#'  
#'  data2011=na.omit(cbind(under5.mort[,22],inf.mort[,22],maesles.imm[,22]))
#'  x<-data2011[,3]
#'  y<-data2011[,2]
#'  plot(x,y,cex=1.2, ylab="infant mortality rate per 1000 live birth", 
#'  xlab="against masles immunized #'  
#'  percentage", main='Projection Depth Trimmed vs. LS regressions')
#'  abline(lm(x~y,data = pension), lwd = 2, col='black') #lm
#'  abline(trimProjReg2d(x,y), lwd = 2,col='red') #trimmed reg
#'  legend("bottomleft",c("LS","TrimReg"),fill=c("black","red"),cex=1.4,bty="n")
#'  
trimProjReg2d<-function(x, y, alpha = 0.1)
{

	yX= cbind(y,x)
	depth = depth(yX,yX,method = "Projection")
	cut = quantile(depth,alpha)

	ycut = y[ depth > cut ]	
	xcut = x[ depth > cut ]
	
	fitcut = lm(ycut~xcut)$coeff 
  new("TrimReg2d",coef = fitcut)
  
}
