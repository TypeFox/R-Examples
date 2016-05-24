#' Ternary diagram
#' 
#' This plot shows the relative proportions of three variables (compositional
#' parts) in one diagramm.  Before plotting, the data are scaled.
#' 
#' The relative proportions of each variable are plotted.
#' 
#' @param x matrix or data.frame with 3 columns
#' @param name names of the variables
#' @param grid if TRUE a grid is plotted additionally in the ternary diagram
#' @param gridCol color for the grid lines
#' @param mcex label size
#' @param line may be set to \dQuote{none}, \dQuote{pca}, \dQuote{regression},
#' \dQuote{regressionconf}, \dQuote{regressionpred}, \dQuote{ellipse},
#' \dQuote{lda}
#' @param robust if line equals TRUE, it dedicates if a robust estimation is
#' applied or not.
#' @param group if line equals \dQuote{da}, it determines the grouping variable
#' @param tol if line equals \dQuote{ellipse}, it determines the parameter for
#' the tolerance ellipse
#' @param \dots further parameters, see, e.g., \code{par()}
#' @author Peter Filzmoser <\email{P.Filzmoser@@tuwien.ac.at}>
#' \url{http://www.statistik.tuwien.ac.at/public/filz/}, Matthias Templ
#' @seealso \code{\link[StatDA]{ternary}}
#' @references C. Reimann, P. Filzmoser, R.G. Garrett, and R. Dutter:
#' Statistical Data Analysis Explained. Applied Environmental Statistics with
#' R. John Wiley and Sons, Chichester, 2008.
#' @keywords multivariate aplot
#' @export
#' @examples
#' 
#' data(arcticLake)
#' ternaryDiag(arcticLake)
#' 
#' data(coffee)
#' x <- coffee[,2:4]
#' grp <- as.integer(coffee[,1])
#' ternaryDiag(x, col=grp, pch=grp)
#' ternaryDiag(x, grid=FALSE, col=grp, pch=grp)
#' legend("topright", legend=unique(coffee[,4]), pch=1:2, col=1:2)
#' 
#' ternaryDiag(x, grid=FALSE, col=grp, pch=grp, line="ellipse", tol=c(0.975,0.9), lty=2)
#' ternaryDiag(x, grid=FALSE, line="pca")
#' ternaryDiag(x, grid=FALSE, col=grp, pch=grp, line="pca", lty=2, lwd=2)
#' 
ternaryDiag <- function(x, name=colnames(x), grid=TRUE, 
		                gridCol=grey(0.6), mcex=1.2, line="none", 
						robust=TRUE, group=NULL, tol=0.975, ...)
{
# Ternary plot
#
# x ... matrix with 3 columns
# name ... names of the variables
# grid ... TRUE if grid should be plotted
# "..." ... further graphical parameters, see par
	
	if(!(line %in% c("none","pca","regression","regressionconf","regressionpred","ellipse","lda"))) stop("Choose a method for function argument *line*, which is implemented.")
	if (is.null(name)) { name <- c("P1", "P2", "P3")}
	if (length(name) != 3){ 
		warning("incorrect length of name. Variable names P1, P2 and P3 are used instead.")
	}
	if (dim(x)[2] > 3){ 
		warning("only the first three parts are used for plotting")
		x <- x[,1:3]
	}
	if (dim(x)[2] < 3){ 
		stop("x must include 3 variables/parts")
	}
	s <- rowSums(x)
	if (any(s <= 0))
		stop("each row of the input `object' must have a positive sum")
	dat <- x/s 
#	dat <- constSum(x)
	
	xp <- dat[,2] + dat[,3]/2
	yp <- dat[,3] * sqrt(3)/2
	
	par(pty="s")
	plot(xp,yp,xlim=c(0,1),ylim=c(0,0.9), 
			frame.plot=FALSE, xaxt="n", yaxt="n", xlab="", ylab="", ...)
	
	segments(0,0,1,0)
	segments(0,0,1/2,sqrt(3)/2)
	segments(1/2,sqrt(3)/2,1,0)
	
	mtext(name[1],side=1, line=-1, at=-0.05,cex=mcex)
	mtext(name[2],side=1, line=-1, at=1.05,cex=mcex)
	text(0.5, 0.9, name[3],cex=mcex)
	
	if(grid)
	{    
		b <- sqrt(c(0.03,0.12,0.27,0.48))
#		segments(c(0.2,0.4,0.6,0.8,0.2,0.4,0.6,0.8,0.1,0.2,0.3,0.4), 
#				 c(rep(0,8), b),
#		         c(seq(0.1,0.9,0.1), c(0.9,0.8,0.7,0.6)),
#				 c(b, rev(b), b), col=gridCol, lty="dashed")
#				 
#		 )
		segments(0.2,0, 0.1,sqrt(0.03), col=gridCol, lty="dashed")
		segments(0.4,0, 0.2,sqrt(0.12), col=gridCol, lty="dashed")
		segments(0.6,0, 0.3,sqrt(0.27), col=gridCol, lty="dashed")
		segments(0.8,0, 0.4,sqrt(0.48), col=gridCol, lty="dashed")
		segments(0.2,0, 0.6,sqrt(0.48), col=gridCol, lty="dashed")
		segments(0.4,0, 0.7,sqrt(0.27), col=gridCol, lty="dashed")
		segments(0.6,0, 0.8,sqrt(0.12), col=gridCol, lty="dashed")
		segments(0.8,0, 0.9,sqrt(0.03), col=gridCol, lty="dashed")
		segments(0.1,sqrt(0.03), 0.9,sqrt(0.03), col=gridCol, lty="dashed")
		segments(0.2,sqrt(0.12), 0.8,sqrt(0.12), col=gridCol, lty="dashed")
		segments(0.3,sqrt(0.27), 0.7,sqrt(0.27), col=gridCol, lty="dashed")
		segments(0.4,sqrt(0.48), 0.6,sqrt(0.48), col=gridCol, lty="dashed")
		
		text(0.5,0.66,"0.8", col=gridCol, cex = 0.6)
		text(0.5,0.49,"0.6", col=gridCol, cex = 0.6)
		text(0.5,0.32,"0.4", col=gridCol, cex = 0.6)
		text(0.5,0.14,"0.2", col=gridCol, cex = 0.6)
		text(0.95,0.21,"0.8", col=gridCol, cex = 0.6, srt = 60)
		text(0.86,0.35,"0.6", col=gridCol, cex = 0.6, srt = 60)
		text(0.75,0.54,"0.4", col=gridCol, cex = 0.6, srt = 60)
		text(0.64,0.72,"0.2", col=gridCol, cex = 0.6, srt = 60)
		text(0.05,0.21,"0.8", col=gridCol, cex = 0.6,srt = 300)
		text(0.14,0.35,"0.6", col=gridCol, cex = 0.6,srt = 300)
		text(0.25,0.54,"0.4", col=gridCol, cex = 0.6,srt = 300)
		text(0.36,0.72,"0.2", col=gridCol, cex = 0.6,srt = 300)
	}

  ## adding lines
  plotTern <- function(x, conf=line, rob=robust){
#	  x <- x[,c(3,2,1)]
	  z <- data.frame(isomLR(x))
#	  z[,2] <- z[,2]*(-1)
	  colnames(z) <- c("x", "y")
	  if(rob) lm1 <- MASS::rlm(y ~ x, data=z, method="MM") else lm1 <- lm(y ~ x, data=z)
#	  new <- data.frame(x = seq(min(z$x), max(z$x), length=100)) 
	  new <- data.frame(x = seq(-30, 30, length=10000)) 
	  if(conf=="regressionpred"){
	    pred.w.plim <- predict(lm1, new, interval="prediction")
		s1 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim[,1]))
		s2 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim[,2]))
		s3 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim[,3]))
     } else if(conf=="regressionconf"){
	    pred.w.plim <- predict(lm1, new, interval="confidence")		
		s1 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim[,1]))
		s2 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim[,2]))
		s3 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim[,3]))
	  } else {
		pred.w.plim <- predict(lm1, new)	
		s1 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim))
		s2 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim))
		s3 <- isomLRinv(data.frame(z1=new$x, z2=pred.w.plim))
	  }
#	  s1 <- invilr(data.frame(z1=new$x, z2=pred.w.plim[,1]))
#	  s2 <- invilr(data.frame(z1=new$x, z2=pred.w.plim[,2]))
##	s2 <- invilr(data.frame(z1=new$x, z2=pred.w.clim[,2]))	
#	  s3 <- invilr(data.frame(z1=new$x, z2=pred.w.plim[,3]))
#	s3 <- invilr(data.frame(z1=new$x, z2=pred.w.clim[,3]))
	  
#	  s1r <- rowSums(s1)
#	  s2r <- rowSums(s2)
#	  s3r <- rowSums(s3)
	  dat1 <- s1#/s1r
	  dat2 <- s2#/s2r
	  dat3 <- s3#/s3r
	  
	  xp1 <- dat1[, 2] + dat1[, 3]/2
	  yp1 <- dat1[, 3] * sqrt(3)/2
	  xp2 <- dat2[, 2] + dat2[, 3]/2
	  yp2 <- dat2[, 3] * sqrt(3)/2
	  xp3 <- dat3[, 2] + dat3[, 3]/2
	  yp3 <- dat3[, 3] * sqrt(3)/2
	  lines(xp1, yp1, xlim = c(0, 1), ylim = c(0, 0.9), #frame.plot = FALSE, 
			  xaxt = "n", yaxt = "n", xlab = "", ylab = "", col=1,
			  lwd=2)
      if(conf %in% c("regressionpred", "regressionconf")){	  
		 lines(xp2, yp2, xlim = c(0, 1), ylim = c(0, 0.9), #frame.plot = FALSE, 
			  xaxt = "n", yaxt = "n", xlab = "", ylab = "", lty=2, col=gray(0.5))
	     lines(xp3, yp3, xlim = c(0, 1), ylim = c(0, 0.9), #frame.plot = FALSE, 
			  xaxt = "n", yaxt = "n", xlab = "", ylab = "", lty=2, col=gray(0.5))
      }
	  legend("topleft", legend=paste(colnames(x)[1], "~", colnames(x)[2], "+", colnames(x)[3]), lwd=2, col="black")
  }
  
  f <- function(x, co="black", lt=1, rob=robust){
	  a <- constSum(x,1)
	  a <- isomLR(x)
	  if(rob){
		  rc <- robustbase::covMcd(a)
		  me <- rc$center
		  acov <- rc$cov
	  } else {
	      me <- apply(a, 2, mean)
	      acov <- var(a)
      }
#	  cov.svd <- svd(acov, nv = 0)
#	  r <- cov.svd[["u"]] %*% diag(sqrt(cov.svd[["d"]]))
	  apca <- princomp(cov=acov)
	  x0 <- -11*apca$loa[1,1]+me[1]
	  y0 <- -11*apca$loa[2,1]+me[2]
	  x1 <- 11*apca$loa[1,1]+me[1]
	  y1 <- 11*apca$loa[2,1]+me[2]
	  
	  s1 <- seq(x0,x1,length=100)
	  s2 <- seq(y0,y1,length=100)
	  
	  s <- data.frame(s1=s1,s2=s2)
	  ss <- isomLRinv(s)
#	ternaryDiag(arcticLake)
	  s1 <- rowSums(ss)
	  dat <- ss/s1
	  xp <- dat[, 2] + dat[, 3]/2
	  yp <- dat[, 3] * sqrt(3)/2
	  lines(xp, yp, xlim = c(0, 1), ylim = c(0, 0.9), #frame.plot = FALSE, 
			  xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...)
	  legend("right", legend="PC 1", lty=1, lwd=2, col="black")
  }
  ### tollerance ellipses:
  dcov <- function(x, tolerance=tol){
	  z <- isomLR(x)
	  dat1 <- drawMahal(z, colMeans(z), cov(z), plot=FALSE, whichlines=tolerance) 
	  for(i in 1:length(tolerance)){
	      e <- isomLRinv(cbind(dat1$mdX[,i], dat1$mdY[,i]))
		  xp1 <- e[, 2] + e[, 3]/2
		  yp1 <- e[, 3] * sqrt(3)/2	  
		  lines(xp1, yp1, xlim = c(0, 1), ylim = c(0, 0.9), #frame.plot = FALSE, 
				  xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...)
      }
  }
  da <- function(x, grp=group){
	  z <- isomLR(x)
	  lev <- levels(factor(grp))
	  if(length(lev) != 2) stop("group must be a factor with exactly two levels")
	  z1 <- z[grp==lev[1],]
	  z2 <- z[grp==lev[2],]
	  n1=nrow(z1)
	  n2=nrow(z2)
	  n=n1+n2
	  p1=n1/n
	  p2=n2/n
	  m1=apply(z1,2,mean)
	  m2=apply(z2,2,mean)
	  S1=cov(z1)
	  S2=cov(z2)
	  Sp=((n1-1)/(n1-1+n2-1))*S1+((n2-1)/(n1-1+n2-1))*S2
	  Sp1=solve(Sp)
	  yLDA=as.numeric(t(m1-m2)%*%Sp1%*%t(z)-as.numeric(1/2*t(m1-m2)%*%Sp1%*%(m1+m2)))-log(p2/p1)
# LDA
	  y1=seq(from=min(z[,1])-1.5,to=max(z[,1])+1.9,by=0.05)
	  y2=seq(from=min(z[,2]),to=max(z[,2])+0.2,by=0.05)
	  y1a=rep(y1,length(y2))
	  y2a=sort(rep(y2,length(y1)))
	  ya=cbind(y1a,y2a)
	  
	  yaLDA <- as.numeric(t(m1-m2)%*%Sp1%*%t(ya)-
					  as.numeric(1/2*t(m1-m2)%*%Sp1%*%(m1+m2)))-log(p2/p1)
	  
	  boundLDA <- abs(yaLDA)<1.5
	  bline <- lowess(y1a[boundLDA],y2a[boundLDA])
#	  bline <- lowess(y1a,y2a)
	  blines <- data.frame(z1=bline$x, z2=bline$y)
	  #k*x+d
      k <- (bline$x[2]-bline$x[1]) / (bline$y[2]-bline$y[1])
	  LINE <- function(p,k){
		  seq(p,0.95,) ## stopped here!
	  }
#	  blines <- data.frame(z1=seq(bline$x[1],bline$x[2],length=100), z2=seq(bline$y[1],bline$y[2],length=100))
	  xblines <- isomLRinv(blines)
	  xp1 <- xblines[, 2] + xblines[, 3]/2
	  yp1 <- xblines[, 3] * sqrt(3)/2	  
	  lines(xp1, yp1, xlim = c(0, 1), ylim = c(0, 0.9), #frame.plot = FALSE, 
			  xaxt = "n", yaxt = "n", xlab = "", ylab = "", ...)
  }
  
  
  if(line == "pca") f(dat, co="black", lt=1)
  if(line == "regression") plotTern(dat[,c(1,2,3)], conf=FALSE)
  if(line == "regressionconf") plotTern(dat[,c(1,2,3)], conf="regressionconf")
  if(line == "regressionpred") plotTern(dat[,c(1,2,3)], conf="regressionpred")
  if(line == "ellipse") dcov(x)
  if(line == "lda"){ 
#	  if(is.null(group)) stop("parameter group must be specified.")
	  da(x, grp=group)
	  grp <- group
	  points(xp[grp==levels(factor(grp))[1]],yp[grp==levels(factor(grp))[1]], pch=1)
	  points(xp[grp==levels(factor(grp))[2]],yp[grp==levels(factor(grp))[2]], pch=4)	  
  }
}


