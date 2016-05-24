risk <- function(f, g, delta = 0, log = TRUE, h = NULL, adaptive = FALSE, res = 50, WIN = NULL, tolerate = FALSE, plotit = TRUE, comment = TRUE){
    
	
	if(class(f)!=class(g)) stop("both 'f' and 'g' must have the same class")
	
	pooled <- NA
	org.cls <- class(f)
	if(comment&&org.cls!="bivden") cat(date())
	
	if(class(f)!="bivden"){
		if(class(f)=="data.frame"){
			if(ncol(f)!=2||ncol(g)!=2) stop("data.frames must both have exactly two columns")
			if(nrow(f)<10) warning("less than 10 case observations!")
			if(nrow(g)<10) warning("less than 10 control observations!")
		} else if(class(f)=="list"){
			if(is.null(f$x)||is.null(f$y)||is.null(g$x)||is.null(g$y)) stop("data lists must both have two components named 'x' and 'y'")
			if(length(f$x)!=length(f$y)) stop("data components 'x' and 'y' are of unequal lengths in 'f'")
			if(length(g$x)!=length(g$y)) stop("data components 'x' and 'y' are of unequal lengths in 'g'")
			if(length(f$x)<10) warning("less than 10 case observations!")
			if(length(g$x)<10) warning("less than 10 control observations!")
			f <- as.data.frame(f)
			g <- as.data.frame(g)
		} else if(class(f)=="matrix"){
			if(ncol(f)!=2||ncol(g)!=2) stop("data matrices must both have exactly two columns")
			if(nrow(f)<10) warning("less than 10 case observations!")
			if(nrow(g)<10) warning("less than 10 control observations!")
			f <- as.data.frame(f)
			g <- as.data.frame(g)
		} else if(class(f)=="ppp"){
			if(is.null(f$x)||is.null(f$y)||is.null(g$x)||is.null(g$y)) stop("data ppp.objects must both have two non-empty components named 'x' and 'y'")
			if(length(f$x)<10) warning("less than 10 case observations!")
			if(length(g$x)<10) warning("less than 10 control observations!")
			if(!identical_windows(f$window,g$window)) stop("'f' and 'g' ppp.objects must have identical 'window' components")
			WIN <- f$window
			f <- data.frame(cbind(f$x,f$y))
			g <- data.frame(cbind(g$x,g$y))
		} else {
			stop("'f' and 'g' must be objects of type 'data.frame', 'list', 'matrix', 'ppp', or 'bivden'")
		}
		
		if(!is.null(h)){
			if(!is.vector(h)){
				stop("'h' must be a numeric vector of length 1 or 2")
			} else {
				if(length(h)!=1&&length(h)!=2){
					stop("'h' must be a numeric vector of length 1 or 2")
				}
			}
		}
		
		if(is.null(WIN)&&class(f)!="ppp") stop("'WIN' must be provided as an object of class 'owin' if 'f' and 'g' are not of class 'bivden' or 'ppp'")
		
		h0p <- OS(list(x=c(f[,1],g[,1]),y=c(f[,2],g[,2])))
		
		if(is.null(h)){
			h01 <- h02 <- OS(list(x=c(f[,1],g[,1]),y=c(f[,2],g[,2])),nstar=sqrt(nrow(f)*nrow(g)))
		} else {
			if(length(h)==1){
				h01 <- h02 <- h[1]
			} else {
				h01 <- h[1]
				h02 <- h[2]
			}
		}
		
		if(adaptive){
			hp <- 0.5*h0p
			hf <- 0.5*h01
			hg <- 0.5*h02
			
			if(comment) cat("\ncalculating pooled density...\n")
			pooled <- bivariate.density(data=list(x=c(f[,1],g[,1]),y=c(f[,2],g[,2])),pilotH=hp,globalH=h0p,adaptive=adaptive,res=res,WIN=WIN,comment=F)
			if(comment) cat("calculating case density...\n")
			f <- bivariate.density(data=f,pilotH=hf,globalH=h01,adaptive=adaptive,res=res,WIN=WIN,gamma=pooled$gamma,comment=F)
			if(comment) cat("calculating control density...\n")
			g <- bivariate.density(data=g,pilotH=hg,globalH=h02,adaptive=adaptive,res=res,WIN=WIN,gamma=pooled$gamma,comment=F)
		} else {
			if(comment) cat("\ncalculating pooled density...\n")
			pooled <- bivariate.density(data=list(x=c(f[,1],g[,1]),y=c(f[,2],g[,2])),pilotH=h0p,adaptive=adaptive,res=res,WIN=WIN,comment=F)
			if(comment) cat("calculating case density...\n")
			f <- bivariate.density(data=f,pilotH=h01,adaptive=adaptive,res=res,WIN=WIN,comment=F)
			if(comment) cat("calculating control density...\n")
			g <- bivariate.density(data=g,pilotH=h02,adaptive=adaptive,res=res,WIN=WIN,comment=F)
		}
	} 
	
	if(delta<0) delta <- 0
	fvec <- as.vector(t(f$Zm))
	gvec <- as.vector(t(g$Zm))
	if(length(fvec)!=length(gvec)) stop("'f' and 'g' densities must be estimated on grids with identical resolutions") 
	if(!identical_windows(f$WIN,g$WIN)) stop("study regions 'f$WIN' and 'g$WIN' appear different - vertices must match!")
	
	if(length(f$hypoH)!=length(g$hypoH)) stop("smoothing approaches of 'f' and 'g' do not match! must both be adaptive or fixed")
	
	edgef <- range(as.vector(f$qhz),na.rm=T)[1]!=range(as.vector(f$qhz),na.rm=T)[2]
	edgeg <- range(as.vector(g$qhz),na.rm=T)[1]!=range(as.vector(g$qhz),na.rm=T)[2]
	
	if(edgef+edgeg==1) stop("edge-correction is inconsistent. both densities must be either edge-corrected or not.")
	
	gamCon <- delta*max(gvec[!is.na(gvec)])
    
    rhohatVec <- (fvec+gamCon)/(gvec+gamCon)
    rhohatVec[rhohatVec==0] <- NA
    if(log) rhohatVec <- log(rhohatVec)
    rhohatM <- matrix(rhohatVec,length(f$X),length(f$Y),byrow=T)
    
	result <- list(rsM=rhohatM,f=f,g=g,log=log,pooled=pooled,P=NA)
	class(result) <- "rrs"
	
	if(tolerate&&org.cls!="bivden"){
		if(comment) cat("running asymptotics for tolerance contours...")
		
		if(res>=50)	reduce <- 50/res
		else reduce <- 1 
		
		tol <- tolerance(result,pooled,reduce=reduce,comment=F)
		result$P <- tol$P
		if(comment) cat("done.\n")
	}
	
	if(plotit){
		plot(result,display="heat",col=heat.colors(12)[12:1],main="",xlab="",ylab="")
		if(tolerate) contour(tol$X,tol$Y,tol$P,levels=0.05,add=T)
	}
	
	if(comment&&org.cls!="bivden") cat("\n\n")
	if(comment&&org.cls!="bivden") cat(date())
    return(result)
}
