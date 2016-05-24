###Wrapper function.
###Current version are not able to control the nlm function, which will be improved in the later versions.

BMControl <- function(print=FALSE, zStepSize=1, logPiTol=3, returnParam=FALSE,
				CILevel=0.95,...)
{
  if(zStepSize<0|logPiTol<0|CILevel>1|CILevel<0)
      stop("Wrong Input !")
	return(list(print=print, 
			zStepSize=zStepSize,
			logPiTol=logPiTol,
			logPiSubTol=logPiTol*2,
			returnParam=returnParam,
			CILevel=CILevel,
			...))
}

BMAnimalTrack <- function(dataList, controlList=BMControl())
{
	
	if(controlList$print)	cat("Step 1. Estimating parameters.\n")
	#Parameter estimation.
	pest <- nlm(nllh.BB.Phi_XY,  c(0, 0), logS2=TRUE,
		gpsList=dataList$glist, print.level=as.numeric(controlList$print))
	
	#In the log scale again. Obtain the Hessian
	pestWh <- nlm(nllh.BB.Phi_XY, pest$est, hessian=TRUE, logS2=TRUE,
			gpsList=dataList$glist, print.level=0)
	
	if(controlList$print)	cat("Step 2. Searching for the grid of numerical Integration. \n")
	piMx <- zSearch(pestWh, dataList$glist, 
			controlList$zStepSize, controlList$logPiTol, controlList$logPiSubTol)
	if(controlList$print)	cat(nrow(piMx), "grid points are used. \n")
	
	
	if(controlList$print)	cat("Step 3. Bayesian Melding for bias correction and uncertainty characterization. \n")
	fbRes <- postMar.BB.Eta(piMx, dataList$XMx, dataList$glist, 
				printK=controlList$print)
	ciConst <- qnorm(1-(1-controlList$CILevel)/2)
	fbRes <- cbind(fbRes, fbRes[,1]-ciConst*sqrt(fbRes[,2]), fbRes[,1]+ciConst*sqrt(fbRes[,2]))
	colnames(fbRes[,3:4]) <- c("CI.lower", "CI.upper")
	
	if (controlList$returnParam)
	{
		return(list(param=exp(pest$est), etaMar=fbRes))
	}
	else
		return(fbRes)
}