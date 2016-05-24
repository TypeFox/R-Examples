plotPower <-
function(x,
											coefNum=1,
											include="All",
											yAxis=seq(from=0, to=100, by=20),
											xAxis=NULL,
											main=NULL,
											legendXY=NULL)
{
	## Checks
	##
	if(!is.element(class(x), c("tpsPower", "ccPower")))
		stop("Error: 'x' is neither a 'tpsPower' or 'ccPower' object")
	##
  if(!is.element(coefNum, 1:length(x$betaTruth)))
  	stop("Error: 'coefNum' is invalid")
  
  ##
  if(!is.element(include, c("All", "TPS", "WL", "PL", "ML", "CC")))
  	stop("Error: 'include' is invalid")
  
  ##
  if(is.null(xAxis))
  {
  	if(class(x) == "ccPower") xAxis <- x$nCC
  	if(class(x) == "tpsPower") xAxis <- x$nII
  }
  if(is.null(main))
  	main <- paste("Power for", colnames(x$betaPower)[coefNum])
  	
  ##
  xlab <- ifelse(class(x) == "ccPower" | include == "CC",
  							"Case-control sample size, n",
  							"Phase II sample size, n")
  plot(range(xAxis), range(yAxis), xlab=xlab,  ylab="Power", main=main, type="n", axes=FALSE)
  axis(1, at=xAxis)
  axis(2, at=yAxis)
  
  ## Power for a ccPower object
  ##
  if(class(x) == "ccPower")
  {
  	powerCC <- x$betaPower[-1, coefNum]
    lines(x$nCC, powerCC, type="o", lty=1, lwd=2, pch=2)
  }

	## Power for a tpsPower object
	##
	if(class(x) == "tpsPower")
	{
  	nLvls   <- length(x$nII)
  	powerCC <- x$betaPower[(1 + (0*nLvls) + c(1:nLvls)), coefNum]
  	powerWL <- x$betaPower[(1 + (1*nLvls) + c(1:nLvls)), coefNum]
  	powerPL <- x$betaPower[(1 + (2*nLvls) + c(1:nLvls)), coefNum]
  	powerML <- x$betaPower[(1 + (3*nLvls) + c(1:nLvls)), coefNum]
    
    ## WL estimator
    if(is.element(include, c("All", "TPS", "WL")))
    	lines(x$nII, powerWL, type="o", lty=1, lwd=2, pch=2)
    
    ## PL estimator
    if(is.element(include, c("All", "TPS", "PL")))
    	lines(x$nII, powerPL, type="o", lty=switch(include, "All"=2, "TPS"=2, "PL"=1), lwd=2, pch=3)
    
    ## ML estimator
    if(is.element(include, c("All", "TPS", "ML")))
    	lines(x$nII, powerML, type="o", lty=switch(include, "All"=3, "TPS"=3, "ML"=1), lwd=2, pch=4)
    
    ## CC estimator
    if(is.element(include, c("All", "CC")))
    	lines(x$nII, powerCC, lty=switch(include, "All"=4, "CC"=1), lwd=2)    
    
    ## Legend
    ##
    if(include == "TPS")
    {
    	if(is.null(legendXY))
    		legendXY <- c(max(xAxis) - ((max(xAxis) - min(xAxis)) * 0.2), min(yAxis) + ((max(yAxis) - min(yAxis)) * 0.2))
    	legend(legendXY[1], legendXY[2], c("WL", "PL", "ML"), pch=c(2:4), lwd=c(2, 2, 2), lty=c(1:3))
    }
    if(include == "All")
    {
    	if(is.null(legendXY))
    		legendXY <- c(max(xAxis) - ((max(xAxis) - min(xAxis)) * 0.2), min(yAxis) + ((max(yAxis) - min(yAxis)) * 0.18))
      legend(legendXY[1], legendXY[2], c("WL", "PL", "ML", "CC"), pch=c(2:4, NA), lwd=c(2, 2, 2, 2), lty=c(1:4))
    }
	}

	##
	invisible()
}
