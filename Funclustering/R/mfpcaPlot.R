
#'
#' This function plots the functional pca. 
#' @title Plot multivariate functional pca
#' 
#' 
#' @param pca is the result of mfpca. In the univariate case mfpcaPlot use the package fda 
#' and will be similar to it's function "plot.pca.fd". 
#' In multivariate functional pca, we will make a graphic window for each dimension.
#' 
#' @param grid specify how to divide the graphics window. grid=c(n,m) divided the widow in to n lines and
#' m columns. If user don't specify grid then he must enter <Enter> to pass to the next graphic.
#' 
#' 			
#' @export 
#' @examples 
#' 
#' # Multivariate
#' # ---------  CanadianWeather (data from the package fda) --------
#' CWtime<- 1:365
#' CWrange<-c(1,365)
#' CWbasis <- create.fourier.basis(CWrange, nbasis=65)
#' harmaccelLfd <- vec2Lfd(c(0,(2*pi/365)^2,0), rangeval=CWrange)
#' 
#' # -- Build the curves ---
#' temp=CanadianWeather$dailyAv[,,"Temperature.C"]
#' CWfd1 <- smooth.basisPar(
#' CWtime, CanadianWeather$dailyAv[,,"Temperature.C"],CWbasis, 
#' Lfdobj=harmaccelLfd, lambda=1e-2)$fd
#' precip=CanadianWeather$dailyAv[,,"Precipitation.mm"]
#' CWfd2 <- smooth.basisPar(
#' CWtime, CanadianWeather$dailyAv[,,"Precipitation.mm"],CWbasis, 
#' Lfdobj=harmaccelLfd, lambda=1e-2)$fd
#' 
#' CWfd=list(CWfd1,CWfd2) 
#' 
#' pca=mfpca(CWfd,nharm=4)
#' mfpcaPlot(pca,grid=c(2,2))
#' 
#' @useDynLib Funclustering

mfpcaPlot <- function(pca,grid=c()) {
	# if we have an univariate functional data, we plot the pca using pca.fd
	if (class(pca)=="pca.fd"){
		#if the user don't specify grid, he must enter <Enter> to pass to the next graphic
		if(length(grid)==0){
			plot.pca.fd(pca)
		}
		# else we divide the graphic window according to the grid parameter
		else{
			dev.new()
			par(mfrow=grid)
			plot.pca.fd(pca)
		}
		
	}
	else{# in the multivariate case
		# dim of the functional data
		dimFd=length(pca)
		for(i in 1:dimFd){
			cat("\n","Dimension ",i," of the functional data","\n")
			#new window for each dimension
			dev.new()#attention to the very large dimensions
			
			#if the user don't specify grid, he must enter <Enter> to pass to the next graphic
			if(length(grid)==0){
				plot.pca.fd(pca[[i]])	
			}
			# else we divide the graphic window according to the grid parameter
			else{
				par(mfrow=grid)
				plot.pca.fd(pca[[i]])	
			}
		}
	}
}
