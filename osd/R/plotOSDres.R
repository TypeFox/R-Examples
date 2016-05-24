
plotOSDres <-
function(OSDobj, type=c("eic","tic","s"), ncomp=vector(), comp.spectra=vector(), mztxt.top=5)
{
	plot.type <- match.arg(type, c("eic","tic","s"))
	Ncomp <- ncomp
	if(plot.type=="eic")
	{	
		if(length(Ncomp)==0) Ncomp <- 1:ncol(OSDobj$C)
		matplot(OSDobj$data, type="l", lty=1, lwd=1.5, col="gray", main="Mixture Resolution (EIC)", ylab="Intensity", xlab="Time (scans)", ylim=c(0,max(OSDobj$data*1.1)))
		matplot(OSDobj$C[,Ncomp], type="l", lty=1, add=T)	
	}
	if(plot.type=="tic")
	{
		if(length(Ncomp)==0) Ncomp <- 1:ncol(OSDobj$C)
		comp.tic <- apply(as.matrix(1:ncol(OSDobj$C)),1,function(x) rowSums(OSDobj$C[,x] %*% t(normalize(OSDobj$S[,x]))))
		mat.TIC <- rowSums(OSDobj$data)
		plot(mat.TIC, type="l", lty=1, lwd=2, col="gray", main="Mixture Resolution (TIC)", ylab="Intensity", xlab="Retention Time (scans)", ylim=c(0,max(mat.TIC*1.1)))
		matplot(comp.tic[,Ncomp], type="l", lty=1, add=T)	
	}
	if(plot.type=="s")
	{
		if(length(Ncomp)==0)
		{ 
			warning("No component was selected for plotting. The first component is shown. Select a component to plot by the 'Ncomp' parameter (in this case, Ncomp=1 was used)")	
			Ncomp <- 1
		}	
		if(length(Ncomp)>1) stop("Components are plotted individually. Select a single number by the 'Ncomp' parameter (for instance, Ncomp=1)")
		if(Ncomp>ncol(OSDobj$S)) stop("'Ncomp' out of range")
		
		empiric.spectra <- OSDobj$S[,Ncomp]
		empiric.spectra <- normalize(empiric.spectra)*1000

		mz.len <- mztxt.top
		if(mz.len!=0)
		{
			if(length(which(empiric.spectra!=0))<5) mz.len <- length(which(empiric.spectra!=0))
			main_mz.empiric <- sort(empiric.spectra, decreasing=T, index.return=T)$ix[1:mz.len]
		}
	
		if(length(comp.spectra)==0)
		{	
			plot(empiric.spectra, type="h", main=paste("Component Spectra #", Ncomp), xlab="Mz", ylab="Intensity", ylim=c(0,1100))
			if(mz.len!=0) text(main_mz.empiric, empiric.spectra[main_mz.empiric]+50, labels=main_mz.empiric, cex=0.8)	
		}else{
			plot(empiric.spectra, type="h", ylim=c(-1100,1100), main=paste("Component Spectra #", Ncomp), xlab="Mz", ylab="Intensity")
			lines(normalize(comp.spectra)*(-1000), type="h", col="purple")
			if(mz.len!=0) text(main_mz.empiric, empiric.spectra[main_mz.empiric]+50, labels=main_mz.empiric, cex=0.8)			
		}
	}
}
