hardThresholding <-
function(xdata, delta, verbose=FALSE, varName=NULL, wavFilter="s8"){

	N <- ncol(xdata)
	n <- nrow(xdata)

	if(n==1){

		if(verbose)
			cat("Univariate reduction dimension\n")

		wt <- wavDWT(xdata, wavelet=wavFilter)
		waveCoefs <- .vectorizeWavelets(wt$data)
		sigest <- mad(wt[[1]])
		
		if(missing(delta))
			delta <- sigest*sqrt(2*log(N))

		idx <- which(abs(waveCoefs)>=delta)
		mht.names <- names(waveCoefs)[idx]

		return(list("mht.names"=mht.names, "estimatedDesign"=waveCoefs[mht.names]))

	}else{

		J <- log2(N)
		design <- t(apply(xdata, MARGIN=1, FUN=function(v) .vectorizeWavelets(wavDWT(v, wavelet=wavFilter)$data)))
		levNames <- numeric(0)

		for(j in J:1){
			K <- 2^(J-j)
			levNames <- c(levNames, paste("d", j, "_", (1:K)-1, sep=""))
		}

		levNames <- c(paste("s", J, sep=""), levNames)
		colnames(design) <- levNames

		if(!is.null(varName)) colnames(design) <- paste(varName, levNames, sep="_")

		sigest <- mad(as.numeric(design[,(2^(J-1)+1):(2^J)]))
		normDesign <- apply(design, FUN=function(v) sqrt(v%*%v), MARGIN=2)

		if(missing(delta)){
			x <- 2*log(N)
			delta <- sigest * sqrt(2*x + 2*sqrt(n * x) + n)
			cat("Automatic threshold ", delta, "\n")
		}else{
			cat("Threshold ", delta, "\n")
		}

		idx <- which(normDesign>=delta)
		mht.names <- names(normDesign)[idx]
		
		if(verbose) cat(length(idx), "selected coefficients using multiple hard-thresholding.\tFilter: ", wavFilter, "\n")

		return(list("mht.names"=mht.names, "estimatedDesign"=design[,mht.names]))
	}
}
