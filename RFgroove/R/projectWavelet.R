projectWavelet <-
function(xdata, wavFilter="s8"){

	N <- ncol(xdata)
	n <- nrow(xdata)
	J <- log2(N)
	design <- t(apply(xdata, MARGIN=1, FUN=function(v) .vectorizeWavelets(wavDWT(v, wavelet=wavFilter)$data)))
	levNames <- numeric(0)

	for(j in J:1){
		K <- 2^(J-j)
		levNames <- c(levNames, paste("d", j, "_", (1:K)-1, sep=""))
	}

	levNames <- c(paste("s", J, sep=""), levNames)
	colnames(design) <- levNames
	design
}
