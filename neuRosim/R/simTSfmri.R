simTSfmri <-
function(design=list(), base=0, nscan=NULL, TR=NULL, SNR=NULL, noise=c("none", "white", "temporal", "low-frequency", "physiological", "task-related", "mixture"), type=c("gaussian","rician"), weights, verbose=TRUE, rho=0.2, freq.low=128, freq.heart=1.17, freq.resp=0.2, vee=1){
	
	if(missing(noise)){
		noise <- "white"
	}
		if(missing(type)){
			type <- "gaussian"
		}
	if(noise=="mixture"){
		if(missing(weights)){
			stop("Weights should be provided with noise=mixture.")
		}
		if(length(weights)!=5){
			stop("Weights vector should have 5 elements.")
		}
		if(sum(weights)!=1){
			stop("The sum of the weights vector should be equal to 1.")
		}
	}

	if(length(design)==0){
		act <- base
		sigma <- mean(act)/SNR
		if(is.null(TR)){
			stop("TR value is missing.")
		}
		if(is.null(nscan)){
			stop("nscan value is missing.")
		}
	} else if(length(design)>1){
		stop("Multiple regions are undefined for time series.")
	} else {
		act <- base + rowSums(specifydesign(design[[1]]$onsets, design[[1]]$durations, design[[1]]$totaltime, design[[1]]$TR, design[[1]]$effectsize, design[[1]]$acc, conv=design[[1]]$hrf, param=design[[1]]$param))
		sigma <- mean(act)/SNR
		TR <- design[[1]]$TR
		nscan <- design[[1]]$totaltime/design[[1]]$TR
	}

	if(noise=="none"){
		n <- 0
	}
	if(noise=="white"){
		n <- c(systemnoise(dim=c(1), sigma=sigma, nscan=nscan, type=type, verbose=verbose, vee=vee))
	}
	if(noise=="temporal"){
		n <- c(temporalnoise(dim=c(1), sigma=sigma, nscan=nscan, rho=rho, verbose=verbose))
	}
	if(noise=="low-frequency"){
		n <- c(lowfreqdrift(dim=c(1), freq=freq.low, nscan=nscan, TR=TR, verbose=verbose))
	}
	if(noise=="physiological"){
		n <- c(physnoise(dim=c(1), sigma=sigma, nscan=nscan, TR=TR, freq.heart=freq.heart, freq.resp=freq.resp, verbose=verbose))
	}
	if(noise=="task-related"){
		n <- c(tasknoise(act.image=act, sigma=sigma, type=type, vee=vee))
	}
	if(noise=="mixture"){
		if(weights[1]==0){
			n.white <- 0
		} else {
			n.white <- c(systemnoise(dim=c(1), sigma=sigma, nscan=nscan, type=type, verbose=verbose, vee=vee))
		}
		if(weights[2]==0){
			n.temp <- 0
		} else {
                       	n.temp <- c(temporalnoise(dim=c(1), sigma=sigma, nscan=nscan, rho=rho, verbose=verbose))
		}
		if(weights[3]==0){
			n.low <- 0
		} else {
			n.low <- c(lowfreqdrift(dim=c(1), freq=freq.low, nscan=nscan, TR=TR, verbose=verbose))
		}
		if(weights[4]==0){
			n.phys <- 0
		} else {
			n.phys <- c(physnoise(dim=c(1), sigma=sigma, nscan=nscan, TR=TR, freq.heart=freq.heart, freq.resp=freq.resp, verbose=verbose))
		}
		if(weights[5]==0){
			n.task <- 0
		} else {
			n.task <- c(tasknoise(act.image=act, sigma=sigma, type=type, vee=vee))
		}
		w <- weights
		n <- (w[1]* n.white + w[2]*n.temp + w[3]*n.low + w[4]*n.phys + w[5]*n.task)/sqrt(sum(w^2))
	}

	fmri.ts <- act + n - mean(n)
	return(fmri.ts)
}

