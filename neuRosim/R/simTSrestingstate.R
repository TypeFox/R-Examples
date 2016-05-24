simTSrestingstate <- function(nscan, base=0, TR, SNR=NULL, noise=c("none", "white", "temporal", "low-frequency", "physiological", "mixture"), type=c("gaussian","rician"), weights, verbose=TRUE, rho=0.2, freq.low=128, freq.heart=1.17, freq.resp=0.2, vee=1){

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
		if(length(weights)!=4){
			stop("Weights vector should have 4 elements.")
		}
		if(sum(weights)!=1){
			stop("The sum of the weights vector should be equal to 1.")
		}
	}

fl <- 0.01                         # Lower limit of the resting state passband (Hz).
fu <- 0.1                          # Upper limit
Ampl <- 2                          # Maximum percentage BOLD signal change
N <- nscan
f <- 1/TR                          # Sampling frequency
NOS <- 2^ceiling(log2(abs(N)))     # Effective number of samples is the next power of 2 following N.
Nl <- round(fl*NOS/f)              # Express frequencies as fractions of the number of samples.
Nu <- round(fu*NOS/f) 

## Spectrum construction
Z <- rep(0,length=NOS) 
Phase <- 2*pi*runif(NOS/2-1,min=0, max=1)-pi                            # Randomize the phases. 
z <- rep(0,length=NOS/2-1) 
z[Nl:Nu] <- complex(real=cos(Phase[Nl:Nu]),imaginary=sin(Phase[Nl:Nu])) # Construct the passband and force Hermitian.
Z[2:(NOS/2)] <- z                                                       # Lower spectrum half
Z[(NOS/2+2):NOS] <- Conj(z[length(z):1])                                # Upper half

x <- Re(fft(Z,inverse=TRUE)/NOS)   # Inverse Fourier transform (output should be real).
x <- Ampl*x/max(x)                 # Normalize x.
x <- base + x[1:N]                        # Retain only the first N samples (slightly distorts the power spectrum).

	if(noise!="none"){  
	  sigma <- mean(x)/SNR
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
	if(noise=="mixture"){
		if(weights[1]==0){
			n.white <- 0
		} else {
			n.white <- c(systemnoise(dim=c(1), sigma=sigma, nscan=nscan, type=type, verbose=verbose))
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
		w <- weights
		n <- (w[1]* n.white + w[2]*n.temp + w[3]*n.low + w[4]*n.phys)/sqrt(sum(w^2))
	}
	ts <- x + n - mean(n)
	return(ts)
}
