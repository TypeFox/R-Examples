BWfilter <-
function(dataset, Cutoff, SampleRate){
	#Cutoff <- Cutoff	# Cutoff frequency (Hz) required
	#SampleRate <- SampleRate <U+00A0> # (e.g. frame rate of video for kinematics)
	Wc <- tan(Cutoff * pi/0.802/SampleRate)
	Wc2 <- Wc^2
	K1 <- sqrt(2)*Wc
	K2 <- Wc2
	a0 <-K2/(1+K1+K2)
	a1 <- 2*a0
	a2 <- a0
	K3 <- a1/K2
	b1 <- -2*a0+K3
	b2 <- 1-2*a0-K3
	fp <- length(dataset) + 1
	
	# First pass
	dat <- dataset
	dat[is.na(dat)] <- 0
	dat[fp] <- 0
	pass1 <- rep(0, fp)
	for(i in 3:fp){
	pass1[i] <- a0*dat[i]+a1*dat[i-1]+a2*dat[i-2]+b1*pass1[i-1]+b2*pass1[i-2]
	
		}
	# Second pass
	pass2 <- rep(0, fp)
	for(i in fp:3){
		pass2[i-2] <- a0*pass1[i-2]+a1*pass1[i-1]+a2*pass1[i]+b1*pass2[i-1]+b2*pass2[i]
	
		}

	return(Fdat <- pass2[1:length(dataset)])


}

