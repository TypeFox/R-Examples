LombScarglePeriodogram <-
function(times, signal, frequencies=NULL){
	N = length(times);
	if(N<2) return(list(frequencies=c(), powers=c()));
	T = times[N] - times[1];
	
	if(is.null(frequencies)){
		#Get average, minimum and maximum time steps and see if time increments are uniform
		meandt	= T/(N-1);
		dtimes	= diff(times);
		uniform = (max(dtimes)-min(dtimes) < 0.0001*meandt);
		
		f0 = 1.0/(N*meandt);
		maxMode = floor(N/2)-1; #upper limit given by Nyquist frequency
		powers = frequencies = (0:maxMode)*f0;

		if(uniform){
			#time steps are uniform, so use classical DFT formula for optimization reasons
			summationMatrix = outer(0:maxMode, 0:(N-1));
			FRe = (cos(-(2.0*pi*summationMatrix)/N)%*%signal)/N;
			FIm = (sin(-(2.0*pi*summationMatrix)/N)%*%signal)/N;		
			powers = T * (FRe*FRe + FIm*FIm);
			
		}else{
		
			#time steps not uniform, so use more complicated Lomb-Scargle periodogram formula
			Xmean = mean(signal);
			powers[1] = Xmean*Xmean*T;
			modes = 2:(maxMode+1);
			
			#time shifts
			summationMatrix = outer(frequencies[modes],times);
			Cs 		= rowSums(cos(4*pi*summationMatrix));
			Ss 		= rowSums(sin(4*pi*summationMatrix));
			taus 	= atan(Ss/Cs)/(4*pi*frequencies[modes]);
			
			# powers
			summationMatrix = apply(outer(frequencies[modes],times),2,'-',taus); 
			Cs	= rowSums((cos(2*pi*summationMatrix))^2);
			Ss	= rowSums((sin(2*pi*summationMatrix))^2);
			FCs = (cos(-2*pi*summationMatrix)%*%(signal-Xmean))/sqrt(2*N*Cs);
			FSs = (sin(-2*pi*summationMatrix)%*%(signal-Xmean))/sqrt(2*N*Ss);
			powers[modes] = T*(FCs*FCs + FSs*FSs);
			

		}
		
	}else{
		# frequencies are given by caller. Calculate LSP only for those frequencies
		powers 				= frequencies; # allocate space for powers
		Xmean 				= mean(signal);
		zerofreqs 			= which(frequencies==0); # treat zero frequency separately
		powers[zerofreqs] 	= Xmean*Xmean*T;
		nonzerofreqs 		= which(frequencies!=0);
		
		#time shifts
		summationMatrix = outer(frequencies[nonzerofreqs],times);
		Cs 		= rowSums(cos(4*pi*summationMatrix));
		Ss 		= rowSums(sin(4*pi*summationMatrix));
		taus 	= atan(Ss/Cs)/(4*pi*frequencies[nonzerofreqs]);
		
		# powers
		summationMatrix = apply(outer(frequencies[nonzerofreqs],times),2,'-',taus); 
		Cs	= rowSums((cos(2*pi*summationMatrix))^2);
		Ss	= rowSums((sin(2*pi*summationMatrix))^2);
		FCs = (cos(-2*pi*summationMatrix)%*%(signal-Xmean))/sqrt(2*N*Cs);
		FSs = (sin(-2*pi*summationMatrix)%*%(signal-Xmean))/sqrt(2*N*Ss);
		powers[nonzerofreqs] = T*(FCs*FCs + FSs*FSs);
	
	}
	
	return(list(frequencies=frequencies, powers=powers));
}
