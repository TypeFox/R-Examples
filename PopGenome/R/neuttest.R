###############################################################################
#
# FUNCTIONS: calcPi, thetaH, fuliD, fuliF, faywuH, tajimaD
#
# A variety of summary tests on simulated data are performed in the 
# functions below. Please find detailed comments above each function 
# declaration
# 
#
# FUNCTION CALLS:	none 	
#			
# PARAMETERS:
#
# RETURN VALUES:	
#
# AUTHOR:			Niroshan Nadarajah <niroshan.nadarajah@uni-duesseldorf.de>
#
# LAST MODIFIED:	10/11/03
#
###############################################################################

## the average number of pairwise differences between sequences (labeled pi)
calcPi <- function(nsam, nsegsites, haplotypes, allelFreq) {
	pi = 0.0
	
	nsamc = as.numeric(nsam)
	nnm1 = nsamc/(nsamc-1)
	for (i in 1:nsegsites) {
		p1 = allelFreq[i]/nsamc
		pi = pi+2*p1*(1-p1)*nnm1
	}
	
	return(pi)
}

## Fay & Wu's ThetaH
thetaH <- function(nsam, nsegsites, haplotypes, allelFreq) {
	pi = 0.0
	
	nsamc = as.numeric(nsam)
	nnm1 = nsamc/(nsamc-1)
	for (i in 1:nsegsites) {
		p1 = allelFreq[i]/nsamc
		pi = p1*p1
	}
	
	return (pi*2/(nsamc*(nsamc-1)))
}


## calculate Fu and Li's D 1993 (with outgroup) 
fuliD <- function (nsam, allelFreq, nsegsites, coef) {
	d <- -10000
	
	if (nsegsites == 0 || coef[1] < 1.5) {
		return(d)	
	}
	
	re <- allelFreq[1]
	an <- coef[1]
	vd <- coef[9]
	ud <- coef[10]
	d  <- ( (nsegsites - an*re) / sqrt((ud*nsegsites) + (vd*nsegsites*nsegsites)) )
	return(d)
}


#double fl_d2(int sample_size,int fr1w,int S, double *coef) /* NO outgroup */
#{
#	double an;
#	int n;
#	double ud2,vd2;
#	int rs;
#	double D2 = -10000;
#        
#        if(S == 0 || *(coef+0) < 1.51) return(-10000);
#	
#	rs = fr1w;
#
#	n = sample_size;
#	an = *(coef+0);
#	
#	vd2 = *(coef+4);
#	ud2 = *(coef+5);
#	D2  = ((double)S/an - (double)rs*(((double)n-(double)1)/(double)n)) /
#	      sqrt(ud2*(double)S + vd2*(double)S*(double)S);
#
#	return D2;
#}


## Fu and Li's F 1993 (with outgroup)
fuliF <- function(nsam, allelFreq, nsegsites, pi, coef) {
	if ( (nsegsites == 0) || (coef[1] < 1.5) ) {
		return(-10000);
	} 
	
	
	re <- allelFreq[1];		
	vf <- coef[11];
	uf <- coef[12];
	
	F  <- (pi - re) / sqrt(uf * nsegsites + vf * nsegsites * nsegsites);
	
	return(F);
}


## Fay and Wu H Test 2000
faywuH <- function (nsam, allelFreq, pi) {
	
	if(nsam < 2) 
		return(-10000);
	
	Th <- 0.0;
	for(i in 1:nsam) {
		Th <- Th + (allelFreq[1]+i)*(i*i);
	} 
	
	Th <- 2/(nsam*(nsam-1));
	
	H <- pi - Th;
	
	return(H);
}


#calc_k <- function(nsam, pi) {
#	comb = (nsam * (nsam-1) / 2);	
#	return(pi/comb);
#}


## 
tajimaD <- function(nsam, nsegsites, pi) {
	a1 = a1f(nsam)
	a2 = a2f(nsam)
	b1 = b1f(nsam)
	b2 = b2f(nsam)
	c1 = c1f(a1, b1)
	c2 = c2f(nsam, a1, a2, b2)
	e1 = e1f(a1, c1)
	e2 = e2f(a1, a2, c2)
	
	return( (pi - (nsegsites/a1)) / sqrt((e1*nsegsites) + ((e2*nsegsites) * (nsegsites-1)))  )
}	


## calculate tajima D coefficients
a1f <- function(nsam) {
	a1 = 0.0
	for (i in 1:nsam) {
		a1 = a1+(1/i)
	}
	return(a1)
}

a2f <- function(nsam){
	a2 = 0.0
	for (i in 1:nsam) {
		a2 = a2 + (1/(i*i))
	}
	return(a2)	
}

b1f <- function(nsam) {
	return( (nsam + 1) / (3*(nsam-1)) )
}

b2f <- function(nsam) {
	return( (2*(nsam*nsam+nsam+3)) / (9*nsam*(nsam-1)) )
}

c1f <- function (a1, b1) {
	return( b1 - (1/a1) )	
}

c2f <- function (nsam, a1, a2, b2) {
	return( b2 - ((nsam+2) / (a1*nsam)) + (a2 / (a1*a1))) 
}

e1f <- function (a1, c1) {
	return( c1/a1 ) 
}

e2f <- function (a1, a2, c2) {
	return( c2 / ( (a1*a1)+a2 ) ) 
}