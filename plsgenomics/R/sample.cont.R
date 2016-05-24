### sample.cont.R  (2014-10)
###
###    Generates design matrix X with correlated block of covariates and a continuous random reponse depening on X through gaussian linear model
###
### Copyright 2014-10 Ghislain DURIF
###
###
### This file is part of the `plsgenomics' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA


sample.cont = function(n, p, kstar, lstar, beta.min, beta.max, mean.H=0, sigma.H, sigma.F, sigma.E, seed=NULL) {
	
	### input
	# n : sample size
	# p : number of covariates
	# kstar : number of latent variables
	# lstar : number of blocks of covariates used to generate the response Y
	# seed : seed for the random generator
	
	### tests on input
	
	if((!is.null(seed)) && (!is.numeric(seed)) && (round(seed)-seed!=0)) {
		stop("Message from sample.cont: seed must be integer")
	}
	
	if((!is.numeric(n)) || (round(n)-n!=0) || (!is.numeric(p)) || (round(p)-p!=0) || (!is.numeric(kstar)) || (round(kstar)-kstar!=0) || (!is.numeric(lstar)) || (round(lstar)-lstar!=0) ) {
		stop("Message from sample.cont: n, p, kstar, lstar must be integer")
	}
	
	if((!is.numeric(mean.H)) || (!is.numeric(sigma.H)) || (!is.numeric(sigma.F)) || (!is.numeric(sigma.E)) ) {
		stop("Message from sample.cont: mean.H, sigma.H, sigma.F, sigma.E are not of valid type")
	}
	
	if((sigma.H<0) || (sigma.F<0) || (sigma.E<0)) {
		stop("Message from sample.cont: sigma.H, sigma.F, sigma.E are not of valid type")
	}
	
	if(n<1) {
		stop("Message from sample.cont: n<1, must be strict positive integer")
	}
	
	if(p<1) {
		stop("Message from sample.cont: p<1, must be strict positive integer")
	}
	
	if(kstar<1) {
		stop("Message from sample.cont: kstar<1, must be strict positive integer")
	}
	
	if(lstar>kstar) {
		stop("Message from sample.cont: kstar<lstar, try to use more blocks in X than it actually exists")
	}
	
	if(p<kstar) {
		stop("Message from sample.cont: p<kstar, more blocks than actual covariates")
	}
	
	### random generation
	if(!is.null(seed)) {
		set.seed(seed)
	}
	
	# block size
	block.size <- p %/% kstar
	last.block.size <- p %/% kstar + p %% kstar
	
	# latent variables
	for(i in 1:kstar) {
		assign(paste0("H", i), rnorm(n, mean=mean.H, sd=sigma.H))
	}
	
	### generation of X
	# split into k blocks
	if(last.block.size==0) {
		block.partition <- rep(1:kstar, each=block.size)
	} else {
		block.partition <- c( rep(1:(kstar-1), each=block.size), rep(kstar, each=last.block.size) ) # the last block has more covariates
	}
	
	# index of the jth column in X determines which latent variables to use to generate it
	X <- matrix(data=NA, nrow=n, ncol=p)
	
	X <- sapply(1:p, function(j) {
		
		F <- rnorm(n, mean=0, sd=sigma.F) # noise of column j
		
		return(get(paste0("H", block.partition[j]), inherits=TRUE) + F) # on utilise Hj suivant l'intervalle trouve
		
	})
	
	### generation of Y (binary)
	Y <- numeric(n)
	
	# selected blocks
	block.sel <- sort(sample.int(kstar, size=lstar))
	
	# linear coefficients: non null in the lstar selected blocks
	sel <- (1:p)[block.partition %in% block.sel]
	
	nosel <- (1:p)[-sel]
	
	p0 <- length(sel)
	
	B <- numeric(p)
	B[sel] <- signif(runif(n=p0, min=beta.min, max=beta.max), digits=2)
	B[nosel] <- rep(0, length.out=p-p0)
	
	E = rnorm(n, mean = 0, sd=sigma.E)
	
	Y = X %*% B + E
	
	Y <- as.matrix(Y)
	
	### output:
	# sel: index of variables used to generate Y
	# nosel: index of unused variables
	
	return(list(X=X, Y=Y, residuals=E, sel=sel, nosel=nosel, B=B, block.partition=block.partition, n=n, p=p, kstar=kstar, lstar=lstar, p0=p0, block.sel=block.sel, beta.min=beta.min, beta.max=beta.max, mean.H=mean.H, sigma.H=sigma.H, sigma.F=sigma.F, seed=seed))
	
}
