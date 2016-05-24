### Defines the class
bces0 <- function(data,dist.c=c("gamma","logn","norm"),dist.e=c("beta","gamma","bern","norm"),
w=1e-6,W=1e-6,n.iter=10000,n.burnin=5000,n.chains=2,robust=TRUE,model.file="model.txt") UseMethod("bces0")

### Call to the actual function
bces0.default <- function(data,dist.c=c("gamma","logn","norm"),dist.e=c("beta","gamma","bern","norm"),
	w=1e-6,W=1e-6,n.iter=10000,n.burnin=5000,n.chains=2,robust=TRUE,model.file="model.txt") {

# BCES0 - Bayesian cost-effectiveness analysis in the presence of structural zero costs
# Based on Baio (2013). Bayesian models for cost-effectiveness analysis in the presence of structural zeros.
#
# data is a named list including values for the following variables:
#  e0,e1,	measure of effectiveness for intervention t=0,1, respectively
#  c0,c1,	costs for intervention t=0,1, respectively
#  X0,X1,	design matrices to estimate the selection model for t=0,1, respectively
#               BCES0 checks that these are centered, and if not does it
#               NB: if X0 and X1 are null (no covariates), assumes a random effect model to
#                   estimate the individual probability of having zero costs
#  H.psi,H.zeta	a vector of fixed hyperparameters for the prior in the positive cost groups.
#               If only one value is passed as argument, then BCEs0 assumes that this is to
#               be used for both treatments being considered
#
# The package implements the following models
# for costs:    Gamma, Log-Normal - these are quite general and probably good enough 
#               to describe most real-life situations.
#               Also considers a Normal model (generally not good, but possibly useful 
#               for transformed data
#
# for benefits: Beta (good to describe QALYs in a 1-year horizon)
#               Gamma (good to describe QALYs in a longer horizon)
#               Bernoulli (good to describe hard measures, such as death)
#               Normal (probably not too generalisable, but could be used for continuous scales)
#
# for the pattern model: if covariates are observed, then runs a simple logistic regression
#               to estimate the probability of observing zero costs. If no covariates are passed
#               as data, then only uses the intercept (marginal probability, assuming no confounders)
#
# w,W           parameters to characterise the degenerate distribution for the structural zeros.
#               w is the mean and W is the standard deviation *on the natural scale* for this 
#               distribution. Default values are w=W=1e-6 (0.000001), which imply effectively a
#               mean of 0 and no variability
#
# robust        if TRUE (default) then implements a Cauchy prior for the pattern model coefficients
#               if FALSE, then implements a vague Normal prior
#
# model.file    string text with the name of the file to which the JAGS model code is saved. The
#               default is "model.txt", but the user can choose a different string
#
# version 23 October 2013

## 0. Sets up the required path and libraries
if(is.element("R2jags",installed.packages()[,1])==F){install.packages("R2jags")}
if(isTRUE(requireNamespace("R2jags",quietly=TRUE))) {

working.dir <- getwd()
# checks that some covariates are available for the selection model
chk <- is.null(data$X0) & is.null(data$X1)
if (chk==TRUE) {		# if not:
	dist.d <- "int"		#  then uses only the intercept (Cauchy prior)
} else {			# alternatively:
	if (robust==TRUE) {
		dist.d <- "cov.cauchy"	# regression model with the observed covariates (Cauchy prior)
	}
	if (robust==FALSE) {
		dist.d <- "cov.norm"	# regression model with the observed covariates (Normal prior)
	}
}

## Now runs JAGS with the required model
# 1. writes the model code to a file using the function "writeModel". This selects the required modules 
#    and then assign the name of the file to the variable "filein"
writeModel(dist.c,dist.e,dist.d,model.file)
filein <- model.file

# 2. Creates variables & check that all works out OK
d0 <- ifelse(data$c0==0,1,0)		# indicator of null cost for t=0
d1 <- ifelse(data$c1==0,1,0)		# indicator of null cost for t=1
n0 <- length(d0)			# sample size for t=0
n1 <- length(d1)			# sample size for t=1

# Defines the data list; this works if there are no covariates (sets 
# parameters for Cauchy priors on the intercepts of the selection model)
dataJags <- list(n0=n0,n1=n1,d0=d0,d1=d1,m.beta0=0,Q.beta0=2.5^-2,m.beta1=0,
		  Q.beta1=2.5^-2,c0=data$c0,c1=data$c1,e0=data$e0,e1=data$e1,
		  H.psi=data$H.psi,H.zeta=data$H.zeta,w=w,W=W)
# If only one value is passed for H.psi,H.zeta, then uses that value for both treatments
if (length(data$H.psi)==1) {
	dataJags$H.psi <- rep(data$H.psi,2)
}
if (length(data$H.zeta)==1) {
	dataJags$H.zeta <- rep(data$H.psi,2)
}
#
if (chk==TRUE) {dataJags$J <- 1}		# needs to specify the number of covariates
						# if none, then J=1 (only the intercept)

# Do this only if there are some "true" covariates (in addition to the intercept)
if (chk==FALSE) {
# first checks that the design matrices have a column of 1s (for the intercept)
    if (class(data$X0)!="matrix") {data$X0 <- as.matrix(data$X0)}
    if (class(data$X1)!="matrix") {data$X1 <- as.matrix(data$X1)}
	if (any(data$X0[,1]!=1)==TRUE) {
		data$X0 <- cbind(rep(1,dim(data$X0)[1]),data$X0)	# if not, adds it
	} 
	if (any(data$X1[,1]!=1)==TRUE) {
		data$X1 <- cbind(rep(1,dim(data$X1)[1]),data$X1)	# if not, adds it
	}
# sets the number of covariates in the pattern model (including the intercept)
	J <- dim(data$X0)[2]

# then checks that the covariates are centered and if not does it
	threshold <- 1e-12	# if a number is within this threshold to 0 then it is effectively 0
	Z0 <- data$X0; Z1 <- data$X1
	for (j in 2:J) {
		if (mean(data$X0[,j])>threshold) {
			Z0[,j] <- scale(data$X0[,j],scale=FALSE)
		}
		if (mean(data$X1[,j])>threshold) {
			Z1[,j] <- scale(data$X1[,j],scale=FALSE)
		}
	}

# finally modifies the data list to include the covariates as well
	dataJags$J <- J
	dataJags$Z0 <- Z0
	dataJags$Z1 <- Z1
	if (robust==FALSE) {
		dataJags$m.beta0 <- rep(0,J)		# m.beta0=Normal mean vector
		dataJags$m.beta1 <- rep(0,J)		# m.beta1=Normal mean vector
		dataJags$Q.beta0 <- .00001*diag(J)	# Q.beta0=Normal precision matrix
		dataJags$Q.beta1 <- .00001*diag(J)	# Q.beta1=Normal precision matrix
	}
}

# 3. Adds elements to the data list, depending on the marginal model for the cost that has been selected
if (dist.c=="logn") {
	# For the log-Normal model needs to define auxiliary variables that are the same as c0 & c1
	# but do not have observed zeros (this does not impact on the model as it is only used in the
	# sub-group where the cost is forced to be 0 anyway
	c0.star <- data$c0; 	c0.star[d0==1] <- .0000001	
	c1.star <- data$c1;	c1.star[d1==1] <- .0000001
	dataJags$c0.star=c0.star
	dataJags$c1.star=c1.star
}

# 4. Defines the parameters vector
params <- c("p","mu.c","mu.e","eta0[1]","lambda0[1]","eta1[1]","lambda1[1]",
	"beta0","beta1","gamma0","gamma1","psi0","psi1")
# the parameters tau0,tau1 are only relevant if benefits are not modelled as Bernoulli
if (dist.e!="bern") {
	params <- c(params,c("tau0","tau1"))
}

# 5. Defines the initial values for the random nodes in the model
#    NB This is a general structure that works for all the model implemented
list.temp <- list(
	beta0=rnorm(dataJags$J,0,1),beta1=rnorm(dataJags$J,0,1),psi0=c(runif(1),NA),
	psi1=c(runif(1),NA),zeta0=c(runif(1),NA),zeta1=c(runif(1),NA),log.tau0=runif(1),
	xi0=runif(1),log.tau1=runif(1),xi1=runif(1)
)
inits <- function(){
	list.temp
}

# 6. Runs JAGS to produce the posterior distributions
n.iter <- n.iter				# number of iterations
n.burnin <- n.burnin				# number of burn-in
n.thin <- floor((n.iter-n.burnin)/500)		# number of thinning so that 1000 iterations are stored
n.chains <- n.chains				# number of Markov chains
mod <- R2jags::jags(dataJags, inits, params, model.file=model.file,
	n.chains=n.chains, n.iter, n.burnin, n.thin,
	DIC=TRUE, working.directory=working.dir, progress.bar="text")

# 7. Defines the output of the function
out <- list(mod=mod,params=params,dataJags=dataJags,inits=inits)
class(out) <- "bces0"
out
}
}


#####
## Print method for objects in the class bces0
print.bces0 <- function(x,...){
print(x$mod,interval=c(.025,.975),digit=3)	# prints the summary results
}

#####
## Traceplot method for objects in the class bces0
plot.bces0 <- function(x,...){
col <- colors()
mdl <- x$mod$BUGSoutput
n.chains <- x$mod$BUGSoutput$n.chains
nodes.to.plot <- c("mu.c[1]","mu.c[2]","mu.e[1]","mu.e[2]","p[1]","p[2]")
xlab <- "Iteration"
ylab <- ""
par(mfrow=c(3,2))
if (n.chains==1) {
	plot(mdl$sims.array[,1,nodes.to.plot[1]],col="blue",t="l",
		xlab=xlab,ylab=ylab,main=expression(paste("Traceplot for ",mu[c0])))

	plot(mdl$sims.array[,1,nodes.to.plot[2]],col="blue",t="l",
		xlab=xlab,ylab=ylab,main=expression(paste("Traceplot for ",mu[c1])))

	plot(mdl$sims.array[,1,nodes.to.plot[3]],col="blue",t="l",
		xlab=xlab,ylab=ylab,main=expression(paste("Traceplot for ",mu[e0])))

	plot(mdl$sims.array[,1,nodes.to.plot[4]],col="blue",t="l",
		xlab=xlab,ylab=ylab,main=expression(paste("Traceplot for ",mu[e1])))

	plot(mdl$sims.array[,1,nodes.to.plot[5]],col="blue",t="l",
		xlab=xlab,ylab=ylab,main=expression(paste("Traceplot for ",p[0])))

	plot(mdl$sims.array[,1,nodes.to.plot[6]],col="blue",t="l",
		xlab=xlab,ylab=ylab,main=expression(paste("Traceplot for ",p[1])))
}
if (mdl$n.chains>=2) {
	cols <- c("blue","red","green","magenta","orange","brown","azure")
	plot(mdl$sims.array[,1,nodes.to.plot[1]],t="l",
		col=col[which(col==cols[1])],xlab=xlab,ylab=ylab,
		main=expression(paste("Traceplot for ",mu[c0])),
		ylim=range(mdl$sims.array[,1:2,nodes.to.plot[1]]))
		for (i in 2:mdl$n.chains) {
			points(mdl$sims.array[,2,nodes.to.plot[1]],t="l",col=col[which(col==cols[2])])
		}

	plot(mdl$sims.array[,1,nodes.to.plot[2]],t="l",
		col=col[which(col==cols[1])],xlab=xlab,ylab=ylab,
		main=expression(paste("Traceplot for ",mu[c1])),
		ylim=range(mdl$sims.array[,1:2,nodes.to.plot[2]]))
		for (i in 2:mdl$n.chains) {
			points(mdl$sims.array[,2,nodes.to.plot[2]],t="l",col=col[which(col==cols[2])])
		}

	plot(mdl$sims.array[,1,nodes.to.plot[3]],t="l",
		col=col[which(col==cols[1])],xlab=xlab,ylab=ylab,
		main=expression(paste("Traceplot for ",mu[e0])),
		ylim=range(mdl$sims.array[,1:2,nodes.to.plot[3]]))
		for (i in 2:mdl$n.chains) {
			points(mdl$sims.array[,2,nodes.to.plot[3]],t="l",col=col[which(col==cols[2])])
		}

	plot(mdl$sims.array[,1,nodes.to.plot[4]],t="l",
		col=col[which(col==cols[1])],xlab=xlab,ylab=ylab,
		main=expression(paste("Traceplot for ",mu[e1])),
		ylim=range(mdl$sims.array[,1:2,nodes.to.plot[4]]))
		for (i in 2:mdl$n.chains) {
			points(mdl$sims.array[,2,nodes.to.plot[4]],t="l",col=col[which(col==cols[2])])
		}

	plot(mdl$sims.array[,1,nodes.to.plot[5]],t="l",
		col=col[which(col==cols[1])],xlab=xlab,ylab=ylab,
		main=expression(paste("Traceplot for ",p[0])),
		ylim=range(mdl$sims.array[,1:2,nodes.to.plot[5]]))
		for (i in 2:mdl$n.chains) {
			points(mdl$sims.array[,2,nodes.to.plot[5]],t="l",col=col[which(col==cols[2])])
		}

	plot(mdl$sims.array[,1,nodes.to.plot[6]],t="l",
		col=col[which(col==cols[1])],xlab=xlab,ylab=ylab,
		main=expression(paste("Traceplot for ",p[1])),
		ylim=range(mdl$sims.array[,1:2,nodes.to.plot[6]]))
		for (i in 2:mdl$n.chains) {
			points(mdl$sims.array[,2,nodes.to.plot[6]],t="l",col=col[which(col==cols[2])])
		}
}
}


#####
writeModel <- function(dist.c,dist.e,dist.d,model.file) {
## Selects the modules in the model code, according to the distributional assumptions

## 1. Selection model
##  a. Standard regression as a function of the observed covariates X0,X1 (Normal prior
sel.mod.cov.norm <- "
model {
# 1. Selection model for c=0
    # a. Intervention t=0
    for (i in 1:n0) {
        d0[i] ~ dbern(pi0.hat[i])
        pi0.hat[i] <- max(.00001,min(.99999,pi0[i]))
        logit(pi0[i]) <- Z0[i,]%*%beta0
    }

    # b. Intervention t=1
    for (i in 1:n1) {
        d1[i] ~ dbern(pi1.hat[i])
        pi1.hat[i] <- max(.00001,min(.99999,pi1[i]))
        logit(pi1[i]) <- Z1[i,]%*%beta1
    }

    # c. Priors on the regression coefficients
    beta0[1:J] ~ dmnorm(m.beta0[],Q.beta0[,])
    beta1[1:J] ~ dmnorm(m.beta1[],Q.beta1[,])

    # d. Average probability of zero cost
    p[1] <- exp(beta0[1])/(1+exp(beta0[1]))   # intervention t=0
    p[2] <- exp(beta1[1])/(1+exp(beta1[1]))   # intervention t=1


"

##  b. Standard regression as a function of the observed covariates X0,X1 (Cauchy prior)
sel.mod.cov.cauchy <- "
model {
# 1. Selection model for c=0
    # a. Intervention t=0
    for (i in 1:n0) {
        d0[i] ~ dbern(pi0.hat[i])
        pi0.hat[i] <- max(.00001,min(.99999,pi0[i]))
        logit(pi0[i]) <- Z0[i,]%*%beta0
    }

    # b. Intervention t=1
    for (i in 1:n1) {
        d1[i] ~ dbern(pi1.hat[i])
        pi1.hat[i] <- max(.00001,min(.99999,pi1[i]))
        logit(pi1[i]) <- Z1[i,]%*%beta1
    }

    # c. Priors on the regression coefficients
    for (j in 1:J) {
        beta0[j] ~ dt(m.beta0,Q.beta0,1)
        beta1[j] ~ dt(m.beta1,Q.beta1,1)
    }
#    beta0[1:J] ~ dmt(m.beta0[],Q.beta0[,],1)
#    beta1[1:J] ~ dmt(m.beta1[],Q.beta1[,],1)

    # d. Average probability of zero cost
    p[1] <- exp(beta0[1])/(1+exp(beta0[1]))   # intervention t=0
    p[2] <- exp(beta1[1])/(1+exp(beta1[1]))   # intervention t=1


"

##  c. Intercept-only model (when no covariates are available)
sel.mod.int <- "
model {
# 1. Selection model for c=0
    # a. Intervention t=0
    for (i in 1:n0) {
        d0[i] ~ dbern(pi0.hat[i])
        pi0.hat[i] <- max(.00001,min(.99999,pi0[i]))
        logit(pi0[i]) <- beta0  #u0[i]	
#        u0[i] ~ dnorm(beta0,tau.u0)      # individual random effect
    }

    # b. Intervention t=1
    for (i in 1:n1) {
        d1[i] ~ dbern(pi1.hat[i])
        pi1.hat[i] <- max(.00001,min(.99999,pi1[i]))
        logit(pi1[i]) <- beta1  
    }

    # c. Priors on the regression coefficients 
    #   NB: uses a Cauchy(0,2.5) distribution to stabilise the intercept
    beta0[J] ~ dt(m.beta0,Q.beta0,1)
    beta1[J] ~ dt(m.beta1,Q.beta1,1)

    # d. Average probability of zero cost
    p[1] <- exp(beta0[1])/(1+exp(beta0[1]))   # intervention t=0
    p[2] <- exp(beta1[1])/(1+exp(beta1[1]))   # intervention t=1


"


## 2. Marginal model for the costs
##  a. Gamma specification
marg.cost.gamma <- "
# 2. Marginal model for the costs
    # a. Intervention t=0
    for (i in 1:n0) {
        c0[i] ~ dgamma(eta0[(d0[i]+1)],lambda0[(d0[i]+1)])
    }
    # Defines the shape parameters for the two components of the mixture
    for (s in 1:2) {
        eta0[s] <- psi0[s]*lambda0[s]
        lambda0[s] <- psi0[s]/pow(zeta0[s],2)
    }

   # b. Intervention t=1
    for (i in 1:n1) {
        c1[i] ~ dgamma(eta1[(d1[i]+1)],lambda1[(d1[i]+1)])
    }
    # Defines the shape parameters for the two components of the mixture
    for (s in 1:2) {
        eta1[s] <- psi1[s]*lambda1[s]
        lambda1[s] <- psi1[s]/pow(zeta1[s],2)
    }

    # Priors on the natural scale parameters 
    psi0[1] ~ dunif(0,H.psi[1]);        zeta0[1] ~ dunif(0,H.zeta[1])
    psi0[2] <- w;                       zeta0[2] <- W
    psi1[1] ~ dunif(0,H.psi[2]);        zeta1[1] ~ dunif(0,H.zeta[2])
    psi1[2] <- w;                       zeta1[2] <- W

    # Weights the cost components by the probability of null costs
    mu.c[1] <- p[1]*psi0[2] + (1-p[1])*psi0[1]
    mu.c[2] <- p[2]*psi1[2] + (1-p[2])*psi1[1]

 
"

##  b. log-Normal specification
marg.cost.logn <- "
# 2. Marginal model for the costs
    # a. Intervention t=0
    for (i in 1:n0) {
        c0.star[i] ~ dlnorm(eta0[(d0[i]+1)],inv.lambda0[(d0[i]+1)])
    }
    # Defines the mean and variance on the log-cost scale for the two components of the mixture
    for (s in 1:2) {
        eta0[s] <- log(psi0[s])-.5*log(1+pow(zeta0[s]/psi0[s],2))
        lambda0[s] <- pow(log(1+pow(zeta0[s]/psi0[s],2)),.5)
        inv.lambda0[s] <- pow(lambda0[s],-2)
    }

    # b. Intervention t=1
    for (i in 1:n1) {
        c1.star[i] ~ dlnorm(eta1[(d1[i]+1)],inv.lambda1[(d1[i]+1)])
    }
    # Defines the mean and variance on the log-cost scale for the two components of the mixture
    for (s in 1:2) {
        eta1[s] <- log(psi1[s])-.5*log(1+pow(zeta1[s]/psi1[s],2))
        lambda1[s] <- pow(log(1+pow(zeta1[s]/psi1[s],2)),.5)
        inv.lambda1[s] <- pow(lambda1[s],-2)
    }

    # Priors on the natural scale parameters
    psi0[1] ~ dunif(0,H.psi[1]) 
    zeta0[1] ~ dunif(0,H.zeta[1])
    psi0[2] <- w	
    zeta0[2] <- W
    psi1[1] ~ dunif(0,H.psi[2])
    zeta1[1] ~ dunif(0,H.zeta[2])
    psi1[2] <- w	
    zeta1[2] <- W

    # Weights the average by the probability of null costs
    mu.c[1] <- p[1]*psi0[2] + (1-p[1])*psi0[1]
    mu.c[2] <- p[2]*psi1[2] + (1-p[2])*psi1[1]


"
##  c. Normal specification
marg.cost.norm <- "
# 2. Marginal model for the costs
    # a. Intervention t=0
    for (i in 1:n0) {
        c0[i] ~ dnorm(eta0[(d0[i]+1)],inv.lambda0[(d0[i]+1)])
    }
    # Defines the mean and variance on the log-cost scale for the two components of the mixture
    for (s in 1:2) {
        eta0[s] <- psi0[s]
        lambda0[s] <- zeta0[s]
        inv.lambda0[s] <- pow(lambda0[s],-2) 
    } 

    # b. Intervention t=1
    for (i in 1:n1) {
        c1[i] ~ dnorm(eta1[(d1[i]+1)],inv.lambda1[(d1[i]+1)])
    }
    # Defines the mean and variance on the log-cost scale for the two components of the mixture
    for (s in 1:2) {
        eta1[s] <- psi1[s]
        lambda1[s] <- zeta1[s]
        inv.lambda1[s] <- pow(lambda1[s],-2) 
    } 

    # Priors on the natural scale parameters
    psi0[1] ~ dunif(0,H.psi[1]) 
    zeta0[1] ~ dunif(0,H.zeta[1])
    psi0[2] <- w
    zeta0[2] <- W
    psi1[1] ~ dunif(0,H.psi[2])
    zeta1[1] ~ dunif(0,H.zeta[2])
    psi1[2] <- w
    zeta1[2] <- W

    # Weights the average by the probability of null costs
    mu.c[1] <- p[1]*psi0[2] + (1-p[1])*psi0[1]
    mu.c[2] <- p[2]*psi1[2] + (1-p[2])*psi1[1]


"


## 3. Conditional model for the effectiveness
##  a. Beta specification
cond.eff.beta <- "
# 3. Conditional model for the benefits
    # a. Intervention t=0
    for (i in 1:n0) {
        e0[i] ~ dbeta(phi0[i]*tau0,(1-phi0[i])*tau0)
        logit(phi0[i]) <- xi0 + gamma0*(c0[i]-mu.c[1])
    }
    gamma0 ~ dnorm(0,0.0001)
    log.tau0 ~ dnorm(0,0.0001)
    tau0 <- exp(log.tau0)
    xi0 ~ dnorm(0,0.0001)
    mu.e[1] <- exp(xi0)/(1+exp(xi0))
    
   # b. Intervention t=1
    for (i in 1:n1) {
        e1[i] ~ dbeta(phi1[i]*tau1,(1-phi1[i])*tau1)
        logit(phi1[i]) <- xi1 + gamma1*(c1[i]-mu.c[2])
    }
    gamma1 ~ dnorm(0,0.0001)
    log.tau1 ~ dnorm(0,0.0001)
    tau1 <- exp(log.tau1)
    xi1 ~ dnorm(0,0.0001)
    mu.e[2] <- exp(xi1)/(1+exp(xi1))
}
"

##  b. Normal specification
cond.eff.norm <- "
# 3. Conditional model for the benefits
    # a. Intervention t=0
    for (i in 1:n0) {
        e0[i] ~ dnorm(phi0[i],tau0)
        phi0[i] <- xi0 + gamma0*(c0[i]-mu.c[1])
    }
    gamma0 ~ dnorm(0,0.0001)
    log.tau0 ~ dnorm(0,0.0001)
    tau0 <- exp(log.tau0)
    xi0 ~ dnorm(0,0.0001)
    mu.e[1] <- xi0
    
   # b. Intervention t=1
    for (i in 1:n1) {
        e1[i] ~ dnorm(phi1[i],tau1)
        phi1[i] <- xi1 + gamma1*(c1[i]-mu.c[2])
    }
    gamma1 ~ dnorm(0,0.0001)
    log.tau1 ~ dnorm(0,0.0001)
    tau1 <- exp(log.tau1)
    xi1 ~ dnorm(0,0.0001)
    mu.e[2] <- xi1
}
"


##  c. Gamma specification
cond.eff.gamma <- "
# 3. Conditional model for the benefits
    # a. Intervention t=0
    for (i in 1:n0) {
        e0[i] ~ dgamma(tau0,kappa0[i])
        kappa0[i] <- tau0/phi0[i]
        phi0[i] <- xi0 + gamma0*(c0[i]-mu.c[1])
    }
    gamma0 ~ dnorm(0,0.0001)
    log.tau0 ~ dnorm(0,0.0001)
    tau0 <- exp(log.tau0)
    xi0 ~ dnorm(0,0.0001)
    mu.e[1] <- xi0
    
   # b. Intervention t=1
    for (i in 1:n1) {
        e1[i] ~ dgamma(tau1,kappa1[i])
        kappa1[i] <- tau1/phi1[i]
        phi1[i] <- xi1 + gamma1*(c1[i]-mu.c[2])
    }
    gamma1 ~ dnorm(0,0.0001)
    log.tau1 ~ dnorm(0,0.0001)
    tau1 <- exp(log.tau1)
    xi1 ~ dnorm(0,0.0001)
    mu.e[2] <- xi1
}
"


##  d. Bernoulli specification
cond.eff.bern <- "
# 3. Conditional model for the benefits
    # a. Intervention t=0
    for (i in 1:n0) {
        e0[i] ~ dbern(phi0[i])
        logit(phi0[i]) <- xi0 + gamma0*(c0[i]-mu.c[1])
    }
    gamma0 ~ dnorm(0,0.0001)
    log.tau0 ~ dnorm(0,0.0001)    # not really needed in this model, but kept for simplicity
    xi0 ~ dnorm(0,0.0001)
    mu.e[1] <- exp(xi0)/(1+exp(xi0))
    
   # b. Intervention t=1
    for (i in 1:n1) {
        e1[i] ~ dbern(phi1[i])
        logit(phi1[i]) <- xi1 + gamma1*(c1[i]-mu.c[2])
    }
    gamma1 ~ dnorm(0,0.0001)
    log.tau1 ~ dnorm(0,0.0001)    # not really needed in this model, but kept for simplicity 
    xi1 ~ dnorm(0,0.0001)
    mu.e[2] <- exp(xi1)/(1+exp(xi1))
}
"


## Combines the required modules in a single model code
mod.sel <- c("cov.cauchy","cov.norm","int")
mod.cost <- c("gamma","logn","norm")
mod.eff <- c("beta","norm","gamma","bern")
lab.sel <- c("Model with covariates (Cauchy prior)","Model with covariates (Normal prior)",
	"Model with intercept only")
lab.cost <- c("Gamma","Log-Normal","Normal")
lab.eff <- c("Beta","Normal","Gamma","Bernoulli")
model <- NULL	# avoids a NOTE when checking the package
for (ms in 1:length(mod.sel)) {
    for (mc in 1:length(mod.cost)) {
        for (me in 1:length(mod.eff)) {
            if (dist.d==mod.sel[ms] & dist.c==mod.cost[mc] & dist.e==mod.eff[me]) {
	        summary.text <- paste0("# ",lab.cost[mc]," marginal model for the cost, ",lab.eff[me],
		    " conditional model for effectiveness\n")
	       txt <- paste0("summary.text,sel.mod.",mod.sel[ms],",marg.cost.",
	                     mod.cost[mc],",cond.eff.",mod.eff[me])
	       cmd <- paste0("model <- paste(",txt,")")
	       eval(parse(text=cmd))
            }
	}
    }
}

filein <- model.file
writeLines(model,con=filein)
}

