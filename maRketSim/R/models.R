# Term structure model

# TODO: Make vasicek.discrete.bond which takes a bond object as input and converts the time variable to the 0,h form that vasicek.discrete() uses 

# Vasicek - Discrete time
# mu is what the short term rate converges to
# sigma is the SD of the forecast (determines step size)
# phi determines mean reversion
# h is the time per period
# N is the number of periods total to run the simulation
vasicek.discrete <- function(mu,sigma,h,phi,r.initial,N) {
	# Input validation
	if(phi<0|phi>=1) { stop("phi must be in [0,1)\n") }
	# Build tree of state probabilities
	step <- sigma*sqrt(-2*log(phi))*sqrt(h)
	tree <- matrix(0,ncol=N,nrow=2*N-1) # Tree will hold our rates
	colnames(tree) <- paste("t",seq(N),sep="")
	r.steps <- seq(r.initial-step*(N-1),r.initial+step*(N-1),step) # Centered around r.initial
	price.steps <- exp(-r.steps) # Price of equivalent zero in future dollars
	rownames(tree) <- paste(round(r.steps,4)*100,"%=",round(price.steps,2),"$",sep="") 
	# Iterate
	tree[N,1] <- 1 # Probability of being at r.initial in time 1 is 100% 
	for(t in seq(N-1)) {
		qup <- .5+(mu-r.steps)*sqrt(h)*sqrt(-log(phi))/(sigma*sqrt(8))
		qup[qup<0] <- 0 # By rules of Vasicek model
		qup[qup>1] <- 1 # By rules of Vasicek model
		qdn <- 1-qup
		cat(qup,"\n")
		for(rw in seq(2,nrow(tree)-1)) {
			cat("Time",t,"row",rw,"\n")
			tree[rw-1,t+1] <- tree[rw,t]*qdn[rw]+tree[rw-1,t+1]
			cat(tree[rw,t]*qdn[rw]+tree[rw-1,t+1],"\n")
			tree[rw+1,t+1] <- tree[rw,t]*qup[rw]+tree[rw+1,t+1]
			cat(tree[rw,t]*qup[rw]+tree[rw+1,t+1],"\n")
		}
	}
	# Check
	if(any(apply(tree,2,sum)!=1)) { stop("Something went wrong; columns of our probability tree don't sum to 1\n") } 
	# Return
	structure(list(prob.tree=tree,rates=r.steps,prices=price.steps),class="vasicek.discrete")
}

# Lambda is a measure of excess return per unit of standard deviation (risk-reward relationship), and is constant across all terms
summary.vasicek.discrete <- function(object,lambda,...) {
	# Input validation
	if(class(object)!="vasicek.discrete") { stop("object must be a vasicek.discrete model") }
	# Build term structure data.frame
	T <- ncol(object$prob.tree)
	N <- nrow(object$prob.tree)
	ts <- data.frame(T=seq(T),YTM=rep(NA,T))
	ts$YTM[1] <- object$rates[(N+1)/2] #Initial rate
	# Back out from futures
	
	
	
	structure(list(term.structure=ts),class="sum.vasicek.discrete")
}

# Plot probability histogram
# t is the time period to plot the distribution in; defaults to the last period in the model
plot.vasicek.discrete <- function(x,t=ncol(x$prob.tree),...) {
	plot(x$prob.tree[,t]~x$rates,...)
}
# Plot term structure
plot.sum.vasicek.discrete <- function(x,...) {
	stop("Not yet completed") 
}
