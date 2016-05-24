# maRketSim - Generation functions, to generate portfolio, account, and market histories from high-level parameters

# The heirarchy:
# core functions calculate values from values
# market objects calculates from core functions
# bond object and functions calculate from market object
# portfolio calculates from bond object and market object
# account calculates from bond object and market objects held within a portfolio


# --- Generate objects from underlying parameters --- #

# Create a bond portfolio using randomly-generated individual bonds according to user-specified parameters
# mkt=mkt object specifying interest rates, dur=portfolio duration, dur.sd=duration spread, n=number of bonds in portfolio
# max.mat = maximum tolerable maturity of any individual issue
# Can pass along arguments to bond(), like par=, f=, etc.
# type="random" is pure random draws, "random.constrained" constrains each successive bond to be closer and closer to the average duration
genPortfolio.bond <- function(n,mkt,dur,dur.sd,max.mat=30,type="random",f=.5,name=as.character(runif(1)),...) {
	# Validate inputs
	if(class(mkt)!="market") { stop("market must be a market object.\n") }
	if(type!="random"&type!="random.constrained") { stop("Tightly constrained portfolio generation not yet supported.  For now it fluctuates around the target duration.\n") }
	# Calculate constants
	max.dur <- findMaxDur()
	# Calculate portfolio
	bonds <- vector("list",length=n) # Create a list to hold the bonds
	for(num in seq(n)) {
		# Keep the maturity under max.mat years
		mat.Q <- max.mat + 1
		while(mat.Q>max.mat) {
			# Draw a requested duration > 0.5 and ensure that it's not greater than the max achievable given market conditions
			requested.dur <- -99
			while(requested.dur<0.5) {
				if(type=="random") {
					requested.dur <- rnorm(1,dur,dur.sd)
				} else if (type=="random.constrained") {
					bonds.current <- bonds[sapply(bonds,function(x) length(x)!=0)]
					if(length(bonds.current)==0) {
						requested.dur <- rnorm(1,dur,dur.sd)
					} else { # Prevent the random walk
						achieved.dur <- summary(portfolio(bonds=bonds.current,mkt=mkt))$portfolio.sum$dur
						scale <- (n-num)/n
						requested.dur <- rnorm(1,(2*dur-achieved.dur),dur.sd*scale) # Get the duration on the opposite side of our target
					} 
				}
			}
			# Create a bond of the requested duration
			cat("requested.dur=",requested.dur,"\n")
			bonds[[num]] <- bond(mkt=mkt,mat=NA,dur=requested.dur,f=.5,...)
			cat("one worked\n")
			# Check that it's not over max.mat
			mat.Q <- bonds[[num]]$mat
		}
	}
	# After we're done, see how many bonds we have to hard constrain in order to meet our duration target exactly while still staying within our maturity constraints
	if(type=="random.constrained") {
		newbond.mats <- c(max.mat+100)
		n.last <- 0 # number of bonds in each iteration
		while(max(newbond.mats)>max.mat | min(newbond.mats)<0) { # iterate until we've spread out the weighting over enough bonds to have realistic bonds
			n.last <- n.last + 1
			achieved.dur <- summary(portfolio(bonds=bonds[seq(n-n.last)],mkt=mkt))$portfolio.sum$dur
			newbonds <- list()
			requested.dur <- (dur-achieved.dur*(n-n.last)/n)/(n.last/n)
			for(m in seq(n.last)) {
				newbonds[[m]] <- bond(mkt=mkt,mat=NA,dur=requested.dur,f=.5,...)
			}
			newbond.mats <- sapply(newbonds,function(x) x$mat)
		}
		# Insert them into our bond list
		for(p in seq(n.last)) {
			bonds[[n-p+1]] <- newbonds[[p]]
		}
	}
	# Return portfolio
	portfolio(bonds=bonds,mkt=mkt,name=name)
}

# Generate a history.market by applying a quoted function to obtain the interest rates for each period
	# i.fxn is a quoted function of the current period t
	# TODO: Enable this to produce a yield.curve for each period instead
genHistory.market <- function(i.fxn,start.t,end.t,f,...) {
	mkts <- list()
	ts <- seq(start.t,end.t,f)
	for(t in ts) {
		i <- eval(i.fxn) # will take into account the current t
		mkts[[length(mkts)+1]] <-  market(mkts=market.bond(i=i,MM.frequency=f), t=t,  ...)
	}
	history.market(mkts)
}


