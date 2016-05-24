# Simulation tools to estimate the distribution of basic log-linear estimates

#' Simulate basic log-linear CRC experiments
#' 
#' Replicate and summarize the generation and log-linear analysis of data sets that are consistent with
#' arbitrary log-linear models 
#' 
#' @param n.grid A vector of positive integers, by default \code{c(100,300,900,2700)}.  Each integer is the number of
#' population units that are observed in a corresponding collection of simulations.
#' @param n.reps The number of replicates for each integer in \code{n.grid}, i.e., for each population size of interest.
#' @param u.vec A vector of log-linear parameters, excluding the intercept term.  The length of the vector and the order
#' 	of its terms must correspond to the column names of the design matrix produced by \code{make.design.matrix(k)},
#'	where \code{k} is the number of lists.
#' @param p0 Optional: a number in \code{(0,1)}, the fraction of the population that is to be undetected.  See details.
#' @param models See \code{\link{lllcrc}}
#' @param ic See \code{\link{lllcrc}}
#' @param cell.adj See \code{\link{lllcrc}}
#' @param averaging \code{\link{lllcrc}}
#' @param fixed.sample.size Logical: If \code{TRUE}, the simulations fix the number of units that are detected, defining the true
#' 	population size such that the number of units detected is equal to its expectation.  If \code{FALSE},
#' 	the observed population size is variable, such that the integers in \code{n.grid} 
#' 	indicate only the expectations of the corresponding simulation sizes.
#' @details \code{u.vec}, together with the constraint that the multinomial probabilities sum to 1, 
#'	uniquely determines the unspecified intercept term.  Specifying \code{p0} overdetermines
#'	the intercept term.  We rectify this overspecification by adjusting all main effects by the same
#'	additive adjustment \code{a}, where the unique value of \code{a} is approximated with numerical methods.
#'	
#'	Once the log-linear terms are fully specified, we perform multinomial draws to simulate a CRC experiment.
#'	We include the zero cell in the multinomial draw only if \code{fixed.sample.size = TRUE}.
#'	
#'	On each replicate, the data log-linear model search according to the parameters \code{models}, 
#'	\code{ic}, \code{cell.adj}, and \code{averaging} produces an estimate of the missing cell. The
#'	main matrix \code{res} of simulation results stores the ratios of the estimated missing cell over
#'	the 'true' missing cell.
#' @return A list of class \code{llsim}, for "log-linear simulations".  The list contains the set of multinomial
#' 	capture pattern probabilities \code{p}, the matrix \code{res} of simulation results, and many of the 
#' 	arguments to the \code{llm.sim}.
#' @author Zach Kurtz
#' @examples
#' \dontrun{
#' ## A basic simulation with four lists.  
#' #	Begin by specifying the vector of log-linear parameters.
#' #	The parameters must match the design matrix:
#' names(make.design.matrix(k=4))
#' u.vec = initialize.u.vec(k=4)
#' u.vec[5:10] = 2
#' ## Run the simulation with an adjustment to the main effects in
#' #	u.vec such that the probability of nondetection is 0.5.
#' sim = llm.sim(n.grid = c(100,300,900,2700), n.reps = 10, u.vec, 
#' 	p0 = 0.5, ic = "BIC", cell.adj = FALSE)
#' # View the results
#' plot(sim)
#' }
#' @export llm.sim
llm.sim = function(n.grid = c(100,300,900,2700), n.reps = 100, u.vec,
	p0 = NULL, models = NULL, ic = "BICpi", cell.adj = TRUE, averaging = FALSE, fixed.sample.size = FALSE)
{
	# Figure out the number of lists
	if(length(u.vec) == 6){k = 3
	}else if(length(u.vec) == 14){k = 4
	}else if(length(u.vec) == 30){k = 5
	}else{ stop("The given u.vec is not compatible with k = 3, 4, or 5")
	}

	# Determine the set of models for model search
	if(is.null(models)) models = make.hierarchical.term.sets(k)
	des = data.matrix(make.design.matrix(k))
	if(!identical(names(u.vec), colnames(des))){
		stop(paste("u.vec must be named with the same names and name order given in the\n",
			"biggest model returned by make.hierarchical.term.sets(k)"))}

	# Determine the multinomial probabilities by u.vec
	if(!is.null(p0)){
		# Compute the number main.adj to add to all main effects
		#	such that the implied intercept term is consistent with p0
		u.vec = zero.inflate(u.vec, p0, k, des)
	}
	p = get.p.from.u(u.vec, des, k)

	# Call the log-linear simulation workhorse
	des = data.frame(des)
	des$c = rep(NA, nrow(des))
	res = matrix(NA, nrow = n.reps, ncol = length(n.grid))
	colnames(res) = paste("n=", as.character(n.grid), sep = "")
	s.grid = n.grid
	# If we're not using a fixed observed sample size, we set the true population size
	#	to satisfy E(observed) = n.grid, approximately
	if(!fixed.sample.size) s.grid = round(n.grid/(1-p$p0)) 
	for(i in 1:length(n.grid)) {
		n = n.grid[i]
		res[,i] = replicate(n.reps, one.llm.sim(size = s.grid[i], k, p, des, 
			models, ic, cell.adj, averaging, fixed.sample.size))
	}

	out = list(p = p, res = res, n.grid = n.grid, u.vec = u.vec, ic = ic, cell.adj = cell.adj, 
		averaging = averaging, fixed.sample.size = fixed.sample.size)
	class(out) = "llsim" # log-linear simulation
	return(out)
}

#' Initialize log-linear parameters
#'
#' A tool for setting up the simulations of \code{\link{llm.sim}}.
#' 
#' @param k The number of lists to be modeled 
#' @return A vector of log-linear parameters, all initialized to zero, corresponding to the columns of
#'	the most general design matrix (but no Rasch terms).
#' @author Zach Kurtz
#' @export initialize.u.vec
initialize.u.vec = function(k)
{
	names.u = colnames(data.matrix(make.design.matrix(k)))
	u.vec = rep(0, length(names.u))
	names(u.vec) = names.u
	return(u.vec)
}

one.llm.sim = function(size, k, p, des, models, ic, cell.adj, averaging, fixed.sample.size = FALSE)
{
	# Multinomial sampling:
	if(fixed.sample.size){
		des$c = rmultinom(1, size, p$p.obs)
		c0 = size*p$p0/(1-p$p0)
	}else{
		p.vec = c(as.numeric(p$p.obs), p$p0)
		draws = rmultinom(1, size, p.vec)
		des$c = draws[-length(draws),]
		c0 = draws[length(draws)]
	}

	# Optionally, the cell adjustment
	if(cell.adj) des$c = des$c + 1/2^(k-1)
	
	# Loglinear modelling:
	icd = ic.all(models, ddat = des, ic, normalized = FALSE)
	if(averaging){
		pred = sum(icd[, "est"] * icd[, "wghts"])
	}else{
		winner = which.min(icd[, "score"])
	    	best.terms = models[[winner]]
		pred = icd[winner, "est"]
	}
	
	# Compute the ratio of the estimated missing cell to the actual or expected missing cell
	return(pred/c0)
}

get.p.from.u = function(u.vec, des, k)
{
	p.obs = t(exp(des %*% u.vec))
	colnames(p.obs) = apply(des[,1:k], 1, paste, collapse = "")
	p0 = saturated.local(p.obs)
	sump = p0+sum(p.obs)
	p.obs = p.obs/sump
	p0 = p0/sump
	return(list(p0=p0, p.obs=p.obs))
}

zero.inflate = function(u.vec, p0, k, des)
{
	loss.a = function(a){
		u.vec[1:k] = u.vec[1:k] + a
		return((get.p.from.u(u.vec, des, k)$p0 - p0)^2)
	}
	u.vec[1:k] = u.vec[1:k] + optimize(f = loss.a,  lower = -10, upper = 10)$minimum
	return(u.vec)
}

#' Plot the output of \code{\link{llm.sim}}
#' 
#' @param x An object of class \code{llsim}
#' @param y.top The upper bound of the plotting window
#' @param probs The interval width, in terms of quantiles
#' @param main Plot title
#' @param ...  Additional parameters to be passed into \code{plot}
#' @author Zach Kurtz
#' @method plot llsim
#' @export
plot.llsim = function(x, y.top = 2, probs = c(0.25, 0.75), main = NULL, ...)
{
	if(is.null(main)) main = paste(nrow(x$res), "replications")
	plot(c(0,0), c(0,0), type = "n", bty = "n", ylim = c(0, y.top), xlim = c(0.5,length(x$n.grid)+0.5),
		ylab = "c0 estimated divided by \"truth\"", xaxt = "n", xlab = "Number of observed units", main = main)
	abline(h = 1, lty = 2)
	for(i in 1:length(x$n.grid)){
		qt = quantile(x$res[,i], probs, na.rm = TRUE)
		mn = mean(x$res[,i], na.rm = TRUE)
		segments(x0 = i, x1 =  i, y0 = qt[1], y1 = qt[2])
		points(x = i, y = mn, pch = 16, cex = 0.8)
		text(x = i, y = 0, labels = colnames(x$res)[i])
	}
}
