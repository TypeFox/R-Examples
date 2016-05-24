#' Tool for bootstrapping
#' 
#' Used by \code{vgam.crc.boots} and \code{lllcrc.boots}.  See these functions
#' for details.
#' 
#' @param boot.list See \code{\link{lllcrc.boots}} or \code{\link{vgam.crc.boots}}.
#' @author Zach Kurtz
#' @references
#' Zwane E and Heijden Pvd (2003). "Implementing the parametric bootstrap
#' in capture-recapture models with continuous covariates." \emph{Statistics &
#' Probability Letters}, \bold{65}, pp. 121-125.
#' @export resample.captures
resample.captures = function(boot.list){ #psi = b.list$psi; b.ydens = b.list$b.ydens
  # location:  core/bootfuns.R
	#nc = sum(boot.list$dat[,"mct"])
	k = nchar(colnames(boot.list$dens)[1])
	cpi0 = boot.list$cpi0 
	dens = boot.list$dens
	dat = boot.list$dat
	nd = nrow(dat)
        # Interpret pi0 as the number of missing units.  But if pi0 is not integer, must somehow round.
	# Here we use bernoulli's to randomly round up or down.
	cpi0.int = floor(cpi0) + rbinom(nd, 1, prob = cpi0%%1) + dat[,"mct"]
	# Append the estimated pi0 with dens
	pmat = cbind(dens, "000"=cpi0/dat[,"mct"]) #, "cpi0.int"=cpi0.int, dat)
	pmat = pmat/rowSums(pmat)

	# Draw multinomial outcomes for ith row:
	ny = 2^k
	draw.cp = function(i) {
		m = sample(1:ny, cpi0.int[i], prob = pmat[i,], replace = TRUE)
		mt = table(m); md = rep(0,ny)
		md[as.numeric(names(mt))] = mt
		return(md)}
	cp.sample = t(sapply(1:nd, draw.cp))

	# Replace the original capture patterns with resampled ones
	y.bits = paste("y", colnames(dens), sep = "")
	dat[,y.bits] = cp.sample[,1:(ny-1)]
	dat[,"mct"] = rowSums(dat[,y.bits])
	return(dat) }


################################################################################################
## Bootstrapping the sampling distribution of the estimated rates of missingness


#' Bootstrap for LLLMs
#' 
#' Implement a single instance of resampling process to simulate a replicate of
#' the \code{lllcrc} model.
#' 
#' \code{vgam.crc} calls this function for each bootstrap replicate.  Many
#' replicates together can be used for empirical confidence interval
#' estimation.  The bootstrap is much like that described by Zwane 2003, but
#' more like Method 2 in Norris and Polluck 1996, as we repeat the process of
#' selecting log-linear terms for each bootstrap iteration, and the multinomial
#' sampling probabilities are based on raw local estimates instead of
#' post-log-linear modeling estimates.
#' 
#' @param boot.list A list of control parameters for the bootstrapping process.
#' See the example under \code{lllcrc-package}
#' @return A list containing two vectors.  The first vector, \code{cpi0}, gives
#' the estimated number of missing units at each point; the second, \code{mct},
#' gives the number of units observed at each point.  (Dividing the first
#' vector by the second gives an estimate rate of missingness)
#' @author Zach Kurtz
#' @references
#' Norris JL and Pollock KH (1996). "Including model uncertainty in
#' estimating variances in multiple capture studies." \emph{Environmental and
#' Ecological Statistics}, \bold{3}, pp. 235-244.
#' @export lllcrc.boots
lllcrc.boots = function(boot.list){ #b.list = bdat.list; imputation.method = "morph.mono"
	# location:  core/bootfuns.R
	if(is.null(boot.list$averaging)) stop("You must define boot.list$averaging as TRUE or FALSE")
	# Resample from the data
	print("lllcrc.boots:  resampling data")
	bdat = resample.captures(boot.list)
	## Weighted data (smoothed) ##
	print("lllcrc.boots:  defining weighted data")
	bsdat = smooth.patterns(dat = as.matrix(bdat), kfrac = boot.list$kfrac, bw = boot.list$bw)
	## Local log-linear models
	print("lllcrc.boots:  fitting local models")
	bloc = apply.ic.fit(ydens = bsdat$hpi, models = boot.list$models, ess = bsdat$ess[,1], mct = bsdat$dat[,"mct"], 
		ic = boot.list$ic, averaging = boot.list$averaging, cell.adj = boot.list$cell.adj, loud = FALSE)
	bloc$bdat = bsdat
	return(bloc)}



#' Bootstrapping for a VGAM CRC model
#' 
#' Implement a single instance of resampling process to simulate a replicate of
#' the \code{vgam.crc} model.
#' 
#' \code{vgam.crc} calls this function for each bootstrap replicate.  Many
#' replicates together can be used for empirical confidence interval
#' estimation.  The bootstrap is much like that described by Zwane 2003, but
#' more like Method 2 in Norris and Polluck 1996, as we repeat the process of
#' selecting log-linear terms for each bootstrap iteration, and the multinomial
#' sampling probabilities are based on raw local estimates instead of
#' post-log-linear modeling estimates.
#' 
#' @param boot.list A list of control parameters for the bootstrapping process.
#' See the example under \code{lllcrc-package}
#' @return A list containing two vectors.  The first vector, \code{cpi0}, gives
#' the estimated number of missing units at each point; the second, \code{mct},
#' gives the number of units observed at each point.  (Dividing the first
#' vector by the second gives an estimate rate of missingness)
#' @author Zach Kurtz
#' @references
#' Norris JL and Pollock KH (1996). "Including model uncertainty in
#' estimating variances in multiple capture studies." \emph{Environmental and
#' Ecological Statistics}, \bold{3}, pp. 235-244.
#' @export vgam.crc.boots
vgam.crc.boots = function(boot.list){ #b.list = bdat.list; imputation.method = "morph.mono"
	# location:  core/bootfuns.R
	# Resample from the data
	print("vgam.crc.boots:  resampling data")
	bdat = resample.captures(boot.list)
	nonzeros = which(bdat$mct>0)
	zeros = which(bdat$mct == 0)
	if(length(zeros)>0) bdat[zeros,c("pi0", "cpi0")] = NA
	nzd = bdat[nonzeros,]
	## Fitting vgam:
	print("vgam.crc.boots:  fitting vgam crc model")
	models = boot.list$models
	n.mod = length(models)
	option.scores = rep(NA, n.mod)
	for(i in 1:n.mod){
		llterms = models[[i]] #loglinear terms
		mod = construct.vgam(sdf = boot.list$sdf, constr.cols = c(models[[i]]),
			constraints = boot.list$cons.mat, dat = nzd)
		option.scores[i] = AICc.vgam(mod)
		}
	llterms = models[[which.min(option.scores)]]

	#Finally, build model:
	mod = construct.vgam(sdf = boot.list$sdf, constr.cols = llterms, constraints = boot.list$cons.mat, dat = nzd)
	fits = fitted(mod)
	bdat[nonzeros,paste("y", names(fits), sep = "")] = fits
	colnames(fits) = substring(colnames(fits),2)
	bdat$pi0[nonzeros] = saturated.local(fits)
	bdat$cpi0 = bdat$pi0*bdat$mct
	out = list(cpi0 = bdat$cpi0, mct = bdat$mct)
	return(out)}

#' Bootstrapping LLMs
#' 
#' For a log-linear model, bootstrap estimates of the corresponding estimates
#' in a way that incorporates the uncertainty due to model selection.
#' 
#' Implements the nonparametric bootstrap of Norris and Pollock (1996), `Method
#' 2'.  The results can be used to create empirical confidence intervals.
#' 
#' @param boot.list A list of control arguments to direct the bootstrap
#' procedure.  See example.
#' @return A vector of bootstrap estimates.
#' @author Zach Kurtz
#' @references
#' Norris JL and Pollock KH (1996). "Including model uncertainty in
#' estimating variances in multiple capture studies." \emph{Environmental and
#' Ecological Statistics}, \bold{3}, pp. 235-244.
#' @export llcrc.flat.boots
llcrc.flat.boots = function(boot.list){ #b.list = bdat.list; imputation.method = "morph.mono"
	# location:  core/bootfuns.R
	if(is.null(boot.list$averaging)) stop("You must define boot.list$averaging as TRUE or FALSE")
	# Resample from the data
	pop = boot.list$pop
	est = boot.list$est
	densi = pop.to.counts(pop$y)
	k = nchar(pop$y[1])
	zero.string = paste(rep(0,k),collapse = "")
	densi[,zero.string] = est
	n = round(sum(densi))
	p = densi/n
	out = rep(NA, boot.list$n.reps)
	for(i in 1:boot.list$n.reps){
		print(i)
		bdens = unlist(sample(names(densi), size = n, prob = p, replace = TRUE))
		dfs = data.frame(matrix(NA, ncol = 1, nrow = n))
		names(dfs) = "y"
		dfs$y = bdens
		dfs = dfs[dfs$y != zero.string, , drop = FALSE]
		est = flat.IC(dfs, rasch = boot.list$rasch, ic = boot.list$ic,
			adjust = boot.list$adjust, averaging = boot.list$averaging)$pred
	out[i] = est
	}
	return(out)
}
