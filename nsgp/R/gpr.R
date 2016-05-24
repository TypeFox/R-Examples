# a file for gaussian process regression
# Two different methods (standard gpr and sample gpr where we use duplicate values to estimate standard deviations)
# provides visualizations


#source("gpr_optim.R")
#source('gp_aux.R')
#source('gp_methods.R')


#' @keywords internal
#' @title GPR posterior
#' @description a simple GPR posterior distribution, no parameter learning
#' @param x obs timepoints
#' @param y obs values
#' @param x.targets target timepoints
#' @param noise noise std, a single value or a vector
#' @param kernelfunc a kernel function of (x1,x2) (returns a matrix)
#' @param derivatives compute also derivatives
#' @return a \code{gpsimple}-object with fields
#' \item{x}{timepoints}
#' \item{mean}{GP mean}
#' \item{cov}{covariance matrix}
#' \item{noisestd}{vector of noise std's}
#' \item{mll}{marginal log likelihood}
#' \item{x.obs}{original observation times}
#' \item{y.obs}{original observation values}
gpr_posterior <- function(x, y, x.targets, noise, kernelfunc, derivatives=FALSE) {
	# center data for easier handling
	adjust = mean(y)
	y = y - adjust
	n = length(y)
	n.targets = length(x.targets)
	
	# noise handling, 
	#  we have noise for both targets and observations
	if (is.null(noise) || is.na(noise)) {
		noise = rep(sd(y)/10, n)   # choose a small constant
	} else if (length(noise) == 1) {
		noise = rep(noise, n)
		noise.targets = rep(noise[1],n.targets)
	} else if (length(noise) == n) {   # interpolate
		noise.targets = .interpolate(x, noise, x.targets)
	} else if (length(noise) == n.targets) {
		noise.targets = noise
		noise = noise.targets[match(x, x.targets)]
	}
	
	# compute kernels
	k.tt <- kernelfunc(x, x)
	k.ts <- kernelfunc(x, x.targets)
	k.st <- kernelfunc(x.targets, x)
	k.ss <- kernelfunc(x.targets, x.targets)
	k.y = k.tt + diag(noise)^2
	k.y.inv <- solve(k.y)
	
	# Calculate the model posterior
	f.mean <- as.vector(k.st%*%k.y.inv%*%y + adjust) # cancel the centering
	f.cov <- k.ss - k.st%*%k.y.inv%*%k.ts
	
	# force f.cov symmetric
	f.cov = as.matrix((f.cov+t(f.cov))/2)
	
	# MLL
	mll = -0.5*t(y)%*%k.y.inv%*%y -0.5*log(det(k.y)) - 0.5*n*log(2*pi)
	
	# derivatives
	#	k.st.deriv = .RBF.deriv(x.targets, x)
	#	k.ts.deriv = .RBF.deriv(x, x.targets)
	#	k.ss.deriv2 = .RBF.deriv2(x.targets, x.targets)
	
	#	f.deriv.mean = k.st.deriv%*%k.y.inv%*%y
	#	f.deriv.cov = k.ss.deriv2 - k.st.deriv%*%k.y.inv%*%k.ts.deriv
	
	res = list('x'=x.targets, "mean"=f.mean, "cov"=f.cov, "noisestd"=noise.targets, 'mll'=mll, x.obs=x, y.obs=y)
	class(res) = 'gpsimple'
	
	return (res)
}


#' @export 
#' @title Perform one-sample GP regression
#' @description Computes the optimal GP model by optimizing the marginal likelihood
#' @details Parameter optimization performed through L-BFGS using analytical gradients
#'  with restarts. The input points \code{x} and output values {y} need to be matching
#'   length vectors. If replicates are provided, they are used to estimate dynamic 
#'   observational noise. 
#' 
#'  The resulting GP model is encapsulated in the return object. The estimated posterior
#'   is in \code{targets$pmean} and \code{targets$pstd} for target points \code{x.targets}.
#'    Use \code{\link{plot.gp}} to visualize the GP.
#' @param x input points
#' @param y output values (same length as \code{x})
#' @param x.targets target points
#' @param noise observational noise (variance), either NULL, a constant scalar or a vector
#' @param nsnoise estimate non-stationary noise from replicates, if possible (default)
#' @param nskernel use non-stationary kernel
#' @param expectedmll use an alternative expected mll optimization criteria
#' @param params gaussian kernel parameters: \code{(sigma.f, sigma.n, l, lmin, c)}
#' @param defaultparams initial parameters for optimization (5-length vector)
#' @param lbounds lower bounds for parameters (5-length vector)
#' @param ubounds upper bounds for parameters (5-length vector)
#' @param optim.restarts restarts in the gradient ascent (default=3)
#' @param derivatives compute also GP derivatives
#' @return A \code{gp}-object (list) containing following elements
#'  \item{targets}{data frame of predictions with points as rows and columns..}
#'  \item{_$x}{points}
#'  \item{_$pmean}{posterior mean of the gp}
#'  \item{_$pstd}{posterior standard deviation of the gp}
#'  \item{_$noisestd}{noises (variance)}
#'  \item{_$mll}{the MLL log likelihood ratio }
#'  \item{_$emll}{the EMLL log likelihood ratio}
#'  \item{_$pc}{the posterior concentration log likelihood ratio}
#'  \item{_$npc}{the noisy posterior concentration log likelihood ratio} 
#'  \item{cov}{learned covariance matrix}
#'  \item{mll}{marginal log likelihood value}
#'  \item{emll}{expected marginal log likelihood value}
#'  \item{kernel}{the kernel matrix used}
#'  \item{ekernel}{the EMLL kernel matrix}
#'  \item{params}{the learned parameter vector:}
#'  \item{_$sigma.f}{kernel variance}
#'  \item{_$sigma.n}{kernel noise}
#'  \item{_$l}{maximum lengthscale}
#'  \item{_$lmin}{minimum lengthscale}
#'  \item{_$c}{curvature}
#'  \item{x}{the input points}
#'  \item{y}{the output values}
#' @seealso \code{\link{gpr2sample}} \code{\link{plot.gp}} 
#' @examples
#' # load example data
#' data(toydata)
#' 
#' \dontrun{can take sevaral minutes
#'  # perform gpr
#'  res = gpr1sample(toydata$ctrl$x, toydata$ctrl$y, seq(0,22,0.1))
#'  print(res)}
#' 
#' # pre-computed toydata model
#' data(toygps)
#' print(toygps$ctrlmodel)
gpr1sample <- function(x, y, x.targets, noise=NULL, nsnoise=TRUE, nskernel=TRUE, expectedmll=FALSE, 
							  params=NULL, defaultparams=NULL, lbounds=NULL, ubounds=NULL, 
							  optim.restarts=3, derivatives=FALSE) {
	# center data
	adjust <- mean(y)
	y <- y - adjust
	n <- length(x)
	
	# save original noise
	noise.given = noise
	
	# add measured points to prediction anyways
	x.targets = unique(sort(c(x,x.targets)))
	n.targets <- length(x.targets)
	
	obj = .handlenoise(x,y,x.targets,noise,nsnoise)
	noise.obs = obj$noise.obs
	noise.targets = obj$noise.targets
	
	if (is.null(params)) {
		# do parameter estimation
		# params are: (sigma.f, sigma.n, l, lmin, c)
		# set good defaultparams based on data
		avg_interval = mean(diff(sort(unique(x))))
		max_interval = max(diff(sort(unique(x))))
		min_interval = min(diff(sort(unique(x))))
		diff_y = max(y) - min(y)
		sd_y = sd(y)
		
		if (is.null(defaultparams)) {
			defaultparams <- c(diff_y/2, sd_y/2, max_interval*1.5, min_interval, 0.1)
		}
		if (is.null(lbounds)) { # set lowerbounds for hyperamaters
			lbounds <- c(0.01, sd_y/50, max_interval*0.60, min_interval/2, 0.01)
		}
		if (is.null(ubounds)) { # set upperbounds for hyperparameters
			ubounds <- c(diff_y*50, sd_y*50, max(x)-min(x), (max(x)-min(x))/2, 3)
		}
		if (!all(is.na(noise.targets))) { # predefined noise: fix sigma.n to 1 and use predefined noises (noise.targets)
			defaultparams[2] <- 1
			lbounds[2] <- 1
			ubounds[2] <- 1
		}
		
		a <- 1
		b <- 0
		# no prior for noise is really needed
#		if (all(is.na(noise.targets))) { # learn a prior for optimization if no noise levels given
#			gamma <- try(fitdistr(tapply(y,x,FUN=var), 'gamma', start=list(shape = 1, rate = 0.1),lower=c(0.001,0.001)), silent=TRUE)
#			if (class(gamma) != 'try-error') {
#				a <- gamma[[1]][[1]]
#				b <- gamma[[1]][[2]]
#			}
#		}
		
		# expected case, slow to optimize, only use 'defaultparams' guess
		if (nskernel == TRUE && expectedmll == TRUE) {
			params <- .optimize_expected(x, y, x.targets, noise.targets, a, b, defaultparams, lbounds, ubounds, optim.restarts)
		# non-stationary case, ok to optimize a bit
		} else if (nskernel == TRUE) {
			params <- .optimize_dynamic(x, y, x.targets, noise.targets, a, b, defaultparams, lbounds, ubounds, optim.restarts) 
		# stationary case, not much need to optimize  
		} else { 
			params <- .optimize_static(x, y, x.targets, noise.targets, a, b, defaultparams, lbounds, ubounds, optim.restarts)
		}
	}
	
	# finally assign params as a list for easier handling
	params <- list("sigma.f"=params[1], "sigma.n"=params[2], "l"=params[3], "lmin"=params[4], "c"=params[5])
	if (is.na(params$c))    { params$c <- Inf  }
	if (is.na(params$lmin)) { params$lmin <- 0 }
	
	# if static noise is learned
	#  -> expand static noise into each target timepoint
	if (all(is.na(noise.targets)) | all(is.null(noise.targets))) { 
		noise.obs = rep(params$sigma.n, length(noise.obs))
		noise.targets = rep(params$sigma.n, length(noise.targets))
	}
	
	# compute kernels
	if (nskernel==TRUE) {
		K = .nsRBF.blocks(x, x.targets, params$sigma.f, params$l, params$lmin, params$c)
		k.tt = K$tt
		k.ts = K$ts
		k.st = K$st
		k.ss = K$ss
	} else {
		K = .RBF.blocks(x, x.targets, params$sigma.f, params$l)
		k.tt = K$tt
		k.ts = K$ts
		k.st = K$st
		k.ss = K$ss
	}
	
	k.y  <- k.tt + diag(noise.obs)^2 # noise var
	k.y.inv <- solve(k.y)
	
	# Calculate the model posterior
	f.mean <- k.st%*%k.y.inv%*%y
	f.cov <- k.ss - k.st%*%k.y.inv%*%k.ts
	k.m <- f.cov + k.ss + 2*diag(noise.targets)^2   # expected likelihood kernel
	k.m.inv <- solve(k.m)
	
	# compute GP derivatives
	f.mean.deriv = NA
	f.std.deriv = NA
	f.cov.deriv = NA
	if (derivatives==TRUE) {
		if (nskernel==TRUE) {
			k.st.deriv = .nsRBF.deriv(x.targets, x, params$sigma.f, params$l, params$lmin, params$c)
			k.ts.deriv = .nsRBF.deriv(x, x.targets, params$sigma.f, params$l, params$lmin, params$c)
			k.ss.deriv2 = .nsRBF.deriv2(x.targets, x.targets, params$sigma.f, params$l, params$lmin, params$c)
		} else {
			k.st.deriv = .RBF.deriv(x.targets, x, params$sigma.f, params$l)
			k.ts.deriv = .RBF.deriv(x, x.targets, params$sigma.f, params$l)
			k.ss.deriv2 = .RBF.deriv2(x.targets, x.targets, params$sigma.f, params$l)
		}
		
		f.mean.deriv = k.st.deriv%*%k.y.inv%*%y
		f.cov.deriv = k.ss.deriv2 - k.st.deriv%*%k.y.inv%*%k.ts.deriv
		f.std.deriv = sqrt(diag(f.cov.deriv))
		f.std.deriv[is.na(f.std.deriv)] = 0.001
	}

	# log likelihood of the whole curve
	wholemll <- -0.5* t(y) %*% k.y.inv %*% y -0.5* log(det(k.y)) - 0.5*n*log(2*pi)
	wholeemll <- -0.5*t(f.mean)%*%k.m.inv%*%f.mean -0.5*log(det(k.m)) - 0.5*n.targets*log(2*pi)
	
	# marginal log likelihoods (MLL) for observed timepoints
	mll = rep(0,n.targets)
	for (t in unique(x)) {
		k.y.slice <- as.matrix(k.y[x==t,x==t])
		mll[x.targets==t] <- -0.5*t(y[x==t])%*%solve(k.y.slice)%*%y[x==t] -0.5*log(det(k.y.slice)) -0.5*sum(x==t)*log(2*pi)
	}
	
	# marginalized expected MLL (EMLL) for all target timepoints
#	emll <- -0.5*f.mean^2/diag(k.m) -0.5*log(diag(k.m)) -0.5*log(2*pi)
	
	# EMLL using a window of 2 timepoints for stability
	emll = rep(0,n.targets)
	emll[1] <- -0.5* t(f.mean[1:2]) %*% solve(k.m[1:2,1:2]) %*% f.mean[1:2] -0.5* log(det(k.m[1:2,1:2])) - 0.5*2*log(2*pi)
	for (i in 2:(n.targets-1)) {
		emll[i] <- mean(c(-0.5* t(f.mean[i:(i+1)]) %*% solve(k.m[i:(i+1),i:(i+1)]) %*% f.mean[i:(i+1)] -0.5* log(det(k.m[i:(i+1),i:(i+1)])) - 0.5*2*log(2*pi),
								-0.5* t(f.mean[(i-1):i]) %*% solve(k.m[(i-1):i,(i-1):i]) %*% f.mean[(i-1):i] -0.5* log(det(k.m[(i-1):i,(i-1):i])) - 0.5*2*log(2*pi)))
	}
	emll[n.targets] <- -0.5* t(f.mean[(n.targets-1):n.targets]) %*% solve(k.m[(n.targets-1):n.targets,(n.targets-1):n.targets]) %*% f.mean[(n.targets-1):n.targets] -0.5* log(det(k.m[(n.targets-1):n.targets,(n.targets-1):n.targets])) - 0.5*2*log(2*pi)
	
	# posterior concentrations (PC) for all target timepoints (with or without noise)
	pc  <- -0.5*log(diag(f.cov)) + log(4*pi)
	npc <- -0.5*log(diag(f.cov) + noise.targets^2) + log(4*pi)
	
	# adjust y-values and posterior back to original scale
	f.mean <- f.mean + adjust
	y = y + adjust
	
	targets <- data.frame(x=x.targets,
								 pmean=f.mean,              # posterior mean
								 pstd=sqrt(diag(f.cov)),    # posterior std
								 noisestd=noise.targets,    # model noise, i.e. \sigma_t
#								 plower95=f.mean - 2*sqrt(diag(f.cov)),  # posterior 95% confidence intervals
#								 pupper95=f.mean + 2*sqrt(diag(f.cov)),
#								 nplower95=f.mean - 2*(sqrt(diag(f.cov) + noise.targets^2)),
#								 npupper95=f.mean + 2*(sqrt(diag(f.cov) + noise.targets^2)),
								 pderiv=f.mean.deriv,
								 pderivstd=f.std.deriv,
								 mll=mll,        # marginal log likelihood at timepoint
								 emll=emll,      # expected marginal log likelihood at timepoint
								 pc=pc,          # posterior concentration at timepoint
								 npc=npc)        # noisy posterior concentration at timepoint
	
	res = list('targets'=targets, 'cov'=f.cov, 'derivcov'=f.cov.deriv, 'mll'=wholemll, 'emll'=wholeemll,
				  'kernel'=k.y, 'ekernel'=k.m, 'params'=params, 'x'=x, 'y'=y, 'noise'=noise.given)
	class(res) = 'gp'
	return (res) 
}

#' @export 
#' @title Performs two-sample GP regression
#'  
#' @description Performs gaussian process regression for two time-series: control and case. 
#'  A third null GP model is learned that assumes both data coming from same process. 
#'  Various likelihood ratios between the null and individual models are estimated to 
#'  distinguish when case and control processes are significantly different. Use \code{\link{plot.gppack}} 
#'  to visualize the models.
#'  
#' @details The control and case do not need have same amount of points. The resulting \code{gppack} object
#'  contains the three learned models and the likelihood ratios along \code{x.targets}.
#'  
#' @param x.ctrl input points (control)
#' @param y.ctrl output values (control)
#' @param x.case input points (case)
#' @param y.case output values (case)
#' @param x.targets target points
#' @param noise.ctrl observational noise
#' @param noise.case observational noise
#' @param nsnoise estimate non-stationary noise function from replicates, if available
#' @param nskernel use non-stationary kernel (default)
#' @param expectedmll use expected MLL optimization criteria
#' @param params.ctrl kernel parameters (control)
#' @param params.case kernel parameters (case)
#' @param defaultparams initial parameters for optimization
#' @param lbounds lower bounds for parameter optimization
#' @param ubounds upper bounds for parameter optimization
#' @param lockatzero estimate a pseudo-observation for time 0
#' @param optim.restarts restarts in the gradient ascent (default=3)
#' @param derivatives compute also GP derivatives
#' @return a \code{gppack}-object that contains
#' \item{ctrlmodel}{the gp-object corresponding to the control data}
#' \item{casemodel}{the gp-object corresponding to the case data}
#' \item{nullmodel}{the gp-object corresponding to the shared null data}
#' \item{ratios}{the log likelihood ratios between the control and case against the null model, contains..}
#' \item{______$mll}{marginal log likelihood ratio}
#' \item{______$emll}{expected marginal log likelihood ratio}
#' \item{______$pc}{log posterior concentration ratio}
#' \item{______$npc}{log noisy posterior concentration ratio}
#'
#' @seealso \code{\link{gpr1sample}} \code{\link{plot.gppack}}
#'
#' @examples
#' # read toy data
#' data(toydata)
#' 
#' \dontrun{can take several minutes
#'  # perform two-sample regression
#'  res = gpr2sample(toydata$ctrl$x, toydata$ctrl$y, toydata$case$x, toydata$case$y, seq(0,22,0.1))
#'  print(res)}
#'
#' # pre-computed model for toydata
#' data(toygps)
#' print(toygps)
gpr2sample <- function(x.ctrl, y.ctrl, x.case=NULL, y.case, # input data
							  x.targets,                           # output times
							  noise.ctrl=NULL, noise.case=NULL,    # observational noise (var)
							  nsnoise = TRUE,                      # noise options 
							  nskernel = TRUE, expectedmll = FALSE,
							  params.ctrl = NULL, params.case = NULL,
							  defaultparams = NULL, lbounds = NULL, ubounds = NULL,
							  lockatzero = FALSE, optim.restarts=3, derivatives=FALSE) {
	# common observations for both case/ctrl
	if (is.null(x.case) || is.na(x.case)) {
		x.case = x.ctrl
	}
	
	# remove missing observations
	x.ctrl = x.ctrl[!is.na(y.ctrl) & !is.null(y.ctrl)]
	y.ctrl = y.ctrl[!is.na(y.ctrl) & !is.null(y.ctrl)]
	x.case = x.case[!is.na(y.case) & !is.null(y.case)]
	y.case = y.case[!is.na(y.case) & !is.null(y.case)]
	
	# construct null model
	x.null = c(x.ctrl,x.case)
	y.null = c(y.ctrl,y.case)
	noise.null = NULL
	
	x.targets = unique(sort(c(x.ctrl,x.case,x.targets)))
	n.targets <- length(x.targets)
	n.ctrl <- length(x.ctrl)
	n.case <- length(x.case)
	n.null = length(x.null)
		
	# try to learn noise from replicates using splines
	#     .. if no replicates, proceed with static noise
	if (nsnoise==TRUE && (is.null(noise.ctrl) || is.null(noise.case))) {
		noise.null = .learnnoise(x.null, y.null, x.targets)
	}
	if (nsnoise==TRUE && is.null(noise.ctrl)) {
		noise.ctrl = .learnnoise(x.ctrl, y.ctrl, x.targets)
	}
	if (nsnoise==TRUE && is.null(noise.case)) {
		noise.case = .learnnoise(x.case, y.case, x.targets)
	}
	
	# compute initial GP to find out the where 0h value is
	# .. add a virtual data point to zero hours with estimated variance
	if ((0 %in% x.ctrl)==FALSE && lockatzero==TRUE) {
		res.init <- gpr1sample(x.ctrl, y.ctrl, x.targets, noise.ctrl, 
									  nskernel=nskernel, expectedmll=expectedmll, 
									  params=params.ctrl, defaultparams=defaultparams, 
									  lbounds=lbounds, ubounds=ubounds, optim.restarts=optim.restarts)
		x.ctrl = c(0, x.ctrl)
		x.case = c(0, x.case)
		x.null = c(0, x.null)
		y.ctrl = c(res.init$targets$pmean[1], y.ctrl)
		y.case = c(res.init$targets$pmean[1], y.case)
		y.null = c(res.init$targets$pmean[1], y.null)
	}
	
	# run GPR's 
	res.ctrl <- gpr1sample(x.ctrl, y.ctrl, x.targets, noise.ctrl, nsnoise=nsnoise, nskernel=nskernel, 
								  expectedmll=expectedmll, params=params.ctrl, defaultparams=defaultparams, 
								  lbounds=lbounds, ubounds=ubounds, optim.restarts=optim.restarts, derivatives=derivatives)
	res.case <- gpr1sample(x.case, y.case, x.targets, noise.case, nsnoise=nsnoise, nskernel=nskernel, 
								  expectedmll=expectedmll, params=params.case, defaultparams=defaultparams, 
								  lbounds=lbounds, ubounds=ubounds, optim.restarts=optim.restarts, derivatives=derivatives)
	res.null <- gpr1sample(x.null, y.null, x.targets, noise.null, nsnoise=nsnoise, nskernel=nskernel, 
								  expectedmll=expectedmll, params=NULL,        defaultparams=defaultparams, 
								  lbounds=lbounds, ubounds=ubounds, optim.restarts=optim.restarts, derivatives=derivatives)
	
	# compute likelihood ratios for comparing the intervals of significant difference
	ratios <- list()
	ratios$pc   <- 0.5*res.ctrl$targets$pc   + 0.5*res.case$targets$pc   - res.null$targets$pc
	ratios$npc  <- 0.5*res.ctrl$targets$npc  + 0.5*res.case$targets$npc  - res.null$targets$npc
	ratios$mll  <- 0.5*res.ctrl$targets$mll  + 0.5*res.case$targets$mll  - res.null$targets$mll
	ratios$emll <- 0.5*res.ctrl$targets$emll + 0.5*res.case$targets$emll - res.null$targets$emll
	
	res = list('ctrlmodel'=res.ctrl, 'casemodel'=res.case, 'nullmodel'=res.null, 'ratios'=ratios)
	class(res) = 'gppack'
	
	return (res)
}
