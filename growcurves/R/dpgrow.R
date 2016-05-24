###################################
## generic for non-session models
###################################
#' @include growthCurve.R
#' @include XZcov.R
#' @include relabel.R
#' @include mcmcPlots.R
#' @include summary_quantiles.R
#' @include dpgrowmm.R
NULL

#' Bayesian semiparametric growth curve models.
#'
#' Employs a Dirichlet Process (DP) prior on the set of by-subject random effect parameters
#' under repeated waves of measurements to allow the number of random effect parameters specified 
#' per subject, \code{q}, to be equal to the number of measurement waves, \code{T}.  
#' Random effects are grouped by subject and
#' all \code{q} parameters receive the DP prior.  The resulting joint marginal 
#' distribution over the data is a DP mixture.
#'
#' @param y A univariate continuous response, specified as an \emph{N x 1} matrix or vector, 
#' where \code{N} captures the number of subject-time cases (repeated subject measures).  
#' Data may reflect unequal number of measures per subject.  Missing occasions are left out as no 
#' \code{NA} values are allowed.
#' @param subject The objects on which repeated measures are conducted that serves 
#' as the random effects grouping factor.  Input as an \emph{N x 1} matrix or vector of 
#' subject-measure cases in either integer or character formt; e.g. 
#' \code{(1,1,1,2,2,3,3,3,...,n,n,n)}, where \code{n} is the total
#'  number of subjects.
#' @param trt An integer or character matrix/vector of length \code{N} 
#' (number of cases) indicating treatment
#'	group assignments for each case.  May also be input as length \code{P} vector, 
#'	where \code{P} is the number of unique subjects, indicating subject group assignment.  
#'	Multiple treatment groups are allowed and if the vector is entered as numeric, 
#'	e.g. \code{(0,1,2,3,..)}, the lowest numbered
#'	group is taken as baseline (captured by global fixed effects).  
#'	If entered in character format,
#'	the first treatment entry is taken as baseline.  
#'	If the are no treatment (vs. control) groups,
#'	then this input may be excluded (set to NULL).
#' @param time A univariate vector of length \code{N}, capturing the time 
#' points associated to each by-subject
#'	measure.  Mav leave blank if only one time point (no repeated measures).
#' @param n.random The desired number of subject random effect terms, \code{q}.  
#' Under \code{option = "dp"} may be set equal to the number of measurement 
#' waves, \code{T}.  The \code{y, trt, time} vectors will together
#'	be used to create both fixed and random effect design matrices.  
#'	The random effects matrix will be of the 
#' 	the form, \code{(1, time, ... , time^(n.random - 1))} (grouped, by \code{subject}). 
#'	This formulation is a growth curve model that allows assessment of 
#'	by-treatment effects and by-client growth curves.
#' @param n.fix_degree The desired polynomial order in time to use for 
#' generating time-based fix effects.
#'	The fixed effects matrix will be constructed as, 
#'	\code{(time, ..., time^(n.fix_degree), trt_1,time*trt_1, ... ,
#'	time^(n.fix_degree)*trt_l, trt_L,..., time^(n.fix_degree)*trt_L)}.
#'	If \code{is.null(n.fix_degree) | n.fix_degree == 0 & is.null(trt)} 
#'	time-by-treatment fixed effects and growth curves are not generated.
#' @param formula Nuisance fixed and random effects may be entered in 
#' \code{formula} with the following format,
#'	\code{y ~ x_1 + x_2*x_3 | z_1*z_2 } as an object of class \code{formula}.  The bar, \code{|}
#'	separates fixed and random effects.  If it
#'	is only desired to enter either fixed or random effects, but not both then the \code{|} may be 
#'	omitted.  Note: the nuisance random effects are assumed to be grouped by subject.  
#'	The fixed and random effects values may change with
#'	each repeated measure; however, within subject growth curves will keep 
#'	constant \code{z} and \code{x} values between
#'	measurement waves.   It is possible to bypass the growth curve construction by 
#'	leaving \code{y, trt, time, n.random, n.fix_degree} 
#'	blank and entering only \code{formula}, instead.  The model output plots, will, however
#'	exclude growth curves in that event.  If a formula is input 
#'	(which requires response, \code{y}) then
#'	the separate entry of \code{y} may be omitted.  If the parameter \code{y} is input, 
#'	it will be over-written by that from \code{formula}.
#' @param random.only A Boolean variable indicating whether the input formula contains 
#' random (for fixed) effects in the case that only
#' one set are entered.  If excluded and \code{formula} is entered without 
#' a \code{|}, \code{random.only} defaults to \code{FALSE}.
#' @param data a \code{data.frame} containing the variables with names as 
#' specified in \code{formula}, including the response, \code{y}.
#' @param n.iter Total number of MCMC iterations.
#' @param n.burn Number of MCMC iterations to discard.  \code{dpgrow} will 
#' return \code{(n.iter - n.burn)} posterior samples.
#' @param n.thin Gap between successive sampling iterations to save.
#' @param shape.dp Shape parameter under a \emph{c ~ G(shape.dp, 1)} 
#' prior on the concentration parameter of the DP (prior
#' on the set of random effects parameters, \emph{b_1, ..., b_n ~ DP(c,G_0)}
#' where \code{n} is the total number of subjects.
#' @param rate.dp Rate parameter under a \emph{c ~ G(shape.dp, rate.dp)} prior on 
#' the concentration parameter of the DP.
#' @param plot.out A boolean variable indicating whether user wants to return plots with output results.  Defaults to \code{TRUE}.
#' @param option Modeling option, of which there are two: 1. \code{dp} places a DP prior on 
#' the set of subject random effects;
#'	2. \code{lgm} places the usual independent Gaussian priors on the set of random effects.
#' @return S3 \code{dpgrow} object, for which many methods are available to return and
#'  view results.  Generic functions applied
#'	to an object, \code{res} of class \code{dpgrow}, includes:
#'	\item{summary(res)}{ returns \code{call}, the function call made to \code{dpgrow} 
#'	and \code{summary.results}, which contains a list of objects that 
#'	include \emph{95\%} credible intervals for each set of sampled parameters, 
#'	specified as (\code{2.5\%}, mean, \emph{97.5\%}, including fixed and random effects. 
#'	Also contains model fit statistics, including \code{DIC} 
#'	(and associated \code{Dbar}, \code{Dhat}, \code{pD}, \code{pV}), as well as the log pseudo 
#'	marginal likelihood (LPML), a leave-one-out fit statistic.  
#'	Note that for \code{option = "dp"}, \code{DIC} is constructed as \code{DIC3} 
#'	(see Celeaux et. al. 2006), where the conditional likehihood evaluated at the 
#'	posterior mode is replaced by the marginal predictive density. 
#'	Lastly, the random and fixed effects design matrices, \code{X, Z}, are returned that 
#'	include both the user input nuisance covariates appended to the time and treatment-based 
#'	covariates constructed by \code{dpgrow}.}  
#'	\item{print(summary(res))}{ prints contents of summary to console.}
#'  \item{plot(res)}{ returns results plots, including the set of subject random 
#'  effects values and credible intervals, a sample
#'	of by-subject growth curves, mean growth curves split by each treatment and control, 
#'	as well as selected trace plots for number of clusters and for precision parameters 
#'	for the likehilood and random effects.  Lastly, a trace plot
#'	for the deviance statistic is also included.}
#'  \item{samples(res)}{ contains (\code{n.iter - n.burn}) posterior sampling 
#'  iterations for every model parameter, including fixed and random
#'	effects.}
#'	\item{resid(res)}{ contains the model residuals.}
#' @note The intended focus for this package are data where both number of subjects and number of 
#' repeated measures are limited.  A DP prior
#'	is placed on the by-subject random effects to borrow strength across subjects for 
#'	each estimation of each subject's growth curve.  The
#'	imposition of the DP prior also allows the resulting posterior distributions 
#'	over the subject random effects to be non-Gaussian.
#'	The \code{dpgrow} function is very similar to \code{dpgrowmm}; 
#'	only the latter includes a separate set of random effects not grouped
#'	by subject (e.g. for treatment dosages allocated to subjects) mapped 
#'	back to subject-time cases through a multiple membership design matrix. 
#'	The \code{dpgrowmult} function generalizes \code{dpgrowmm} by allowing 
#'	more than one multiple membership effects term. 
#'	See Savitsky and Paddock (2011) for detailed model constructions.
#' @keywords model
#' @seealso \code{\link{dpgrowmm}}
#' @examples 
#' \dontrun{
#' ## extract simulated dataset
#' library(growcurves)
#' data(datsim)
#' ## attach(datsim)
#' ## run dpgrow mixed effects model; returns object of class "dpgrow"
#' shape.dp	= 4
#' res		= dpgrow(y = datsim$y, subject = datsim$subject, 
#'			trt = datsim$trt, time = datsim$time,
#'			n.random = datsim$n.random, 
#'			n.fix_degree = 2, n.iter = 10000, 
#'			n.burn = 2000, n.thin = 10,
#'			shape.dp = shape.dp, option = "dp")
#' plot.results	= plot(res) ## ggplot2 plot objects, including growth curves
#' summary.results = summary(res) ## parameter credible intervals,  fit statistics
#' samples.posterior = samples(res) ## posterior sampled values
#' }
#' @aliases dpgrow
#' @aliases dpgrow.default
#' @author Terrance Savitsky \email{tds151@@gmail.com} Susan Paddock \email{paddock@@rand.org}
#' @references 
#'	S. M. Paddock and T. D. Savitsky (2012) Bayesian Hierarchical Semiparametric Modeling of 
#'	Longitudinal Post-treatment Outcomes from Open-enrollment Therapy Groups, submitted to: JRSS 
#'	Series A (Statistics in Society).
#' @references
#' 	T. D. Savitsky and S. M. Paddock (2011) Visual Sufficient Statistics for Repeated Measures data 
#' 	with growcurves for R, submitted to: Journal of Statistical Software.
#' @export dpgrow 
dpgrow			<- function(y, subject, trt, time, n.random, n.fix_degree, formula, 
                      random.only, data, n.iter, n.burn, n.thin, 
					            shape.dp, rate.dp, plot.out, option)
					            UseMethod("dpgrow")

################################################
## default dispatch method for mm-session models
################################################
dpgrow.default		<- function(y = NULL, subject, trt = NULL, time = NULL, n.random = NULL, n.fix_degree = NULL, formula = NULL, random.only = FALSE, 
					data = NULL, n.iter, n.burn, n.thin = 1,
					shape.dp = 1, rate.dp = 1, plot.out = TRUE, option = "dp")
{ ## start function dpgrow.default

  ############################
  ## check inputs
  ############################
  ## model choices
  if( option != "dp" & option != "lgm" )
  {
	stop("You must pick 1 of 2 modeling options, c('dp','lgm')")
  }

  # data choices
  if( is.null(subject) ) stop("must input 'subject' vector that links subjects to cases (of length equal to the number of cases)")


  if(is.null(n.fix_degree)) 
  {
	if( !is.null(time) & !is.null(trt) ) ## user wants growth curve 
	{
		n.fix_degree = length(unique(time)) - 1
		warning("Since 'n.fix_degree' not input, assumed it is equal to maximum number of unique values in 'time' to generate fixed effects.")
	} ## else, the user wants time-based random effects, but no time-by-treatment based fixed effects - so no growth curve
  }else{
	if( n.fix_degree == 0 ) n.fix_degree <- NULL
  }

  if( !is.null(y) ) 
  {
	if( length(subject) != length(y) ) stop("y and subject must be input in subject-time case format")
  }

  if( !is.null(trt) )
  {
	if( length(subject)  != length(trt) )
  	{
        	if( length(trt) == length( unique(subject)) ) ## input in subject, rather than case format
		{
			dat.trt		= data.frame(cbind(unique(subject), trt))
			names(dat.trt)	= c("subject","trt")
			subj.mat	= as.data.frame(subject)
			names(subj.mat)	= "subject"
			dat.trt		= merge(subj.mat,dat.trt,by="subject",all.x=T)
			trt		= dat.trt$trt ## now in case format
			
		}else{ 
			stop("the 'subject' and 'trt' vectors should have length = number of (subject-repeated measures) cases")
		}
	}
  }else{ ## is.null(trt)
	trt = matrix(0, length(subject), 1)
  }
 
  ## data choices - test for formula content in the case random.only == NULL
  if( !is.null(formula) )
  {
  	cov = as.character(formula)[[3]]
  	two.part = grep('\\|',cov)
  	not2part.test	= !length(two.part)  ## true if NOT 2part
  	if( not2part.test == TRUE )
  	{
		if( is.null(random.only) )
		{
			stop("The formula is only 1 part - either fixed or random effects - but not both, so must input a boolean value for 'random.only'")
		}else{ 
			if( random.only == FALSE ) ## user inputs no nuisance random effects
			{
				if( is.null(n.random) ) ## user also inputs no time-based random effects
				{
					stop("Data must include random effects; either input in 'formula' and 'data' or generated by 'time' and 'n.random'.")
				}
			}
		}
  	}
  }else{ ## is.null(formula) == TRUE
	if( is.null(n.random) ) ## user also doesn't input any time-based random effects
	{
		stop("Data must include random effects; either input in 'formula' and 'data' or generated by 'time' and 'n.random'.")
	}
  }

  if( is.null(data) & is.null(time) )
  {
	stop("Input data must be supplied to run model; e.g. (subject,time,trt,n.random) for growth curve and/or 'data' for nuisance covariates.")
  }

  if( !is.null(time) )
  {
	if( any(is.na(time)) | any(is.na(trt)) ) stop("No NA's allowed in c(time,trt) vectors")
  }

  if(!is.null(data))
  {
	if( any( is.na(data) ) ) stop("No NA's allowed in 'data' matrix")
	if( nrow(data) != length(subject) ) stop("Input data.frame must contain number of rows equal to number of subject-measure cases")
  }

  if(any( is.na(subject) ) ) stop("Subject vector not allowed to contain NA's")

  if( is.null(y) & is.null(data) ) stop("Response must be input, either through vector input, 'y', or through 'formula' and 'data'")

  if( is.null(n.random) & !is.null(time) ) stop("Must input 'n.random', number of random effects, to construct growth curve random effects")

  #########################################################################
  ## run mixed effects model engine and produce posterior samples and plots
  #########################################################################

  #############################################################################
  ## re-cast subject identifier inputs to be sequential - subject, subj.aff, trt
  #############################################################################
  ## subject
  start		<- 1
  out		<- relabel(label.input = subject, start)
  subject	<- out$label.new
  o		<- order(subject) ## use later to place X, Z, map.subject, map.trt in contiguous order of subject
  subjecti.u	<- out$labeli.u
  map.subject	<- out$dat.label ## colnames = c("label.input","label.new")

  ## trt
  start		<- 0
  out		<- relabel(label.input = trt, start)
  trt		<- out$label.new
  trti.u	<- out$labeli.u
  map.trt	<- out$dat.label

  #################################################################
  ## some subject, session, case lengths for use in subsetting and loops
  #################################################################
  Ncase			= length(subject)
  Nsubject 		= length(unique(subject))
  Nlevel		= length(unique(trt))
  iter.keep		= floor( (n.iter - n.burn)/n.thin )
  if(is.null(n.random)) n.random = min( length(unique(time)),4 ) ## max number of random effects is q = 4, which produces global cubic fit
  if(!is.null(time)) n.waves = length(unique(time)) ## number of measurement waves - used for growth curve generation with nuisance covariates

  ##################################################################
  ## construct fixed and random effect design matrices
  ##################################################################
  out 	<- XZcov(time = time , trt = trt, trt.lab = trti.u, subject = subject, n.random = n.random, n.fix_degree = n.fix_degree, formula = formula, 
		random.only = random.only, data = data) ## re-ordering to contiguous subject for X and Z is contained in the function XZcov
  X	<- out$X
  X.c   <- out$X.c
  X.n   <- out$X.n
  Z	<- out$Z
  Z.n	<- out$Z.n
  Z.c   <- out$Z.c
  if( !is.null(out$y) ) 
  {
	y 	<- out$y  ## over-writes possible duplicative input of y by user (since must be in formula).
  }else{ ## out$y is null, so user separately entered
	y	<- y[o] ## re-order y by subject to ensure subject is in contiguous order
  }

  ## reorder remaining objects to subject (in contiguous fashion) where entries indexed by case
  subject		<- subject[o]
  map.subject		<- map.subject[o,]
  map.trt		<- map.trt[o,]
  time			<- time[o] ## used for growth curve plotting

  ## capture number of fixed effects
  Nfixed		= ncol(X)
  Nrandom		= ncol(Z)

  ################################################################
  ## conduct posterior sampling and capture results
  ################################################################
  option = tolower(option)
  if(option == "dp") ## DP
  {
	print("Your chosen option = dp")
  	res 		= dpPost(y, X, Z, subject, n.iter, n.burn, n.thin, shape.dp, rate.dp)
  }else{ ## option == "lgm"
	print("Your chosen option = lgm")
	res		= lgmPost(y, X, Z, subject, n.iter, n.burn, n.thin)
  }
 
  ##################################################################
  ## summary (short-hand) results
  ##################################################################
  summary.results			<- summary_quantiles(model.output = res, Nfixed = Nfixed, Nrandom = Nrandom, Nsubject = Nsubject)
  summary.results$X			<- X
  summary.results$Z			<- Z
  summary.results$map.subject		<- map.subject
  summary.results$time			<- time ## not used in accessor functions; just reporting back to user to let them know that sorted by subject
  summary.results$map.trt		<- map.trt
  summary.results$model 		<- option
  summary.results$n.fix_degree		<- n.fix_degree

  residuals		= colMeans(res$Residuals)

 if( (!is.null(time) & length(unique(time)) > 1) & !is.null(n.fix_degree) )
 {
   ###################################################################
   ## growth curves
   ###################################################################
   
   ## generate growth curves with associated identifiers for plotting
   T			= 10 ## produces sufficiently smooth plot
   min.T		= min(time)
   max.T		= max(time)
   if(n.thin == 1)
   {
	n.thin.gc	= 10
   }else{
	n.thin.gc	= 1
   }

   if( is.null(X.n) & is.null(Z.n)  ) ## no nuisance covariates; only time-based covariates.
   {
   	gc.plot		= growthCurve(y.case = y, B = res$B, Alpha = res$Alpha, Beta = res$Beta, trt.case = trt, trt.lab = trti.u, subject.case = subject, 
				subject.lab = subjecti.u, T = T, min.T = min.T, max.T = max.T, n.thin = n.thin.gc, time.case = time, n.fix_degree = n.fix_degree)
   }else{ ## other fixed effects besides time-based covariates. Note: Either X.n or Z.n may be NULL (but not both), which is handled in the growthCurve function
	gc.plot		= growthCurve(y.case = y, B = res$B, Alpha = res$Alpha, Beta = res$Beta, X.n = X.n, Z.n = Z.n, 
				trt.case = trt, trt.lab = trti.u, subject.case = subject, subject.lab = subjecti.u, T = T, min.T = min.T, max.T = max.T, n.thin = n.thin, 
				n.waves = n.waves, time.case = time, n.fix_degree = n.fix_degree, Nrandom = n.random)
			## memo: if have nuisance covariates, need input of Nrandom to construct time-based random effects since Nrandom > n.random
   }

 } ## end conditional statement on creating growth curves

 
 if(plot.out == TRUE)
 {
   ##################################################################
   ## plots
   ################################################################## 
   ## memo: if(option == "lgm") then summary_quantiles excludes M, which means is.null(res$M)

   plot.results = mcmcPlots(subjecti.u = subjecti.u, bmat.summary = summary.results$bmat.summary,  
				M = res$M, Taub = res$Taub, Taue = res$Taue, Deviance = res$Deviance)
   	 
 } #end conditional statement on whether to plot 

   
 ##################################################################
 ## function output
 ##################################################################
 if( plot.out == TRUE  )
 {
    if( (!is.null(time) & length(unique(time)) > 1) & !is.null(n.fix_degree) ) ## growth curves are plotted from time-based covariates
    {
    	plot.results$p.gcall = gc.plot$p.gcall; plot.results$p.gcsel = gc.plot$p.gcsel
    	if(option == "dp") ## DP
    	{
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, M = res$M, S = res$optPartition[[3]], devres = res$devres, 
			Num = res$Num, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]], bigSmin = res$bigSmin, Residuals = res$Residuals,
			Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, 
			plot.results = plot.results, residuals = residuals, dat.growthCurve = gc.plot$plot.dat, dat.gcdata = gc.plot$dat.data)
    	}else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, devres = res$devres, devres3 = res$devres3,
		      	Residuals = res$Residuals, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, 
			plot.results = plot.results, residuals = residuals, dat.growthCurve = gc.plot$plot.dat, dat.gcdata = gc.plot$dat.data)
    	}
    }else{ ## is.null(time) == TRUE
   	if(option == "dp") ## DP
    	{
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, M = res$M, S = res$optPartition[[3]], devres = res$devres, 
			Num = res$Num, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]], bigSmin = res$bigSmin, Residuals = res$Residuals, Tau.e = res$Taue, Tau.b = res$Taub, 
			summary.results = summary.results, plot.results = plot.results, residuals = residuals)
    	}else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, devres = res$devres, devres3 = res$devres3,
		      	Residuals = res$Residuals, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, 
			plot.results = plot.results, residuals = residuals)
    	} ## end conditional statement on choice	
   } ## end conditional statement on whether is.null(time)
 }else{ ## plot.out = FALSE
   if(option == "dp") ## DP
   {
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, M = res$M, S = res$optPartition[[3]], devres = res$devres, 
			Num = res$Num, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]], bigSmin = res$bigSmin, Residuals = res$Residuals,
			Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, residuals = residuals)
   }else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, devres = res$devres, devres3 = res$devres3,
		      	Residuals = res$Residuals, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, residuals = residuals)
   }
 } ## end conditional statement on plot.out

 ##
 ## return list output for dpgrow.default()
 ##

 resot$call		<- match.call()
 resot$Nrandom     	<- ncol(resot$summary.results$Z)
 resot$Nsubject		<- length(unique(subject))
 resot$subject		<- unique(subjecti.u) ## will employ for labeling B with user input subject labels
 class(resot)		<- c("dpgrow")
 return(resot)


} #end function dpgrow.default()

#####################################################
## .Call statements to C++ functions
#####################################################
#' Run a Bayesian mixed effects model for by-subject random effects with DP prior
#'
#' An internal function to \code{\link{dpgrow}}
#'
#' @export
#' @aliases dpPost
#' @param y An \emph{N x 1} response (of subject-measure cases)
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix.  Assumed grouped by \code{subjects}
#' @param subjects An \emph{N x 1} set of subject identifiers
#' @param niter The number of MCMC iterations
#' @param nburn The number of MCMC burn-in iterations to discard
#' @param nthin The step increment of MCMC samples to return
#' @param shapealph The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'	The rate parameter is set of \code{1}.
#' @param ratebeta The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. Default value is \code{1}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrow}}
dpPost = function (y, X, Z, subjects, niter, nburn, nthin, shapealph, ratebeta) {
    stopifnot(nrow(X) == nrow(Z))
    stopifnot(length(y) == nrow(X))
    res <- .Call("DPre", y, X, Z, subjects, niter, nburn, nthin, shapealph, ratebeta, package = "growcurves")
} ## end function dpPost


#' Run a Bayesian mixed effects model for by-subject random effects with an independent Gaussian prior
#'
#' An internal function to \code{\link{dpgrow}}
#'
#' @export
#' @aliases lgmPost
#' @param y An \emph{N x 1} response (of subject-measure cases)
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix.  Assumed grouped by \code{subjects}
#' @param subjects An \emph{N x 1} set of subject identifiers
#' @param niter The number of MCMC iterations
#' @param nburn The number of MCMC burn-in iterations to discard
#' @param nthin The step increment of MCMC samples to return
#'	The rate parameter is set of \code{1}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrow}}
lgmPost = function (y, X, Z, subjects, niter, nburn, nthin) {
    stopifnot(nrow(X) == nrow(Z))
    stopifnot(length(y) == nrow(X))
    res <- .Call("lgm", y, X, Z, subjects, niter, nburn, nthin, package = "growcurves")
} ## end function lgmPost


####################################
## accessor methods
####################################
#' S3 functions of dpgrow
#'
#' produces quantile summaries for model parameters
#'
#' @param object A \code{dpgrow} object
#' @param ... Ignored
#' @export 
#' @method summary dpgrow
#' @aliases summary.dpgrow
summary.dpgrow <- function(object,...)
{
  res 		<- list(call = object$call, summary.results = object$summary.results)
  class(res) 	<- "summary.dpgrow"
  return(res)
}

#' Print summary statistics for sampled model parameters
#'
#' provides credible intervals of sampled parameters for 
#' \code{dpgrow} object
#'
#' @param x A \code{dpgrow} object
#' @param ... Ignored
#' @export 
#' @method print summary.dpgrow
#' @noRd
print.summary.dpgrow <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCredible Intervals and Fit Statistics\n")
  print(x$summary.results)
}

#' Produce samples of MCMC output
#'
#' provides posterior sampled values for every model parameter of a
#' \code{dpgrow} object
#'
#' @param object A \code{dpgrow} object
#' @param ... Ignored
#' @export samples dpgrow
#' @aliases samples.dpgrow
#' @method samples dpgrow
#' @aliases samples.dpgrow
samples.dpgrow <- function(object,...)
{
  B				<- as.data.frame(object$B)
  names(B)			<- paste(rep(1:object$Nrandom, each = object$Nsubject), rep(object$subject, object$Nrandom), sep=".") ## 1.1, 1.2, ...., 1.299
  Beta				<- as.data.frame(object$Beta)
  names(Beta)			<- colnames(object$summary.results$X)

  if(object$summary.results$model == "dp")
  {
	res 		<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, B = B, 
					Residuals = object$Residuals, M = object$M, S = object$S, Num.per.cluster = object$Num, bigSmin = object$bigSmin, phat = object$phat, ordscore = object$ordscore,
					Tau.b = object$Tau.b, Tau.e = object$Tau.e)
  }else{ ## lgm
	res 		<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, B = B, 
					Residuals = object$Residuals, Tau.b = object$Tau.b, Tau.e = object$Tau.e)
  }

  if( !is.null(object$dat.growthCurve) ) ## Add growth curve data set if user chooses growth curve option
  {
	res$dat.growthCurve = object$dat.growthCurve
  }
 
  class(res) 	<- "samples.dpgrow"
  return(res)
}

#' Produce model plots
#'
#' Builds model plots, including MCMC trace plots, analysis of subject effects and subject growth curves
#'
#' @param x A \code{dpgrow} object
#' @param plot.out A \code{boolean} object.  If \code{TRUE}, plots are rendered.  In either case, plots are stored
#' @param ... Ignored
#' @export 
#' @return res a list object of class \code{plot.dpgrow} of two items:
#'	\item{plot.results}{	\code{ggplot2} plot objects.  See \code{\link{mcmcPlots}}. }
#'	\item{dat.growcurve}{	A \code{data.frame} containing fields \code{c("fit","time","subject","trt")}
#'		with \code{P*T} rows, where \code{P} is the length of \code{subject} and \code{T = 10} are the number of in-subject
#'		predictions for each subject.  This object may be used to construct additional growth curves using - see \code{\link{growplot}}.} 
#'      \item{dat.gcdata}{	A \code{data.frame} containing fields  \code{c("fit","time","subject","trt")} with \code{N} rows, where \code{N} are the
#'		number of subject-time cases.  This object contains the actual data for all subjects used to co-plot with predicted growth curves.}
#' @method plot dpgrow
#' @aliases plot.dpgrow
plot.dpgrow <- function(x, plot.out = TRUE, ...)
{
  if(plot.out == TRUE)
  {
  	l.pr 	= length(x$plot.results)
  	for(i in 1:l.pr)
  	{
		dev.new()
		print(x$plot.results[[i]])
  	}
  }
  res			<- list(plot.results = x$plot.results, dat.growcurve = x$dat.growthCurve, dat.gcdata = x$dat.gcdata)
  class(res)		<- "plot.dpgrow"
  return(res)
}


