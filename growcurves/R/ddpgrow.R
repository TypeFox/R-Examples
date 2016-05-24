###################################
## generic for non-session models
###################################
#' @include growthCurve.R
#' @include XZcov.R
#' @include relabel.R
#' @include ddpMCMCplots.R
#' @include ddp_quantiles.R
#' @include dpgrowmm.R
NULL

#' Bayesian semiparametric growth curve models.
#'
#' Employs an anova Dependent Dirichlet Process (DDP) prior on the set of by-subject random effect parameters
#' with dependence indexed by multiple membership effects under repeated waves of measurements to allow the 
#' number of random effect parameters specified per subject, \code{q},
#' to be equal to the number of measurement waves, \code{T}.  Random effects are grouped by subject and
#' all \code{q} parameters receive the DP prior.  The resulting joint marginal distribution over the data is a 
#' DP mixture.
#'
#' @param y A univariate continuous response, specified as an \emph{N x 1} matrix or vector, where \code{N}
#'	captures the number of subject-time cases (repeated subject measures).  Data may reflect unequal
#'	number of measures per subject.  Missing occasions are left out as no \code{NA} values are allowed.
#' @param subject The objects on which repeated measures are conducted that serves as the random effects
#'	grouping factor.  Input as an \emph{N x 1} matrix or vector of subject-measure cases in either
#'	integer or character formt; e.g. \code{(1,1,1,2,2,3,3,3,...,n,n,n)}, where \code{n} is the total
#'	number of subjects.
#' @param trt An integer or character matrix/vector of length \code{N} (number of cases) indicating treatment
#'	arm assignments for each case.  May also be input as length \code{n} vector, where \code{n} is
#'	the number of unique subjects, indicating subject arm assignment.  Multiple treatment arms
#'	are allowed and if the vector is entered as numeric, e.g. \code{(0,1,2,3,..)}, the lowest numbered
#'	arm is taken as baseline (captured by global fixed effects).  If entered in character format,
#'	the first treatment entry is taken as baseline.  If the are no treatment (vs. control) arm,
#'	then this vector is composed of a single value for all entries.
#' @param time A univariate vector of length \code{N}, capturing the time points associated to each by-subject
#'	measure.  
#' @param n.random The desired number of subject random effect terms, \code{q}.  May
#'	be set equal to the number of measurement waves, \code{T}.  The \code{y, trt, time} vectors will together
#'	be used to create both fixed and random effect design matrices.  The random effects matrix will be of the 
#' 	the form, \code{(1, time, ... , time^(n.random - 1))} (grouped, by \code{subject}). 
#'	This formulation is a growth curve model that allows assessment of by-treatment effects and by-client growth curves.
#' @param n.fix_degree The desired polynomial order in time to use for generating time-based fix effects.
#'	The fixed effects matrix will be constructed as, 
#'	\code{(time, ..., time^(n.fix_degree), trt_1,time*trt_1, ... ,time^(n.fix_degree)*trt_l, trt_L,..., time^(n.fix_degree)*trt_L)}.
#' @param formula Nuisance fixed and random effects may be entered in \code{formula} with the following format,
#'	\code{y ~ x_1 + x_2*x_3 | z_1*z_2 }.  The bar, \code{|}, separates fixed and random effects.  If it
#'	is only desired to enter either fixed or random effects, but not both then the \code{|} may be omitted.  Note:
#'	the nuisance random effects are assumed to be grouped by subject.  The fixed and random effects valules may change with
#'	each repeated measure; however, within subject growth curves will keep constant \code{z} and \code{x} values between
#'	measurement waves.   It is possible to bypass the growth curve construction by leaving \code{y, trt, time, n.random} 
#'	blank and entering only \code{formula}, instead.  The model output plots, will, however
#'	exclude growth curves in that event.  If a formula is input and a response, \code{y}, is included, then
#'	the parameter input \code{y} may be omitted.  If the \code{y} input is included, it will be over-written by that from \code{formula}.
#' @param random.only A Boolean variable indicating whether the input formula contains random (for fixed) effects in the case that only
#'	one set are entered.  If excluded and \code{formula} is entered without a \code{|}, \code{random.only} defaults to 
#'	\code{FALSE}.
#' @param data A \code{data.frame} containing the variables named in \code{formula}.
#' @param dosemat An \code{n x (sum(numdose)+1)} \code{matrix} object that maps \code{subjects} to treatment dosages.  The first column should be an
#'			intercept column (filled with 1's).  If there is only a single treatment arm, then the number of columns in dosemat should be \code{sum(numdose)}.
#'			There is always a leave-one-out dosage for \code{dosemat}.  For multiple treatment arms, the null (0) treatment is the one left out.
#' @param numdose A vector object containing the number of dosages for each treatment.  So the length should be the same as \code{typetreat}.
#' @param typetreat A vector object specifying the prior formulation for each treatment.   The choices for prior formulations are
#' 			\code{c("car","mvn","ind")}.
#' @param labt An optional vector object (of the same length as \code{typetreat}) providing user names for each treatment.  The names are
#'		used in plot objects.  If \code{NULL}, then the numerical order of treatment entries are used.
#' @param Omega A list object of length equal to the number of treatments chosen with the \code{"car" \%in\% typetreat}.
#'	List element \code{m} contains an \emph{numdose[m] x numdose[m]} numeric matrix to encode the CAR adjacency matrix,
#'	where \code{numdose[m]} is the number of dosages receiving the multivariate CAR prior. 
#'	This input is required only under \code{"car" \%in\% typetreat}. 
#' @param n.iter Total number of MCMC iterations.
#' @param n.burn Number of MCMC iterations to discard.  \code{ddpgrow} will return \code{(n.iter - n.burn)} posterior samples.
#' @param n.thin Gap between successive sampling iterations to save.
#' @param shape.dp Shape parameter under a \emph{c ~ G(shape.dp, rate.dp)} prior on the concentration parameter of the DP (prior
#'	on the set of random effects parameters, \emph{b_1, ..., b_n ~ DP(c,G_0)}
#'	where \code{n} is the total number of subjects.
#' @param rate.dp Rate parameter under a \emph{c ~ G(shape.dp, rate.dp)} prior on the concentration parameter of the DP.
#' @param M.init Scalar value capturing number of initial subject clusters to kick-off MCMC chain.
#' @param plot.out A boolean variable indicating whether user wants to return plots with output results.  Defaults to \code{TRUE}.
#' @return S3 \code{dpgrow} object, for which many methods are available to return and view results.  Generic functions applied
#'	to an object, \code{res} of class \code{dpgrow}, includes:
#'	\item{summary(res)}{ returns \code{call}, the function call made to \code{dpgrow} and \code{summary.results}, which contains
#'			a list of objects that include \emph{95\%} credible intervals for each set of sampled parameters, 
#'			specified as (\code{2.5\%}, mean, \emph{97.5\%}, including fixed and random effects. 
#'			Also contains model fit statistics, including \code{DIC} (and associated \code{Dbar}, \code{Dhat}, \code{pD}, 
#'			\code{pV}), as well as the log pseudo marginal likelihood (LPML), a leave-one-out fit statistic.  
#'			Note that \code{DIC} is constructed as \code{DIC3} (see Celeaux et. al. 2006), where the 
#'			conditional likehihood evaluated at the posterior mode is replaced by the marginal predictive density. 
#'			Lastly, the random and fixed effects design matrices, \code{X, Z}, are returned that include 
#'			both the user input nuisance covariates appended to the time and treatment-based covariates constructed 
#'			by \code{dpgrow}.}  
#'	\item{print(summary(res))}{ prints contents of summary to console.}
#'      \item{plot(res)}{ returns results plots, including the set of subject random effects values and credible intervals, a sample
#'			of by-subject growth curves, mean growth curves split by each treatment and control, as well as selected trace plots
#'			for number of clusters and for precision parameters for the likehilood and random effects.  Lastly, a trace plot
#'			for the deviance statistic is also included.}
#'  	\item{samples(res)}{ contains (\code{n.iter - n.burn}) posterior sampling iterations for every model parameter, including fixed and random
#'			effects.}
#'	\item{resid(res)}{ contains the model residuals.}
#' @note The intended focus for this package are data where both number of subjects and number of repeated measures are limited.  
#'	This function places a DDP prior on the set of subject effects.  This means that the unnknown (random) prior on subject effects
#'	is indexed by the subject dosage patterns across one or more treatments.  
#' @keywords model
#' @seealso \code{\link{dpgrowmult}}, \code{\link{dpgrowmm}}, \code{\link{dpgrow}}
#' @examples 
#' \dontrun{
#' ## extract simulated dataset
#' library(growcurves)
#' data(datddpsim)
#' ## attach(datddpsim)
#' ## run dpgrow mixed effects model; returns object of class "ddpgrow"
#' shape.dp	= 4
#' res		= ddpgrow(y = dat$y, subject = dat$subject, 
#'			trt = dat$trt, time = dat$time,
#'			typetreat = c("mvn","car","ind","car"), 
#'			numdose = dat$numdose, 
#'			labt = dat$labt, dosemat = dat$dosemat, 
#'			Omega = dat$Omega, n.random = dat$n.random, 
#'			n.fix_degree = 2, n.iter = 10000, n.burn = 2000, 
#'			n.thin = 10, shape.dp = 1)
#' plot.results	= plot(res) ## ggplot2 plot objects, including growth curves
#' summary.results = summary(res) ## parameter credible intervals,  fit statistics
#' samples.posterior = samples(res) ## posterior sampled values
#' }
#' @aliases ddpgrow
#' @aliases ddpgrow.default
#' @author Terrance Savitsky \email{tds151@@gmail.com} Susan Paddock \email{paddock@@rand.org}
#' @references 
#'	T. D. Savitsky and S. M. Paddock (2011) Bayesian Hierarchical Semiparametric Modeling of Longitudinal Post-treatment Outcomes from Open-enrollment Therapy Groups, submitted to: JRSS Series A (Statistics in Society).
#' @references
#' 	T. D. Savitsky and S. M. Paddock (2011) Visual Sufficient Statistics for Repeated Measures data with growcurves for R, submitted to: Journal of Statistical Software.
#' @export ddpgrow 
ddpgrow			<- function(y, subject, trt, time, n.random, n.fix_degree, formula, random.only, data, dosemat, numdose, typetreat,
					labt, Omega, n.iter, n.burn, n.thin, shape.dp, rate.dp, M.init, plot.out)
					UseMethod("ddpgrow")

################################################
## default dispatch method for mm-session models
################################################
ddpgrow.default		<- function(y = NULL, subject, trt = NULL, time = NULL, n.random = NULL, n.fix_degree = NULL, formula = NULL, random.only = FALSE, data = NULL,
					dosemat, numdose, typetreat = NULL, labt = NULL, Omega = NULL, n.iter, n.burn, n.thin = 1,
					shape.dp = 1, rate.dp = 1, M.init = NULL, plot.out = TRUE)
{ ## start function dpgrow.default

  ############################
  ## check inputs
  ############################
  ## model choices
  if( is.null(typetreat) )
  {
	typetreat = rep("mvn",length(numdose))
  }else{
	typetreat = tolower(typetreat)
	if( length(typetreat) != length(numdose) ) stop("Length of vector 'typetreat' must equal length of 'numdose'.")
  	if( !("car" %in% typetreat) & !("mvn" %in% typetreat) & !("ind" %in% typetreat) )
  	{
		stop("You must pick from among the following 3 'typetreat' options for the treatments base measures, c('car','mvn','ind')")
  	}
  }

  # data choices
  if( is.null(dosemat) )
  {
	stop("Must input P x (sum(numdose)+1) 'dosemat' matrix to map effects to subjects.")
  }else{
	if( nrow(dosemat) != length(unique(subject)) ) stop("Number of rows of dosemat must equal total number of unique subjects.")
	if( !identical(dosemat[,1],rep(1,nrow(dosemat))) ) ## first column not an intercept
	{
		if( ncol(dosemat) == sum(numdose) ) ## left off intercept
		{
			dosemat = cbind(rep(1,nrow(dosemat)),dosemat)
			warning("Added an intercept column to 'dosemat' because it appears to have been left off.")
		}
	}
	if( (ncol(dosemat) - 1) != sum(numdose) ) stop("Number of columns of dosemat should be sum(numdose) + 1 -- total dosages + intercept.")
  }

  if( is.null(labt) ) ## if no user treatment labels entered, fill label vector as numeric sequence to use for plotting.
  {
	labt = 1:length(typetreat)
  }else{
	if( length(labt) != length(typetreat) ) stop("Length of treatment labels 'labt' must equal length of vector 'typetreat'.")
  }

  if( is.null(subject) ) stop("must input 'subject' vector that links subjects to cases (of length equal to the number of cases)")

  if( is.null(M.init) ) M.init = 10 ## avoid excessive start-up run time


  if(is.null(n.fix_degree)) 
  {
	if( !is.null(time) ) ## user wants growth curve 
	{
		n.fix_degree = length(unique(time)) - 1
		warning("Since 'n.fix_degree' not input, assumed it is equal to maximum number of unique values in 'time' to generate fixed effects.")
	}
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
  	if( not2part.test == TRUE & is.null(random.only) )
  	{
		stop("The formula is only 1 part - either fixed or random effects - but not both, so must input a boolean value for 'random.only'")
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

  ## trtcov - treatment covariates for anova random effects - need a numeric set for ddpeffectsplot
  start		<- 1
  out		<- relabel(label.input = labt, start)
  map.trtcov	<- out$dat.label

  ## create numt, typet as numeric for C++ runs
  numt 				<- numdose ## just a shorter name
  typet 			<- vector("numeric", length(numt))
  typet[typetreat == "car"] 	<- 1
  typet[typetreat == "mvn"] 	<- 2
  typet[typetreat == "ind"] 	<- 3

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
  if("car" %in% typetreat) ## create omega.plus (from Omega)
  {
	ncar		<- length(typet[typet == 1])
	omega.plus 	<- vector("list", ncar)
	for(m in 1:ncar)
	{
		omega.plus[[m]]	= rowSums(Omega[[m]])
	}
  }else{ ## !("car" %in% typetreat), so create dummy Omega and omega.plus
 	Omegamat 	= matrix(0, 2, 2)
	omegaplus	= matrix(0, 1, 2)
	Omega		= list(Omega = Omegamat)
	omega.plus	= list(omega.plus = omegaplus)
  }
  print(paste("Your chosen set of treatment base distributions  = ", paste(typetreat,collapse=" "), sep = ""))
  res 		= ddpPost(y, X, Z, subject, dosemat, numt, typet, Omega, omega.plus, n.iter, n.burn, n.thin, shape.dp, rate.dp, M.init)
 
  ##################################################################
  ## summary (short-hand) results
  ##################################################################
  summary.results			<- ddp_quantiles(model.output = res, dosemat = dosemat, Nfixed = Nfixed, Nrandom = Nrandom, Nsubject = Nsubject, typet = typet)
  summary.results$X			<- X
  summary.results$Z			<- Z
  summary.results$map.subject		<- map.subject
  summary.results$time			<- time ## not used in accessor functions; just reporting back to user to let them know that sorted by subject
  summary.results$map.trt		<- map.trt
  summary.results$map.trtcov		<- map.trtcov
  summary.results$typet			<- typet
  summary.results$numt			<- numt
  summary.results$n.fix_degree		<- n.fix_degree

  residuals				= colMeans(res$Residuals)

 if( !is.null(time) & length(unique(time)) > 1 )
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
   	gc.plot		= growthCurve(y.case = y, B = res$Theta, Alpha = res$Alpha, Beta = res$Beta, trt.case = trt, trt.lab = trti.u, subject.case = subject, 
				subject.lab = subjecti.u, T = T, min.T = min.T, max.T = max.T, n.thin = n.thin.gc, time.case = time, n.fix_degree = n.fix_degree)
   }else{ ## other fixed effects besides time-based covariates. Note: Either X.n or Z.n may be NULL (but not both), which is handled in the growthCurve function
	gc.plot		= growthCurve(y.case = y, B = res$Theta, Alpha = res$Alpha, Beta = res$Beta, X.n = X.n, Z.n = Z.n, 
				trt.case = trt, trt.lab = trti.u, subject.case = subject, subject.lab = subjecti.u, T = T, min.T = min.T, max.T = max.T, n.thin = n.thin, 
				n.waves = n.waves, time.case = time, n.fix_degree = n.fix_degree, Nrandom = Nrandom)
			## memo: if have nuisance covariates, need input of Nrandom to construct time-based random effects since Nrandom > n.random
   }

 } ## end conditional statement on creating growth curves

 
 if(plot.out == TRUE)
 {
   ##################################################################
   ## plots
   ################################################################## 
   ## memo: if(!("car" %in% typet)) then ddp_quantiles excludes alphacar.summary and taucar.summary, which means is.null(summary.results$taucar.summary) == TRUE
   ## and nrow(res$CAR_Q[[1]]) == 0.  (Note: is.null(CAR_Q) == FALSE as it was a defined field object in ddp.cpp).
   ## similarly, if( !("mvn" %in% typet) ) then ddp_quantiles excludes pmvn.mean, so is.null(summary.results$pmvn.mean) == TRUE

   plot.results = ddpMCMCplots(subjecti.u = subjecti.u, labt = labt, typet = typet, numt = numt, theta.summary = summary.results$theta.summary, lambda.mean = summary.results$lambda.mean,
					pmvn.mean = summary.results$pmvn.mean, taucar.summary = summary.results$taucar.summary, alphacar.summary = summary.results$alphacar.summary,
					Taucar = res$CAR_Q[[2]], Alphacar = res$CAR_Q[[1]], tauind.summary = summary.results$tauind.summary, 
					Tauind = res$Tauind, M = res$M, Taue = res$Taue, Deviance = res$Deviance)
   	 
 } #end conditional statement on whether to plot 

   
 ##################################################################
 ## function output
 ##################################################################
 if( plot.out == TRUE  )
 {
    if( !is.null(time)	& length(unique(time)) > 1 ) ## growth curves are plotted from time-based covariates
    {
    	plot.results$p.gcall = gc.plot$p.gcall; plot.results$p.gcsel = gc.plot$p.gcsel
	resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, Theta = res$Theta, devres = res$devres, 
			Num = res$Num, M = res$M, S = res$optPartition[[3]], C = res$C, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]], bigSmin = res$bigSmin,
		      	Residuals = res$Residuals, Tau.e = res$Taue, Pmvn = res$Pmvn, Alphacar = res$CAR_Q[[1]], 
			Taucar = res$CAR_Q[[2]], Tauind = res$Tauind, Lambda = res$Lambda, DoseEffects = res$DoseEffects,
			summary.results = summary.results,  plot.results = plot.results, residuals = residuals, 
			dat.growthCurve = gc.plot$plot.dat, dat.gcdata = gc.plot$dat.data)
    }else{ ## is.null(time) == TRUE
    	resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, Theta = res$Theta, devres = res$devres, 
			Num = res$Num, M = res$M, S = res$optPartition[[3]], C = res$C, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]], bigSmin = res$bigSmin,
		      	Residuals = res$Residuals, Tau.e = res$Taue, Pmvn = res$Pmvn, Alphacar = res$CAR_Q[[1]], 
			Taucar = res$CAR_Q[[2]], Tauind = res$Tauind, Lambda = res$Lambda, DoseEffects = res$DoseEffects,
			summary.results = summary.results,  plot.results = plot.results, residuals = residuals)	
   } ## end conditional statement on whether is.null(time)
 }else{ ## plot.out = FALSE
    	resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, Theta = res$Theta, devres = res$devres, 
			Num = res$Num, M = res$M, S = res$optPartition[[3]], C = res$C, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]], bigSmin = res$bigSmin,
		      	Residuals = res$Residuals, Tau.e = res$Taue, Pmvn = res$Pmvn, Alphacar = res$CAR_Q[[1]], 
			Taucar = res$CAR_Q[[2]], Tauind = res$Tauind, Lambda = res$Lambda, DoseEffects = res$DoseEffects,
			summary.results = summary.results,  residuals = residuals)	
 } ## end conditional statement on plot.out


 if( !any(typet == 1) ) 
 {
	resot$Alphacar 		<- NULL
	resot$Taucar 		<- NULL
 }
 if( !any(typet == 2) ) {resot$Pmvn <- NULL}
 if( !any(typet == 3) ) {resot$Tauind <- NULL}
 resot		<- resot[!sapply(resot, is.null)]
 

 ##
 ## return list output for dpgrow.default()
 ##

 resot$call		<- match.call()
 resot$Nrandom     	<- ncol(resot$summary.results$Z)
 resot$Nsubject		<- length(unique(subject))
 resot$subject		<- unique(subjecti.u) ## will employ for labeling B with user input subject labels
 class(resot)		<- c("ddpgrow")
 return(resot)


} #end function ddpgrow.default()

#####################################################
## .Call statements to C++ functions
#####################################################
#' Run a Bayesian mixed effects model for by-subject random effects with DDP prior
#'
#' An internal function to \code{\link{ddpgrow}}
#'
#' @export
#' @aliases ddpPost
#' @param y An \emph{N x 1} response (of subject-measure cases)
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix.  Assumed grouped by \code{subjects}
#' @param subject An \emph{N x 1} set of subject identifiers
#' @param dosemat An \emph{P x T} Anova or Multiple Membership design matrix linking treatment dosages to subjects
#'		where \emph{T} is the total number dosages across all treatments + 1 for an intercept.
#'		This formulation assumes there is a hold-out dose for each treatment.  e.g. the null dosage.
#' @param numt A numeric vector of length equal to the number of treatments that contains the number of dosages for each treatment.
#' @param typet A numeric vector of length equal to the number of treatments that contains the base distribution for each treatment.
#'		\code{1 = "car"}, \code{2 = "mvn"}, \code{3 = "ind"}
#' @param Omega A list object of length equal to the number of treatments with \code{"car"} selected for base distribution.
#'		Each entry is an \code{numt[m] x numt[m]} numeric CAR adjacency matrix for the dosages of treatment \code{m}.
#' @param omegaplus A list object of length equal to the number of treatments under \code{"car"} containing numeric vectors
#'		that are rowSums of the corresponding matrix element in \code{Omega}.
#' @param n.iter The number of MCMC iterations
#' @param n.burn The number of MCMC burn-in iterations to discard
#' @param n.thin The step increment of MCMC samples to return
#' @param shapealph The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#' @param ratebeta The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. 
#' @param M.init Initial MCMC chain scalar value for number of by-subject clusters.  If excluded defaults to \code{length(unique(subject))}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{ddpgrow}}
ddpPost = function(y,X,Z,subject,dosemat,numt,typet,Omega,omegaplus,n.iter,n.burn,n.thin,shapealph,ratebeta,M.init) {
    stopifnot(nrow(X) == nrow(Z))
    stopifnot(length(y) == nrow(X))
    res <- .Call("DDP", y,X,Z,subject,dosemat,numt,typet,Omega,omegaplus,n.iter,n.burn,n.thin,shapealph,ratebeta,M.init, package = "growcurves")
} ## end function ddpPost


####################################
## accessor methods
####################################
#' S3 functions of dpgrow
#'
#' produces quantile summaries for model parameters
#'
#' @param object A \code{ddpgrow} object
#' @param ... Ignored
#' @export 
#' @method summary ddpgrow
#' @aliases summary.ddpgrow
summary.ddpgrow <- function(object,...)
{
  res 		<- list(call = object$call, summary.results = object$summary.results)
  class(res) 	<- "summary.ddpgrow"
  return(res)
}

#' Print summary statistics for sampled model parameters
#'
#' provides credible intervals of sampled parameters for 
#' \code{ddpgrow} object
#'
#' @param x A \code{ddpgrow} object
#' @param ... Ignored
#' @export 
#' @method print summary.ddpgrow
#' @noRd
print.summary.ddpgrow <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCredible Intervals and Fit Statistics\n")
  print(x$summary.results)
}

#' Produce samples of MCMC output
#'
#' provides posterior sampled values for every model parameter of a
#' \code{ddpgrow} object
#'
#' @param object A \code{ddpgrow} object
#' @param ... Ignored
#' @export samples ddpgrow
#' @aliases samples.ddpgrow
#' @method samples ddpgrow
#' @aliases samples.ddpgrow
samples.ddpgrow <- function(object,...)
{
  Theta				<- as.data.frame(object$Theta)
  names(Theta)			<- paste(rep(1:object$Nrandom, each = object$Nsubject), rep(object$subject, object$Nrandom), sep=".") ## 1.1, 1.2, ...., 1.299
  Beta				<- as.data.frame(object$Beta)
  names(Beta)			<- colnames(object$summary.results$X)

  typet  		<- object$summary.results$typet
  res 			<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, Theta = Theta, DoseEffects = object$DoseEffects,
					Residuals = object$Residuals, M = object$M, S = object$S, C = object$C, Num.per.cluster = object$Num, 
					bigSmin = object$bigSmin, phat = object$phat, ordscore = object$ordscore, Alphacar = object$Alphacar, Taucar = object$Taucar,
					Pmvn = object$Pmvn, Tauind = object$Tauind, Lambda = object$Lambda, Tau.e = object$Tau.e)
  if( !any(typet == 1) ) 
  {
	res$Alphacar 		<- NULL
	res$Taucar 		<- NULL
  }
  if( !any(typet == 2) ) {res$Pmvn <- NULL}
  if( !any(typet == 3) ) {res$Tauind <- NULL}
  res		<- res[!sapply(res, is.null)]

  if( !is.null(object$dat.growthCurve) ) ## Add growth curve data set if user chooses growth curve option
  {
	res$dat.growthCurve = object$dat.growthCurve
  }
 
  class(res) 	<- "samples.ddpgrow"
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
#'	\item{plot.results}{	\code{ggplot2} plot objects.  See \code{\link{ddpMCMCplots}}. }
#'	\item{dat.growcurve}{	A \code{data.frame} containing fields \code{c("fit","time","subject","trt")}
#'		with \code{P*T} rows, where \code{P} is the length of \code{subject} and \code{T = 10} are the number of in-subject
#'		predictions for each subject.  This object may be used to construct additional growth curves using - see \code{\link{growplot}}.} 
#'      \item{dat.gcdata}{	A \code{data.frame} containing fields  \code{c("fit","time","subject","trt")} with \code{N} rows, where \code{N} are the
#'		number of subject-time cases.  This object contains the actual data for all subjects used to co-plot with predicted growth curves.}
#' @method plot ddpgrow
#' @aliases plot.ddpgrow
plot.ddpgrow <- function(x, plot.out = TRUE, ...)
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
  class(res)		<- "plot.ddpgrow"
  return(res)
}