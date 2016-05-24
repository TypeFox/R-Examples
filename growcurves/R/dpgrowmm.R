###################################
## generic for mm-session models
###################################
#' @include growthCurve.R
#' @include XZcov.R
#' @include relabel.R
#' @include mcmcPlots.R
#' @include summary_quantiles.R
NULL

#' Bayesian semiparametric growth curve models with employment of multiple membership random effects.
#'
#' Employs a Dirichlet Process (DP) prior on the set of by-subject random effect parameters
#' under repeated waves of measurements to allow the number of random effect parameters specified per subject, \code{q},
#' to be equal to the number of measurement waves, \code{T}.  Random effects are grouped by subject and
#' all \code{q} parameters receive the DP prior.  An additional set of random effects are included that 
#' utilize a different grouping factor; e.g. treatment(s) exposure or dosage.  These additional random
#' effects are mapped back to subjects through a multiple membership weight matrix.
#'
#' @param y A univariate continuous response, specified as an \emph{N x 1} matrix or vector, where \code{N}
#'	captures the number of subject-time cases (repeated subject measures).  Data may reflect unequal
#'	number of measures per subject.  Missing occasions are left out as no \code{NA} values are allowed.
#' @param subject The objects on which repeated measures are conducted that serves as the random effects
#'	grouping factor.  Input as an \emph{N x 1} matrix or vector of subject-measure cases in either
#'	integer or character format; e.g. \code{(1,1,1,2,2,3,3,3,...,n,n,n)}, where \code{n} is the total
#'	number of subjects.
#' @param trt An integer or character matrix/vector of length \code{N} (number of cases) indicating treatment
#'	group assignments for each case.  May also be input as length \code{P} vector, where \code{P} is
#'	the number of unique subjects, indicating subject group assignment.  Multiple treatment groups
#'	are allowed and if the vector is entered as numeric, e.g. \code{(0,1,2,3,..)}, the lowest numbered
#'	group is taken as baseline (captured by global fixed effects).  If entered in character format,
#'	the first treatment entry is taken as baseline.  If the are no treatment (vs. control) groups,
#'	then this vector may be excluded (set to NULL).
#' @param time A univariate vector of length \code{N}, capturing the time points associated to each by-subject
#'	measure.  Mav leave blank if only one time point (no repeated measures).
#' @param n.random The desired number of time-indexed subject random effect terms, \code{q}.  Since a DP prior is used on subject effects,
#'	may be set equal to the number of measurement waves, \code{T}.  The \code{y, trt, time} vectors will together
#'	be used to create both fixed and random effect design matrices.  The random effects matrix will be of the 
#' 	the form, \code{(1, time, ... , time^(n.random - 1))} (grouped, by \code{subject}). 
#'	This formulation is a growth curve model that allows assessment of by-treatment effects and by-client growth curves.
#' @param n.fix_degree The desired polynomial order in time to use for generating time-based fix effects.
#'	The fixed effects matrix will be constructed as, 
#'	\code{(time, ..., time^(n.fix_degree), trt_1,time*trt_1, ... ,time^(n.fix_degree)*trt_l, trt_L,..., time^(n.fix_degree)*trt_L)}.
#'	This formulation is a growth curve model that allows assessment of by-treatment effects and by-client growth curves.
#'	If \code{is.null(n.fix_degree) | n.fix_degree == 0 & is.null(trt)} time-by-treatment fixed effects and growth curves are not generated.
#' @param formula Nuisance fixed and random effects may be entered in \code{formula} with the following format,
#'	\code{y ~ x_1 + x_2*x_3 | z_1*z_2 } as an object of class \code{formula}.  The bar, \code{|}, separates fixed and random effects.  If it
#'	is only desired to enter either fixed or random effects, but not both then the \code{|} may be omitted.  Note:
#'	the nuisance random effects are assumed to be grouped by subject.  The fixed and random effects values may change with
#'	each repeated measure; however, within subject growth curves will keep constant \code{z} and \code{x} values between
#'	measurement waves.   It is possible to bypass the growth curve construction by leaving \code{y, trt, time, n.random, n.fix_degree} 
#'	blank and entering only \code{formula}, instead.  The model output plots, will, however
#'	exclude growth curves in that event.  If a formula is input (which requires response, \code{y}) then
#'	the separate entry of \code{y} may be omitted.  If the parameter \code{y} is input, it will be over-written by that from \code{formula}.
#' @param random.only A Boolean variable indicating whether the input formula contains random (for fixed) effects in the case that only
#'	one set are entered.  If excluded and \code{formula} is entered without a \code{|}, \code{random.only} defaults to 
#'	\code{FALSE}.
#' @param data a \code{data.frame} containing the variables with names as specified in \code{formula}, including the response, \code{y}.
#' @param Omega An \emph{S x S} numerical matrix object to encode the CAR adjacency matrix for random effects mapped through multiple membership,
#'	where \code{S} is the number of effects mapped to subjects through the multiple membership construction.  
#'	This input applies only to \code{option = "mmcar"}.
#' @param group	A numeric or character vector of length \code{S}, providing group identifiers for each of \code{S} multiple membership effects
#'	(e.g. \code{(1,1,1,2,2,...)}.  If excluded, it is assumed there is a single group.
#' @param subj.aff A \emph{n.aff x 1} vector subset of \code{subject} composed with unique subject identifiers that are linked to the multiple
#'	membership effects; e.g. one or more treatment cohorts.  If all subjects are to receive the mapping of multiple membership effects,
#'	\code{subj.aff} may be left blank.
#' @param W.subj.aff An \emph{n.aff x S} numeric matrix that maps a set of random effects to affected subjects, where \code{P.aff} is the length
#'	of the unique subjects to whom the multiple membership random effects applies.  It is assumed that the row order is the same
#'	as the order of \code{subj.aff} (or \code{unique(subject)} if \code{subj.aff} is not input).  If \code{W.subj.aff} is a multiple membership
#'	weight matrix, then the rows will sum to 1.  The form and therefore, interpretation of output is dependent on form of input; for example, 
#'	the rows of \code{W.subj.aff} may include indicators for whether each of \code{S} treatment dosages are linked to a given \code{subject}.
#' @param multi A boolean scalar input that when set to \code{TRUE} indicates the each of the \code{S} MM effects is multivariate.
#'		Leave blank if univariate multiple membership effects are desired.
#'		It is assumed that the associated design matrix is equal to the number of time-indexed random effects.
#'		For example, the time-indexed (non-nuisance) random effects design matrix, Z = (1,time,time^{n.random-1}) is also
#'		used to compose an inner product with each row of the \code{N x n.random} MM product, \code{W * U}, where 
#'		\code{W} is case-expanded MM design matrix and \code{U} is the \code{S x n.random} set of multivariate MM effects.
#' @param n.iter Total number of MCMC iterations.
#' @param n.burn Number of MCMC iterations to discard.  \code{dpgrow} will return \code{(n.iter - n.burn)} posterior samples.
#' @param n.thin Gap between successive sampling iterations to save.
#' @param strength.mm Sets both the shape and rate parameter for a \code{tau_{gamma} ~ G(strength.mm,strength.mm)} prior on the precision parameter of either a
#'	CAR (\emph{gamma ~ CAR(tau_gamma)}) or independent (\emph{gamma ~ N(0,tau^(-1)I_S}) prior on the set of \code{S} 
#'	multiple membership effects.
#' @param shape.dp Shape parameter under a \emph{c ~ G(shape.dp, 1)} prior on the concentration parameter of the DP (prior
#'	on the set of random effects parameters, \emph{b_1, ..., b_n ~ DP(c,G_0)}
#'	where \code{n} is the total number of subjects.
#' @param rate.dp Rate parameter under a \emph{c ~ G(shape.dp, rate.dp)} prior on the concentration parameter of the DP.
#' @param plot.out A boolean variable indicating whether user wants to return plots with output results.  Defaults to \code{TRUE}.
#' @param option Modeling option, of which there are three: 1. \code{mmcar} places a CAR prior on the set of multiple membership effects;
#'				2. \code{mmi} places the usual independent Gaussian priors on the set of multiple membership effects.
#'				3. \code{mmigrp} employs a set of independent Gaussian priors, but with a common mean parameter
#'					for each sub-group of multiple membership effects sharing a common group identifier.
#'					(e.g. treatment groups that disjointly divide therapy sessions from Savitsky and Paddock (2011))
#' @return S3 \code{dpgrowmm} object, for which many methods are available to return and view results.  Generic functions applied
#'	to an object, \code{res} of class \code{dpgrow}, includes:
#'	\item{summary(res)}{returns \code{call}, the function call made to \code{dpgrowmm} and \code{summary.results}, which contains
#'					a list of objects that include \emph{95\%} credible intervals for each set of sampled parameters, 
#'					specified as (\code{2.5\%}, mean, \emph{97.5\%}, including fixed and random effects. 
#'					Also contains model fit statistics, including \code{DIC} (and associated \code{Dbar}, \code{Dhat}, \code{pD}, 
#'					\code{pV}), as well as the log pseudo marginal likelihood (LPML), a leave-one-out fit statistic.  
#'					Note that \code{DIC} is constructed as \code{DIC3} (see Celeaux et. al. 2006), where the 
#'					conditional likehihood evaluated at the posterior mode is replaced by the marginal predictive density. 
#'					Lastly, the random and fixed effects design matrices, \code{X, Z}, are returned that include 
#'					both the user input nuisance covariates appended to the time and treatment-based covariates constructed 
#'					by \code{dpgrowmm}.}  
#'	\item{print(summary(res))}{prints contents of summary to console.}
#'      \item{plot(res)}{returns results plots, including the set of subject random effects values and credible intervals, a sample
#'					of by-subject growth curves, mean growth curves split by each treatment and control, as well as selected trace plots
#'					for number of clusters and for precision parameters for the likehilood and random effects.  Lastly, a trace plot
#'					for the deviance statistic is also included.}
#'  	\item{samples(res)}{contains (\code{n.iter - n.burn}) posterior sampling iterations for every model parameter, including fixed and random
#'					effects.}
#'	\item{resid(res)}{contains the model residuals.}
#' @note The intended focus for this package are data where both number of subjects and number of repeated measures are limited.  A DP prior
#'	is placed on the by-subject random effects to borrow strength across subjects for each estimation of each subject's growth curve.  The
#'	imposition of the DP prior also allows the resulting posterior distributions over the subject random effects to be non-Gaussian.
#'	The \code{dpgrow} function is very similar to \code{dpgrowmm}; only the latter includes a separate set of random effects not grouped
#'	by subject (e.g. for treatment dosages allocated to subjects) mapped back to subject-time cases through a multiple membership design matrix. 
#'	The \code{dpgrowmult} function generalizes \code{dpgrowmm} by allowing more than one multiple membership effects term. 
#'	See Savitsky and Paddock (2011) for detailed model constructions.
#' @keywords model
#' @seealso \code{\link{dpgrow}}, \code{\link{dpgrowmult}}
#' @examples 
#' \dontrun{
#' ## extract simulated dataset
#' library(growcurves)
#' data(datsim)
#' ## attach(datsim)
#' ## Model with DP on clients effects, but now INCLUDE session random effects
#' ## in a multiple membership construction communicated with the N x S matrix, W.subj.aff.
#' ## Returns object, res.mm, of class "dpgrowmm".
#' shape.dp	= 3
#' strength.mm	= 0.001
#' res.mm	= dpgrowmm(y = datsim$y, subject = datsim$subject, 
#'			trt = datsim$trt, time = datsim$time, 
#'			n.random = datsim$n.random, 
#'			n.fix_degree = 2, Omega = datsim$Omega, 
#'			group = datsim$group, 
#'			subj.aff = datsim$subj.aff,
#'			W.subj.aff = datsim$W.subj.aff, 
#'			n.iter = 10000, n.burn = 2000, n.thin = 10,
#'			shape.dp = shape.dp, rate.dp = rate.dp, 
#'			strength.mm = strength.mm, option = "mmcar")
#' plot.results		= plot(res.mm) ## ggplot2 plot objects,
#' summary.results	= summary(res.mm) ## credible intervals and fit statistics
#' samples.posterior	= samples(res.mm) ## posterior sampled values
#' }
#' @aliases dpgrowmm
#' @aliases dpgrowmm.default
#' @author Terrance Savitsky \email{tds151@@gmail.com} Susan Paddock \email{paddock@@rand.org}
#' @references 
#'	S. M. Paddock and T. D. Savitsky (2012) Bayesian Hierarchical Semiparametric Modeling of Longitudinal Post-treatment Outcomes from Open-enrollment Therapy Groups, submitted to: JRSS Series A (Statistics in Society).
#' @references
#' 	T. D. Savitsky and S. M. Paddock (2012) Visual Sufficient Statistics for Repeated Measures data with growcurves for R, submitted to: Journal of Statistical Software.
#' @export dpgrowmm 
dpgrowmm			<- function(y, subject, trt, time, n.random, n.fix_degree, formula, random.only, data, Omega, group, subj.aff, W.subj.aff, multi, n.iter, 
					n.burn, n.thin, strength.mm, shape.dp, rate.dp, plot.out, option)
					UseMethod("dpgrowmm")

################################################
## default dispatch method for mm-session models
################################################
dpgrowmm.default		<- function(y = NULL, subject, trt = NULL, time = NULL, n.random = NULL, n.fix_degree = NULL, formula = NULL, 
						random.only = FALSE, data = NULL, Omega = NULL, group = NULL, subj.aff = NULL, W.subj.aff, multi = FALSE,
					 	n.iter, n.burn, n.thin = 1, strength.mm = 0.1, shape.dp = 1, rate.dp = 1, plot.out = TRUE, option = "mmi")
{ ## start function dpgrowmm.default

  ############################
  ## check inputs
  ############################
  ## model choices
  if( option != "mmcar" & option != "mmi" & option != "mmigrp")
  {
	stop("You must pick 1 of 3 modeling options, c('mmcar','mmi','mmigrp')")
  }

  if( is.null(Omega) & option == "mmcar" ) ## inconsistencies
  {
	stop("You must specify the (S = number of effects) x S adjacency matrix, Omega with option = 'mmcar'.")
  }

  ## data choices
  if( is.null(W.subj.aff) ) 
  {
  	stop("Must input multiple membership design matrix, 'W.subj.aff'. ")
  }else{ ## !is.null(W.subj.aff)
	if( !is.null(subj.aff) )
	{
		if( length(setdiff(subj.aff,subject)) != 0 ) stop("Input vector, subj.aff, must contain same subject name format as input vector, subject")
		if( nrow(W.subj.aff) != length(subj.aff) ) stop("Number of rows of W.subj.aff must equal length of subj.aff, the affected clients")
		if( length(unique(subj.aff)) != length(subj.aff) ) ## already know that nrow(W.subj.aff) == length(subj.aff)
  		{
			warning("Vector 'subj.aff' should contain number of unique subjects, not subject-time cases.  Function will shrink to unique values.")
			## shrink W.subj.aff (within subject) to unique rows
			tmp		= as.data.frame(cbind(subj.aff,W.subj.aff))
			tmp		= unique(tmp)
			W.subj.aff	= tmp[,-1]
			rm(tmp)
			## shrink subj.aff to unique values
			subj.aff 	= unique(subj.aff)  ## ensure of length Nsubject, not Ncase
  		}
	}else{ ## subj.aff is NULL, so mm random effects apply to all clients
		subj.aff 	= unique(subject) ## all subjects are affected
		if( nrow(W.subj.aff) == length(subject) ) ## input in case format, rather than subject
  		{
			warning("Rows of W.subj.aff should be equal to length of affected unique subjects, not number of subject-time cases.  
				Function will use unique ID's in subject vector to shrink rows of W.subj.aff.")
			tmp		= as.data.frame(cbind(subject,W.subj.aff))
			tmp		= unique(tmp)
			W.subj.aff	= tmp[,-1]
			rm(tmp)
  		}else{ ## W.subj.aff is not in case format and is.null(subj.aff)
			if( nrow(W.subj.aff) != length(unique(subject)) ) stop("If leave subj.aff blank, then number of rows of W.subj.aff must equal length of unique subjects.")
		}
	}
  }

  if ( is.null(group) ) ## always ensure a value for 'group', even for option = "mmi"
  {
	group = matrix(1,ncol(W.subj.aff),1) ## set equal to 1's if no grouping
  }


  if( option == "mmcar" )
  {
	if( ncol(W.subj.aff) != ncol(Omega) ) stop("Number of columns of 'W.sub.aff' and 'Omega' must be equal (to the number of mult-mbrship random effects)")
	if( ncol(Omega) != nrow(Omega) ) stop("Omega must be a square matrix dimensioned by number of multi-mbrship random effects")
	if( ncol(Omega) != length(group) ) stop("The length of 'group' must equal the row and column dimensions of Omega - number of mult- mbrship random effects")
  }

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

	if( !is.null(data) ){ warning("If you want to include supplemental or nuisance predictors in the 'data' input, must also include the 'formula' input to tell the model their structure.")}
  }

  if( is.null(data) & is.null(time) )
  {
	stop("Input data must be supplied to run model; e.g. (subject,time,trt,n.random) for growth curve and/or 'data' for nuisance covariates.")
  }


  if(  any(is.na(subject)) | any(is.na(W.subj.aff)) )
  {
	stop("NA's aren't allowed in any data cells - subject or W.subj.aff")
  }

  if( !is.null(time) )
  {
	if( any(is.na(time)) | any(is.na(trt)) ) stop("No NA's allowed in c(time,trt) vectors")
        if(length(time) != length(subject)) stop("Must input vector time in subject-time case format")
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


  ####################################################################################
  ## re-cast subject identifier inputs to be sequential - subject, subj.aff, trt, group
  ####################################################################################
  ## subject
  start			<- 1
  out			<- relabel(label.input = subject, start)
  subject		<- out$label.new
  o			<- order(subject) ## use later to place X, Z, map.subject, map.trt in contiguous order of subject
  subjecti.u		<- out$labeli.u
  map.subject		<- out$dat.label ## colnames = c("label.input","label.new"), in case format

  ## subj.aff - capture as strict subset of subject.  If user doesn't enter, subj.aff set equal to unique(subject) during user input testing, above.
  ## don't re-order subj.aff b/c correponds to rows in W.subj.aff.  only change labels
  subjaff.input 	<- subj.aff
  o.aff			<- 1:length(subj.aff) ## will use to maintain order of input
  tmp			<- data.frame(o.aff,subjaff.input)
  names(tmp)		<- c("order","label.input")
  smap.u		<- unique(map.subject) ## in subject format
  smap.u		<- subset(smap.u, smap.u$label.input %in% subjaff.input)
  tmp			<- merge(tmp,smap.u, by = "label.input", all.x = T, sort =  FALSE)
  subj.aff		<- tmp$label.new[order(tmp$order)] ## result is a strict subset of subject, but in subject, not case format
  rm(tmp,smap.u)

  
  ## trt
  start		<- 0
  out		<- relabel(label.input = trt, start)
  trt		<- out$label.new
  trti.u	<- out$labeli.u
  map.trt	<- out$dat.label

  ## group
  start		<- 1
  out		<- relabel(label.input = group, start)
  group		<- out$label.new
  groupi.u	<- out$labeli.u
  map.grp	<- out$dat.label
  

  #################################################################
  ## some subject, session, case lengths for use in subsetting and loops
  #################################################################
  Ncase			= length(subject)
  Nsubject 		= length(unique(subject))
  Nsession		= ncol(W.subj.aff)  ## this is true for both univariate and multivariate session effects
  Nsubj.aff		= length(subj.aff)
  Nlevel		= length(unique(trt))
  if(is.null(group)) group = matrix(1,length(subject))
  G 			= length(unique(group))
  iter.keep		= floor( (n.iter - n.burn)/n.thin )
  if(is.null(n.random)) n.random = min( length(unique(time)),4 ) ## max number of random effects is q = 4, which produces global cubic fit
  if(!is.null(time)) n.waves = length(unique(time)) ## number of measurement waves - used for growth curve generation with nuisance covariates
 
  ##################################################################
  ## construct fixed and random effect design matrices - ordered by subject so that subject is contiguous
  ##################################################################
  out 	<- XZcov(time = time , trt = trt, trt.lab = trti.u, subject = subject, n.random = n.random, n.fix_degree = n.fix_degree, formula = formula, 
		random.only = random.only, data = data)
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
  Nrandom		= ncol(Z) ## counts both number of time-indexed and nuisance random effects per subject

  ## set design matrix for multivariate mm effects equal to time-indexed random effects design matrix
  if( multi == TRUE)
  {
	Nmv	<- n.random
	H 	<- Z.c
  }else{ ## univariate MM effects
	Nmv	<- 1
  }

  #################################################################
  ## re-cast inputs to matrices
  #################################################################
  ## Expand W.subj.aff of rows subj.aff to Ncase rows, W.case
  if( length(subj.aff) < length(unique(subject)) )
  {
  	W.subj				= matrix(0, nrow = Nsubject, ncol = Nsession)
  	## subj.u			= unique(subject) # of length total number of subjects
  	## subj.ntrt			= setdiff(subj.u,subj.aff)
  	W.subj[subj.aff,]			= W.subj.aff
  	## W.subj[subj.ntrt,]		= 0
  }else{ ## all subjects should be linked to multiple membership affects
	W.subj				= W.subj.aff  ## may be input by user as a data.frame object
	W.subj				= as.matrix(W.subj) ## need W.subj to be a matrix object in growthCurves to differentiate from dpgrowmult object, where it is a list
  } ## end conditional statement on whether exists subj.aff subset of subject

  W.case			= W.subj[subject,]
  W.case			= as.matrix(W.case) 
  W.subj.aff			= as.matrix(W.subj.aff)
  y				= as.matrix(y,Ncase,1)

  ################################################################
  ## conduct posterior sampling and capture results
  ################################################################
  if( Nmv == 1 ) ## univariate MM effects
  {
  	if(option == "mmcar")
  	{
		print("Your chosen option = mmcar")
		omega.plus	= rowSums(Omega)
  		res 		= mmCplusDpPost(y, X, Z, W.case, W.subj.aff, Omega, omega.plus, group, subject, n.iter, n.burn, n.thin, strength.mm, shape.dp, rate.dp)
  	}else{
		if(option == "mmi")
		{
			print("Your chosen option = mmi")
			res 		= mmIplusDpPost(y, X, Z, W.case, W.subj.aff, subject, n.iter, n.burn, n.thin, strength.mm, shape.dp, rate.dp)
		}
		else{ ## option == "mmigrp"
			print("Your chosen option = mmigrp")
			M			= matrix(0,length(group),G) ## Design matrix for session mean parameters in MM(I)
			for(i in 1:G)
			{
				M[group == i,i]	= 1
			}
			res 		= mmIgroupDpPost(y, X, Z, W.case, W.subj.aff, M, subject, n.iter, n.burn, n.thin, strength.mm, shape.dp, rate.dp)
		}				
  	}
  }else{ ## multivariate MM effects
	if( option %in% c("mmi","mmigrp") )
	{
		print("Your chosen option = mmi for multivariate MM effects")
		Omega	 	= matrix(0, Nsession, Nsession)
		omega.plus	= matrix(0, Nsession, 1)
		corsess		= 0 ## correlations in prior scale matrix for wishart prior on precision matrix for effects order
		typemm		= 0 ## "mmi"
		res 		= mmCmvplusDpPost(y, X, Z, H, W.case, W.subj.aff, Omega, omega.plus, group, subject, n.iter, n.burn, n.thin, strength.mm, corsess, shape.dp, rate.dp, typemm)
		stopifnot( ncol(res$U) == (Nmv*Nsession) )
	}else{ ## option == "mmcar"
		print("Your chosen option = mmcar for multivariate MM effects")
		omega.plus	= rowSums(Omega)
		corsess		= 0 ## correlations in prior scale matrix for wishart prior on precision matrix for effects order
		typemm		= 1 ## "mmcar"
		res 		= mmCmvplusDpPost(y, X, Z, H, W.case, W.subj.aff, Omega, omega.plus, group, subject, n.iter, n.burn, n.thin, strength.mm, corsess, shape.dp, rate.dp, typemm)
		stopifnot( ncol(res$U) == (Nmv*Nsession) )
	}

  }


  ##################################################################
  ## summary (short-hand) results
  ##################################################################
  summary.results			<- summary_quantiles(model.output = res, Nfixed = Nfixed, Nrandom = Nrandom, Nsubject = Nsubject, Nsubj.aff = Nsubj.aff, Nmv = Nmv, Nsession = Nsession)
  if( Nmv == 1 ) {summary.results$rhotauu.summary <- NULL} ## this will remove the element from the list object that is always returned from 'summary.results'
  summary.results$X			<- X
  summary.results$Z			<- Z
  summary.results$map.subject		<- map.subject
  summary.results$time			<- time ## not used in accessor functions; just reporting back to user to let them know that sorted by subject
  summary.results$map.trt		<- map.trt
  summary.results$map.grp		<- map.grp
  summary.results$model 		<- option
  summary.results$n.fix_degree		<- n.fix_degree
  summary.results$Nmv			<- Nmv

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
  
   if( is.null(X.n) & is.null(Z.n) ) ## no nuisance covariates
   {
   	gc.plot		= growthCurve(y.case = y, B = res$B, Alpha = res$Alpha, Beta = res$Beta, U = res$U, aff.clients = subj.aff, W.subj = W.subj, 
					trt.case = trt, trt.lab = trti.u, subject.case = subject, subject.lab = subjecti.u,
					T = T, min.T = min.T, max.T = max.T, n.thin = n.thin.gc, time.case = time, n.fix_degree = n.fix_degree)
   }else{ ## other fixed effects besides time-based covariates. Note: Either X.n or Z.n may be NULL (but not both), which is handled in the growthCurve function
	gc.plot		= growthCurve(y.case = y, B = res$B, Alpha = res$Alpha, Beta = res$Beta, U = res$U, aff.clients = subj.aff, W.subj = W.subj, X.n = X.n, Z.n = Z.n, 
				trt.case = trt, trt.lab = trti.u, subject.case = subject, subject.lab = subjecti.u, T = T, min.T = min.T, max.T = max.T, n.thin = n.thin, 
				n.waves = n.waves, time.case = time, n.fix_degree = n.fix_degree, Nrandom = n.random)
			## memo: if have nuisance covariates, need input of Nrandom = n.random to construct time-based random effects since Q > Nrandom
   }

 } ## end conditional statement on creating growth curves

   
 if(plot.out == TRUE)
 {
   ##################################################################
   ## plots
   ##################################################################

   plot.results = mcmcPlots(subjecti.u = subjecti.u, subj.aff = subj.aff, subjaff.input = subjaff.input, bmat.summary = summary.results$bmat.summary, group = group, 
				groupi.u = groupi.u, u.summary = summary.results$u.summary, Nmv = Nmv, mm.summary = summary.results$mm.summary, 
				M = res$M, Tauu = res$Tauu, Taub = res$Taub, Taue = res$Taue, Deviance = res$Deviance)

 } #end conditional statement on whether to plot 

   
 ##################################################################
 ## function output
 ##################################################################
 if(plot.out == TRUE )
 {
   if( (!is.null(time) & length(unique(time)) > 1) & !is.null(n.fix_degree) ) ## a set of growth curves were generated from time-based covariates
   {
   	plot.results$p.gcall = gc.plot$p.gcall; plot.results$p.gcsel = gc.plot$p.gcsel
   	if( (option != "mmcar") & (option != "mmi") & (multi == FALSE) )
   	{
		
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, Eta = res$Eta,
			M = res$M, S = res$optPartition[[3]], Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, plot.results = plot.results,
			residuals = residuals, dat.growthCurve = gc.plot$plot.dat, dat.gcdata = gc.plot$dat.data) ##, Tau.eta = res$Taueta
   	}else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, 
			M = res$M, S = res$optPartition[[3]], Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, plot.results = plot.results,
			residuals = residuals, dat.growthCurve = gc.plot$plot.dat, dat.gcdata = gc.plot$dat.data)
   	}
   }else{ ## is.null(time) = TRUE
     	if((option != "mmcar") & (option != "mmi") & (multi == FALSE) )
     	{
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, Eta = res$Eta,
			M = res$M, S = res$optPartition[[3]], Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, Tau.eta = res$Taueta, summary.results = summary.results, plot.results = plot.results, residuals = residuals)
     	}else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, 
			M = res$M, S = res$optPartition[[3]], Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, plot.results = plot.results, residuals = residuals)
   	}  ## end conditional statement on choice	
   } ## end conditional statement on whether is.null(time)
 }else{ ## plot.out = FALSE
   if((option != "mmcar") & (option != "mmi"))
   {
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, Eta = res$Eta,
			M = res$M, S = res$optPartition[[3]], Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, 
			Tau.b = res$Taub, Tau.eta = res$Taueta, summary.results = summary.results, residuals = residuals)
   }else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, 
			M = res$M, S = res$optPartition[[3]], Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, residuals = residuals)
   }
 } ## end conditional statement on plot.out

 ##
 ##  Add unique return item for multivariate objects
 ##
 
 if( Nmv > 1 ) { resot$Rhotau.u = res$Rhotauu }

 ##
 ## return list output for dpgrowmm.default()
 ##

 resot$call		<- match.call()
 resot$Nrandom     	<- ncol(resot$summary.results$Z)
 resot$Nsubject		<- length(unique(subject))
 resot$subject		<- unique(subjecti.u) ## will employ for labeling B with user input subject labels
 class(resot)		<- c("dpgrowmm")
 return(resot)


} ### end function dpgrowmm.default()


#####################################################
## .Call statements to C++ functions
#####################################################
#' Bayesian mixed effects model with a DP prior on by-subject effects and CAR prior on a set of multiple membership effects
#'
#' An internal function to \code{\link{dpgrowmm}}
#'
#' @export mmCplusDpPost 
#' @aliases mmCplusDpPost mmC
#' @param y An \emph{N x 1} response (of subject-measure cases)
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix.  Assumed grouped by \code{subjects}
#' @param Wcase An \emph{N x 1} multiple membership weight matrix to map supplemental random effects
#' @param Wsubject An \emph{P.aff x S} multiple membership weight matrix with rows equal to number of unique affected subjects
#' @param Omega An \emph{S x S} unnormalized adjacency matrix with entries equal to 1 where two effects communicate
#'	and 0, otherwise.  Diagonal elements are zero
#' @param omegaplus \emph{S x 1} vector of row sums of \code{Omega}
#' @param groups \emph{S x 1} vector of group identifiers for each effect.  Effects within each group communicate.
#' 	Effects don't communicate across groups.
#' @param subjects An \emph{N x 1} set of subject identifiers
#' @param niter The number of MCMC iterations
#' @param nburn The number of MCMC burn-in iterations to discard
#' @param nthin The step increment of MCMC samples to return
#' @param strength.mm The shape and rate parameters for the \eqn{\Gamma} prior on the CAR precision parameter, \eqn{\tau_{\gamma}}
#' @param shapealph The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'	The rate parameter is set of \code{1}.
#' @param ratebeta The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. Default value is \code{1}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrowmm}}
mmCplusDpPost = function (y, X, Z, Wcase, Wsubject, Omega, omegaplus, groups, subjects, niter, nburn, nthin, strength.mm, shapealph, ratebeta) {
    stopifnot(nrow(X) == nrow(Z))
    stopifnot(nrow(Wcase) == nrow(Z))
    stopifnot(length(y) == nrow(X))
    stopifnot(length(omegaplus) == nrow(Omega))
    res <- .Call("mmCplusDP", y, X, Z, Wcase, Wsubject, Omega, omegaplus, groups, subjects, niter, nburn, nthin, strength.mm, shapealph, ratebeta, package = "growcurves")
}

#' Bayesian mixed effects model with a DP prior on by-subject effects and use of group means for multiple membership effects
#'
#' An internal function to \code{\link{dpgrowmm}}
#'
#' @export mmIgroupDpPost
#' @aliases mmIgroupDpPost mmIgroup
#' @param y An \emph{N x 1} response (of subject-measure cases)
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix.  Assumed grouped by \code{subjects}
#' @param Wcase An \emph{N x 1} multiple membership weight matrix to map supplemental random effects
#' @param Wsubject An \emph{P.aff x S} multiple membership weight matrix with rows equal to number of unique affected subjects
#' @param M An \emph{S x G} design matrix mapping (G) group means to the multiple membership effects
#'     Posterior samples are centered on each iteration to identify the global mean parameter.
#' @param subjects An \emph{N x 1} set of subject identifiers
#' @param niter The number of MCMC iterations
#' @param nburn The number of MCMC burn-in iterations to discard
#' @param nthin The step increment of MCMC samples to return
#' @param strength.mm The shape and rate parameters for the \eqn{\Gamma} prior on the CAR precision parameter, \eqn{\tau_{\gamma}}
#' @param shapealph The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'	The rate parameter is set of \code{1}.
#' @param ratebeta The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. Default value is \code{1}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{dpgrowmm}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrowmm}}
mmIgroupDpPost = function (y, X, Z, Wcase, Wsubject, M, subjects, niter, nburn, nthin, strength.mm, shapealph, ratebeta) {
    stopifnot(nrow(X) == nrow(Z))
    stopifnot(nrow(Wcase) == nrow(Z))
    stopifnot(length(y) == nrow(X))
    res <- .Call("mmIgroupDP", y, X, Z, Wcase, Wsubject, M, subjects, niter, nburn, nthin, strength.mm, shapealph, ratebeta, package = "growcurves")
}

#' Bayesian mixed effects model with a DP prior on by-subject effects and zero mean independent Gaussian priors on multiple membership effects
#'
#' An internal function to \code{\link{dpgrowmm}}
#'
#' @export mmIplusDpPost
#' @aliases mmIplusDpPost mmI
#' @param y An \emph{N x 1} response (of subject-measure cases)
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix.  Assumed grouped by \code{subjects}
#' @param Wcase An \emph{N x 1} multiple membership weight matrix to map supplemental random effects
#' @param Wsubject An \emph{P.aff x S} multiple membership weight matrix with rows equal to number of unique affected subjects
#' @param subjects An \emph{N x 1} set of subject identifiers
#' @param niter The number of MCMC iterations
#' @param nburn The number of MCMC burn-in iterations to discard
#' @param nthin The step increment of MCMC samples to return
#' @param strength.mm The shape and rate parameters for the \eqn{\Gamma} prior on the CAR precision parameter, \eqn{\tau_{\gamma}}
#' @param shapealph The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'	The rate parameter is set of \code{1}.
#' @param ratebeta The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. Default value is \code{1}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{dpgrowmm}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrowmm}}
mmIplusDpPost = function (y, X, Z, Wcase, Wsubject, subjects, niter, nburn, nthin, strength.mm, shapealph, ratebeta) {
    stopifnot(nrow(X) == nrow(Z))
    stopifnot(nrow(Wcase) == nrow(Z))
    stopifnot(length(y) == nrow(X))
    res <- .Call("mmIplusDP", y, X, Z, Wcase, Wsubject, subjects, niter, nburn, nthin, strength.mm, shapealph, ratebeta, package = "growcurves")
}

#' Bayesian mixed effects model with a DP prior on by-subject effects and CAR prior on a multivariate set of multiple membership effects
#'
#' An internal function to \code{\link{dpgrowmm}}
#'
#' @export mmCmvplusDpPost 
#' @aliases mmCmvplusDpPost mmCmv
#' @param y An \emph{N x 1} response (of subject-measure cases)
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix.  Assumed grouped by \code{subjects}
#' @param H Multivariate MM effects design matrix.
#' @param Wcase An \emph{N x 1} multiple membership weight matrix to map supplemental random effects
#' @param Wsubject An \emph{P.aff x S} multiple membership weight matrix with rows equal to number of unique affected subjects
#' @param Omega An \emph{S x S} unnormalized adjacency matrix with entries equal to 1 where two effects communicate
#'	and 0, otherwise.  Diagonal elements are zero
#' @param omegaplus \emph{S x 1} vector of row sums of \code{Omega}
#' @param groups \emph{S x 1} vector of group identifiers for each effect.  Effects within each group communicate.
#' 	Effects don't communicate across groups.  Not used under "mmi" prior (though input is required).
#' @param subjects An \emph{N x 1} set of subject identifiers
#' @param niter The number of MCMC iterations
#' @param nburn The number of MCMC burn-in iterations to discard
#' @param nthin The step increment of MCMC samples to return
#' @param strength.mm The shape and rate parameters for the \eqn{\Gamma} prior on the CAR precision parameter, \eqn{\tau_{\gamma}}.
#' @param corsess A single value to set the prior correlations among the multivariate \code{q = ncol(H)} orders for the MM effects.
#'		where \eqn{\tau_{\gamma}} is replaced by the \code{q x q}, \eqn{\Lambda}.
#' @param shapealph The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'	The rate parameter is set of \code{1}.
#' @param ratebeta The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. Default value is \code{1}.
#' @param typemm An indicator the prior formulation specified for the multivariate MM effects term.
#'		Set \code{typemm = 0} for \code{"mmi"} and \code{typemm = 1} for \code{"mmcar"}.
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{dpgrow}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrowmm}}
mmCmvplusDpPost = function (y, X, Z, H, Wcase, Wsubject, Omega, omegaplus, groups, subjects, niter, nburn, nthin, strength.mm, corsess, shapealph, ratebeta, typemm) {
    stopifnot(nrow(X) == nrow(Z))
    stopifnot(nrow(Wcase) == nrow(Z))
    stopifnot(length(y) == nrow(X))
    stopifnot(length(omegaplus) == nrow(Omega))
    res <- .Call("mmCmvplusDP", y, X, Z, H, Wcase, Wsubject, Omega, omegaplus, groups, subjects, niter, nburn, nthin, strength.mm, corsess, shapealph, ratebeta, typemm, package = "growcurves")
}
		
####################################
## accessor methods
####################################
#' S3 functions of dpgrowmm
#'
#' produces quantile summaries for model parameters
#' for an \code{dgprowmm} object.
#'
#' @param object A \code{dpgrowmm} object
#' @param ... Ignored
#' @export 
#' @method summary dpgrowmm
#' @aliases summary.dpgrowmm 
summary.dpgrowmm <- function(object,...)
{
  res 		<- list(call = object$call, summary.results = object$summary.results)
  class(res) 	<- "summary.dpgrowmm"
  return(res)
}

#' Print summary statistics for sampled model parameters
#'
#' provides credible intervals of sampled parameters for 
#' \code{dpgrowmm} object
#'
#' @param x A \code{dpgrowmm} object
#' @param ... Ignored
#' @export 
#' @method print summary.dpgrowmm
#' @noRd
print.summary.dpgrowmm <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCredible Intervals and Fit Statistics\n")
  print(x$summary.results)
}

#' Produce MCMC samples for model parameters
#'
#' provides posterior sampled values for every model parameter of a
#' \code{dpgrowmm} object
#'
#' @param object A \code{dpgrowmm} object
#' @param ... Ignored
#' @export samples dpgrowmm
#' @return res list object of class \code{samples.dpgrowmm}, \code{samples.dpgrowmult}, or \code{samples.dpgrow}
#' @export samples 
#' @method samples dpgrowmm
#' @aliases samples.dpgrowmm samples
samples 		<- 	function(object,...){
				  if(is.null(attr(object, "class"))){
    				  print("object must be of class dpgrowmm, dpgrowmult, or dpgrow")
  				}
  				else  UseMethod("samples", object)
}

samples.dpgrowmm <- function(object,...)
{
  B				<- as.data.frame(object$B)
  names(B)			<- paste(rep(1:object$Nrandom, each = object$Nsubject), rep(object$subject, object$Nrandom), sep=".") ## 1.1, 1.2, ...., 1.299
  Beta				<- as.data.frame(object$Beta)
  names(Beta)			<- colnames(object$summary.results$X)

  if( !is.null(object$U) )
  {
	Nmv				<- object$summary.results$Nmv
	if(object$summary.results$model == "mmigrp")
	{
  		res 		<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, B = B, Gamma = object$U, Eta = object$Eta,
				Residuals = object$Residuals, M = object$M, S = object$S, Num.per.cluster = object$Num, bigSmin = object$bigSmin, phat = object$phat, ordscore = object$ordscore,
				Tau.gamma = object$Tau.u, Tau.b = object$Tau.b, Tau.e = object$Tau.e)
	}else{
		res 		<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, B = B, Gamma = object$U, 
				Residuals = object$Residuals, M = object$M, S = object$S, Num.per.cluster = object$Num, bigSmin = object$bigSmin, phat = object$phat, ordscore = object$ordscore,
				Tau.gamma = object$Tau.u, Tau.b = object$Tau.b, Tau.e = object$Tau.e)
	}
	if( Nmv > 1 ) { res$Rhotau.gamma = object$Rhotau.u }
  }else{ ## model contains no session effects, U
	if(object$summary.results$model == "dp")
	{
		res 		<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, B = B, 
					Residuals = object$Residuals, M = object$M, S = object$S, Num.per.cluster = object$Num, bigSmin = object$bigSmin, phat = object$phat, ordscore = object$ordscore,
					Tau.b = object$Tau.b, Tau.e = object$Tau.e)
	}else{ ## lgm
		res 		<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, B = B, 
					Residuals = object$Residuals, Tau.b = object$Tau.b, Tau.e = object$Tau.e)
	}
  }

  if( !is.null(object$dat.growthCurve) ) ## Add growth curve data set if user chooses growth curve option
  {
	res$dat.growthCurve = object$dat.growthCurve
  }	 
  class(res) 	<- "samples.dpgrowmm"
  return(res)
}

#' Produce model plots
#'
#' Builds model plots, including MCMC trace plots, analysis of session effects and subject growth curves
#'
#' @param x A \code{dpgrowmm} object
#' @param plot.out A \code{boolean} object.  If \code{TRUE}, plots are rendered.  In either case, plots are stored
#' @param ... Ignored
#' @export 
#' @method plot dpgrowmm
#' @return res a list object of class \code{plot.dpgrowmm} of 3 items:
#'	\item{plot.results}{	\code{ggplot2} plot objects; see \code{\link{mcmcPlots}}. }
#'	\item{dat.growcurve}{	A \code{data.frame} containing fields \code{c("fit","time","subject","trt")}
#'		with \code{P*T} rows, where \code{P} is the length of \code{subject} and \code{T = 10} are the number of in-subject
#'		predictions for each subject.  This object may be used to construct additional growth curves using - see \code{\link{growplot}}.} 
#'      \item{dat.gcdata}{	A \code{data.frame} containing fields  \code{c("fit","time","subject","trt")} with \code{N} rows, where \code{N} are the
#'		number of subject-time cases.  This object contains the actual data for all subjects used to co-plot with predicted growth curves.}
#' @aliases plot.dpgrowmm 
plot.dpgrowmm <- function(x, plot.out = TRUE, ...)
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
  res		<- list(plot.results = x$plot.results, dat.growcurve = x$dat.growthCurve, dat.gcdata = x$dat.gcdata)
  class(res)	<- "plot.dpgrowmm"
  return(res)
}




  


