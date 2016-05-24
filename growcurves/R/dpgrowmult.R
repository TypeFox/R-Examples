###################################
## generic for mm-session models
###################################
#' @include growthCurve.R
#' @include XZcov.R
#' @include relabel.R
#' @include mcmcPlots.R
#' @include summary_quantiles.R
NULL

#' Bayesian semiparametric growth curve models under employment of more-than-\code{1} multiple membership random effects (block) term.
#'
#' Employs a Dirichlet Process (DP) prior on the set of by-subject random effect parameters
#' under repeated waves of measurements to allow the number of random effect parameters specified per subject, \code{q},
#' to be equal to the number of measurement waves, \code{T}.  Random effects are grouped by subject and
#' all \code{q} parameters receive the DP prior.  Additional sets of possibly more than \code{1} multiple membership effect terms
#' are included, each with a separate weight/design matrix that maps the effects back to clients.  A variety
#' of prior formulations are available for the effects in each multiple membership term.
#'
#' @param y A univariate continuous response, specified as an \emph{N x 1} matrix or vector, where \code{N}
#'	captures the number of subject-time cases (repeated subject measures).  Data may reflect unequal
#'	number of measures per subject.  Missing occasions are left out as no \code{NA} values are allowed.
#' @param subject The objects on which repeated measures are conducted that serves as the random effects
#'	grouping factor.  Input as an \emph{N x 1} matrix or vector of subject-measure cases in either
#'	integer or character formt; e.g. \code{(1,1,1,2,2,3,3,3,...,n,n,n)}, where \code{n} is the total
#'	number of subjects.
#' @param trt An integer or character vector of length \code{N} (number of cases) indicating treatment
#'	group assignments for each case.  May also be input as length \code{P} vector, where \code{P} is
#'	the number of unique subjects, indicating subject group assignment.  Multiple treatment groups
#'	are allowed and if the vector is entered as numeric, e.g. \code{(0,1,2,3,..)}, the lowest numbered
#'	group is taken as baseline (captured by global fixed effects).  If entered in character format,
#'	the first treatment entry is taken as baseline.  If the are no treatment (vs. control) groups,
#'	then this vector may be excluded (set to NULL).
#' @param time A univariate vector of length \code{N}, capturing the time points associated to each by-subject
#'	measure.  Mav leave blank if only one time point (no repeated measures).
#' @param n.random The desired number of subject random effect terms, \code{q}.  Since a DP prior is used on client effects,
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
#' @param Omega A list object of length equal to the number of multiple membership (MM) effect terms chosen with the \code{"mmcar"} option.
#'	List element \code{i} contains an \emph{S[i] x S[i]} numeric matrix to encode the CAR adjacency matrix,
#'	where \code{S[i]} is the number of effects mapped to subjects for list component \code{i}.  
#'	This input is required only under \code{option = "mmcar"}.
#' @param group	A list object of length equal to the number of MM terms chosen with prior formulation options \code{"mmcar"} or \code{"mmigrp"}.
#'	List element \code{i} contains a numeric or character vector of length \code{S[i]}, providing group identifiers for each of \code{S[i]} effects in term 
#'	\code{i}. (e.g. \code{(1,1,1,2,2,...)}.  If there is only a single group for term [i], this element should be loaded with an \code{S[i] x 1} vector of a single value.
#' @param subj.aff A list object of length equal to the number of total MM terms.  List element \code{i} contains a
#'	\code{Paff[i] x 1} vector subset of \code{subject} composed with unique subject identifiers that are linked to the effects in term \code{i};
#'	 e.g. one or more treatment cohorts.  \code{P.aff[i]} is the length of the unique subjects linked to the effects in MM term \code{i}.
#'	 If all subjects are to receive the mapping of multiple membership effects then \code{Paff[i]} should contain a list of the unique subjects in \code{subject}.
#' @param W.subj.aff A list object of length equal to the number of MM terms.  List element \code{i} contains 
#'	a \emph{P.aff[i] x S[i]} numeric matrix that maps a set of random effects to affected subjects (\code{subj.aff[[i]]}).  
#'	It is assumed that the row order is the same as the order of \code{subj.aff[[i]]}.  If \code{W.subj.aff[[i]]} is a multiple membership
#'	weight matrix, then the rows will sum to 1, though this is not required.   The rows of \code{W.subj.aff} may alternatively be formulated
#'	with indicators for whether each of \code{S} treatment dosages are linked to a given \code{subject}.
#' @param n.iter Total number of MCMC iterations.
#' @param n.burn Number of MCMC iterations to discard.  \code{dpgrow} will return \code{(n.iter - n.burn)} posterior samples.
#' @param n.thin Gap between successive sampling iterations to save.
#' @param strength.mm Sets both the shape and rate parameter for a \code{tau_{gamma} ~ G(strength.mm,strength.mm)} prior on the precision parameter of either a
#'	CAR (\emph{gamma ~ CAR(tau_gamma)}) or independent (\emph{gamma ~ N(0,tau_gamma^(-1)I_S}) prior on the set of \code{S} 
#'	multiple membership effects. Defaults to \code{strength.mm = 0.01}.
#' @param shape.dp Shape parameter under a \emph{c ~ G(shape.dp, 1)} prior on the concentration parameter of the DP (prior
#'	on the set of random effects parameters, \emph{b_1, ..., b_n ~ DP(c,G_0)}
#'	where \code{n} is the total number of subjects.
#' @param rate.dp Rate parameter under a \emph{c ~ G(shape.dp, rate.dp)} prior on the concentration parameter of the DP.
#' @param plot.out A boolean variable indicating whether user wants to return plots with output results.  Defaults to \code{TRUE}.
#' @param option A character vector of length equal to the total number of multiple membership terms that supplies the prior formulation choice
#'	for each term.  The elements of \code{option} are confined to choose from among \code{c("mmcar","mmi","mmigrp","mmdp")}.
#'	Any element of this choice set may be selected multiple times as desired.  For example, to add 3 multiple membership terms
#'	with effects in the first term under a DP prior, the second term under a CAR prior and the third also under a CAR prior, the
#'	entry would be, \code{option = c("mmdp","mmcar","mmcar")}.  The corresponding list entries should conform this choice for \code{option}.
#'	The order of sampled effect values returned conforms to this order of input in \code{option} (and the corresponding \code{subj.aff} and \code{W.subj.aff}).
#' @param ulabs A vector of the same length as \code{option} containing desired labels for each multiple membership term.  
#'	These label values are employed in returned plot objects.  If left blank, \code{ulabs} is set to a sequential number vector starting at \code{1}.
#' @return S3 \code{dpgrowmult} object, for which many methods are available to return and view results.  Generic functions applied
#'	to an object, \code{res} of class \code{dpgrow}, includes:
#'	\item{summary(res)}{returns \code{call}, the function call made to \code{dpgrowmult} and \code{summary.results}, which contains
#'					a list of objects that include \emph{95\%} credible intervals for each set of sampled parameters, 
#'					specified as (\code{2.5\%}, mean, \emph{97.5\%}, including fixed and random effects. 
#'					Also contains model fit statistics, including \code{DIC} (and associated \code{Dbar}, \code{Dhat}, \code{pD}, 
#'					\code{pV}), as well as the log pseudo marginal likelihood (LPML), a leave-one-out fit statistic.  
#'					Note that \code{DIC} is constructed as \code{DIC3} (see Celeaux et. al. 2006), where the 
#'					conditional likehihood evaluated at the posterior mode is replaced by the marginal predictive density. 
#'					Lastly, the random and fixed effects design matrices, \code{X, Z}, are returned that include 
#'					both the user input nuisance covariates appended to the time and treatment-based covariates constructed 
#'					by \code{dpgrowmult}.}  
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
#'	The \code{dpgrowmult} function generalizes \code{dpgrowmm} by allowing more than one multiple membership effects term.  
#' @keywords model
#' @seealso \code{\link{dpgrowmm}}
#' @examples 
#' \dontrun{
#' ## extract simulated dataset
#' library(growcurves)
#' data(datsimmult)
#' ## Model with DP on clients effects, but now INCLUDE session random effects
#' ## in a multiple membership construction communicated with the N x S matrix, W.subj.aff.
#' ## Returns object, res.mm, of class "dpgrowmm".
#' shape.dp	= 3
#' res.mult	= dpgrowmult(y = datsimmult$y, subject = datsimmult$subject, 
#'			trt = datsimmult$trt, time = datsimmult$time, 
#'			n.random = datsimmult$n.random, Omega = datsimmult$Omega, 
#'			group = datsimmult$group, 
#'			subj.aff = datsimmult$subj.aff, 
#'			W.subj.aff = datsimmult$W.subj.aff, n.iter = 10000, 
#'			n.burn = 2000, n.thin = 10, shape.dp = shape.dp, 
#'			option = c("mmi","mmcar"))
#' plot.results	= plot(res.mult) ## ggplot2 plot objects, including growth curves
#' summary.results = summary(res.mult) ## parameter credible intervals, fit statistics
#' samples.posterior = samples(res.mult) ## posterior sampled values
#' }
#' @aliases dpgrowmult
#' @aliases dpgrowmult.default
#' @author Terrance Savitsky \email{tds151@@gmail.com} Susan Paddock \email{paddock@@rand.org}
#' @references
#'	S. M. Paddock and T. D. Savitsky (2012) Bayesian Hierarchical Semiparametric Modeling of Longitudinal Post-treatment Outcomes from Open-enrollment Therapy Groups, invited re-submission to: JRSS Series A (Statistics in Society).
#' @references
#' 	T. D. Savitsky and S. M. Paddock (2012) Visual Sufficient Statistics for Repeated Measures data with growcurves for R, submitted to: Journal of Statistical Software.
#' @export dpgrowmult 
dpgrowmult			<- function(y, subject, trt, time, n.random, n.fix_degree, formula, random.only, data, Omega, group, subj.aff, W.subj.aff, n.iter, 
					n.burn, n.thin, strength.mm, shape.dp, rate.dp, plot.out, option, ulabs)
					UseMethod("dpgrowmult")

################################################
## default dispatch method for mm-session models
################################################
dpgrowmult.default		<- function(y = NULL, subject, trt = NULL, time = NULL, n.random = NULL, n.fix_degree = NULL, formula = NULL, 
						random.only = FALSE, data = NULL, Omega = NULL, group = NULL, subj.aff, W.subj.aff, 
					 	n.iter, n.burn, n.thin = 1, strength.mm = 0.35, shape.dp = 1, rate.dp = 1, plot.out = TRUE, option, ulabs = NULL)
{ ## start function dpgrowmult.default

  ############################
  ## check inputs
  ############################
  ## model choices
  if( is.null(option) ) stop("Must enter character vector of prior choices of length equal to number of multiple membership terms.")

  option 	= tolower(option)
  choices 	= c("mmcar","mmi","mmigrp","mmdp")
  if( length(setdiff(option,choices)) > 1)
  {
	stop("Your option choices must be confined to within, c('mmcar','mmi','mmigrp','mmdp')")
  }

  if( (is.null(Omega) | is.null(group)) & ("mmcar" %in% option) ) ## inconsistencies
  {
	stop("You must specify list objects Omega (containing (S = number of effects) X S adjacency matrices) and 'group' (S X 1 vectors of group identifiers) for each 'mmcar' term in 'option'.")
  }

  if ( is.null(group) & ("mmigrp" %in% option | "mmcar" %in% option) )
  {
	stop("You must input a list object 'group' with a vector entry for each term of length number of effects for that term for each 'mmcar' and 'mmigrp' element in 'option'. Set values to 1 if only a single group.")
  }

  ## data choices 
  nty = length(option)


  ## W.subj.aff is a list of P.aff[i] x S[i] mm weight matrices.  subj.aff is a list of P.aff[i] x 1 vectors with entries drawn from 'subject'

  ## checking to see that each subj.aff vector drawn from the subject vector
  if( is.null(subj.aff) )
  {
	stop("Must input list object of vectors, subj.aff, that identifies subset of subjects affected by each multiple membership term in 'option'.")
  }else{
	if( length(subj.aff) != length(option) ) stop("There must be a vector entry in 'subj.aff' for each corresponding element in 'option'.")
  	for( i in 1:nty ) ## make sure that subj.aff is a strict subset of subject
  	{
		if( length(setdiff(subj.aff[[i]],subject)) != 0 ) stop("Each vector in list object, 'subj.aff', must only contain elements from 'subject'.")
  	}
  }

  ## checking to see that each element of W.subj.aff entered in subject, not case, format.
  if( is.null(W.subj.aff) ) 
  {
  	stop("Must input list object of multiple membership design matrices, 'W.subj.aff', of length number of elements in 'option'. ")
  }else{ ## !is.null(W.subj.aff)
        if( length(W.subj.aff) != length(option) ) stop("There must be a matrix entry in 'W.subj.aff' for each corresponding element in 'option'.")
	for(i in 1:nty)
	{
		if( nrow(W.subj.aff[[i]]) == length(subject) )  ## input in case format, rather than subject format
		{
			warning("Rows of W.subj.aff should be equal to length of unique affected subjects, not number of subject-time cases.  
				Function will use unique ID's in subject vector to shrink rows of W.subj.aff.")
			tmp		= as.data.frame(cbind(subject,W.subj.aff[[i]])) ## note: only want shrinkage happening within subject id
			tmp		= unique(tmp)
			W.subj.aff[[i]]	= tmp[,-1]
			rm(tmp)
		}
                if( length(subj.aff[[i]]) > length(unique(subj.aff[[i]])) ) ## ensure of length Nsubject, not Ncase
		{
			warning("Each vector 'subj.aff' should contain number of unique subjects, not subject-time cases.  Function will shrink to unique values.")
			subj.aff[[i]] 		= unique(subj.aff[[i]])  
		}
		if( nrow(W.subj.aff[[i]]) != length(subj.aff[[i]]) )   ## subj.aff now in subject format.  Ensuring correspondence b/t W.subj.aff and subj.aff
		{
			stop("Number of entries in each vector element of 'subj.aff' must equal number of rows in corresponding matrix element of 'W.subj.aff'.")
		}
	}
  }

  ## labels
  if( is.null(ulabs) )
  {
	ulabs = 1:length(option)
  }else{
  	if( length(ulabs) != length(option) ) stop("Length of 'ulabs' name vector must be same as 'option' vector.")
  }

  ## Omega is a list of S[i] x S[i] adjacency matrices associated to option == "mmcar"
  ## Create list of S[i] x 1 vectors, omega.plus, from rowSums of Omega
  ## group is a list of S[i] x 1 group identifiers for effects under options == c("mmcar","mmigrp")
  ## for option == "mmcar", create vector 'ngs', that sums number of unique groups in applicable group[[i]]
  if( "mmcar" %in% option | "mmigrp" %in% option )
  {
	poscar 		= which(option == "mmcar")
	numcar 		= length(poscar)
	if( numcar > 0 ) ## allow for possibility "mmcar" not selected as an option
	{
		omega.plus	= vector("list",numcar) ## row sums of Omega[[i]]
		for(i in 1:numcar)
		{
			if( ncol(W.subj.aff[[poscar[i]]]) != ncol(Omega[[i]]) ) stop("\nNumber of columns of respective matrices in 'W.subj.aff' and 'Omega' must be equal (to the number of effects) under the 'mmcar' option choice.\n")
			if( ncol(Omega[[i]]) != nrow(Omega[[i]]) ) stop("Each adjacency matrix in Omega must be a square matrix dimensioned by number of multi-mbrship random effects under the 'mmcar' option choice.")
			omega.plus[[i]]		= rowSums(Omega[[i]])
		
		}
	}

	posicar = which(option == "mmcar" | option == "mmigrp")
	numicar = length(posicar)
	if( length(group) != numicar ) stop("Number of elements in list object, 'group', must be equal to number of option choices equal to 'mmcar' and 'mmigrp'.")

        ngs	= vector("numeric", numcar) ## vector of number of groups for CAR
	countcar = 1; 
	for( i in 1:numicar )
	{
		if( ncol(W.subj.aff[[posicar[i]]]) != length(group[[i]]) ) stop("The length of each vector in 'group' must equal the number of columns of each matrix in 'W.subj.aff' - the number of effects.")
		if(option[posicar[i]] == "mmcar")
		{
			if( ncol(Omega[[countcar]]) != length(group[[i]]) ) stop("The length of each vector in 'group' must equal the row and column dimensions of the corresponding adjacency matrix, Omega, under option choice 'mmcar'.")
			ngs[countcar] 		= length(unique(group[[i]]))
			countcar 		= countcar + 1
		}
	}
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
        	if( length(trt) == length(unique(subject)) ) ## input in subject, rather than case format
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

  ## subj.aff - capture as strict subset of subject
  ## don't re-order subj.aff[[i]] b/c correponds to rows in W.subj.aff[[i]].  only change labels
  subjaff.input 			<- vector("list",nty)
  smap.u				<- unique(map.subject) ## in subject format
  for( i in 1:nty )
  {
  	subjaff.input[[i]] 		<- subj.aff[[i]]
	o.aff.i				<- 1:length(subj.aff[[i]]) ## will use to maintain order of input
  	tmp				<- data.frame(o.aff.i,subjaff.input[[i]])
  	names(tmp)			<- c("order","label.input")
  	smap.u.i			<- subset(smap.u, smap.u$label.input %in% subjaff.input[[i]])
  	tmp				<- merge(tmp,smap.u.i, by = "label.input", all.x = T, sort =  FALSE)
  	subj.aff[[i]]			<- tmp$label.new[order(tmp$order)] ## result is a strict subset of subject, but in subject, not case format
  }
  rm(tmp,smap.u)
  
  ## trt
  start		<- 0
  out		<- relabel(label.input = trt, start)
  trt		<- out$label.new
  trti.u	<- out$labeli.u
  map.trt	<- out$dat.label

  ## group - need a list of length nty for mcmcPlots function.  list elements set equal to 1's, except for those associated to "mmcar" or "mmigrp"
  group.i 	<- vector("list",nty) ## create a group vector for mcmcPlots, but set equal to 1's for options !%in% c("mmcar","mmigrp")
  numt		<- vector("numeric", nty)
  for( i in 1:nty)
  {
	numt[i]		= ncol(W.subj.aff[[i]])
	group.i[[i]] 	= matrix(1,numt[i],1) ## just set group labels all equal to 1
  }

  if( "mmcar" %in% option | "mmigrp" %in% option )
  {
	for( i in 1:numicar )
	{
  		start			<- 1
  		outi			<- relabel(label.input = group[[i]], start)
  		group[[i]]		<- outi$label.new  ## list dimensioned to numicar
		group.i[[posicar[i]]]	<- as.matrix(outi$dat.label$label.input) ## to capture those members of group.i associated to "mmcar" or "mmi" with length equal to number of sessions
	}
   } ## end conditional statement on whether "mmcar" or "mmigrp" %in% option

  ## for option == "mmigrp", create S[i] x ng[i] matrix verson of group[[i]] for each time this option is chosen and wrap it into a list element, Mmats.
  if("mmigrp" %in% option)
  {
	posigrp 	= which(option == "mmigrp")
	numigrp 	= length(posigrp)
	Mmats 		= vector("list", numigrp)
	countigrp	= 1
	for( i in 1:numicar ) ## cycle through both "mmcar" and "mmigrp" option choices since 'group' input for both.
        {
		dimi  			= length(group[[i]]) ## number of sessions
		ngi   			= length(unique(group[[i]])) ## number of groups
		if(option[i] == "mmigrp")
  		{
			mmat 			= matrix(0,dimi,ngi)
			for(j in 1:ngi)
			{
				mmat[group[[i]] == j, j]	= 1
			}
			Mmats[[countigrp]] 	= mmat
			countigrp 		= countigrp + 1
		} ## end conditional statment on whether current option is "mmigrp"
	} ## end loop over options = c("mmcar","mmigrp")
  } ## end conditional statment of whether "mmigrp" %in% option
  

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
  X		<- out$X
  X.c   	<- out$X.c
  X.n   	<- out$X.n
  Z		<- out$Z
  Z.n		<- out$Z.n
  Z.c  		<- out$Z.c
  if( !is.null(out$y) ) 
  {
	y 	<- out$y  ## over-writes possible duplicative input of y by user (since must be in formula).
  }else{ ## out$y is null, so user separately entered
	y	<- y[o] ## re-order y by subject to ensure subject is in contiguous order
  }
  y	= as.matrix(y,Ncase,1)

  ## reorder remaining objects to subject (in contiguous fashion) where entries indexed by case
  subject		<- subject[o]
  map.subject		<- map.subject[o,]
  map.trt		<- map.trt[o,]
  time			<- time[o] ## used for growth curve plotting

  ## capture number of fixed effects
  Nfixed		= ncol(X)
  Nrandom		= ncol(Z)

  #################################################################
  ## re-cast inputs to matrices
  #################################################################
  ## Expand W.subj.aff of rows subj.aff to Ncase rows, W.case
  W.subj 	= vector("list",nty)
  W.case	= vector("list",nty)
  for(i in 1:nty)
  {
	W.subj.aff[[i]]				= as.matrix(W.subj.aff[[i]])
  	W.subj[[i]]				= matrix(0, nrow = Nsubject, ncol = numt[i]) ## non-affected subjects have rows of 0. Remember, 'subject' starts at '1'.
  	W.subj[[i]][subj.aff[[i]],]		= W.subj.aff[[i]]
	W.subj[[i]]				= as.matrix(W.subj[[i]])
	W.case[[i]]				= W.subj[[i]][subject,] ## inflate to rows equal to number of subject-time cases
	W.case[[i]]				= as.matrix(W.case[[i]])
  } ## end loop on creating P x S[i] (W.subj) and N x S[i] (W.case) mm weight matrices

  ################################################################
  ## conduct posterior sampling and capture results
  ################################################################
  print(paste("Your chosen set of MM term priors  = ", paste(option,collapse=" "), sep = ""))

  ## construct numeric version of option vector for input to C++ engine
  typet 	= vector("numeric",nty)
  for( i in 1:nty )
  {
	## build numeric typet C++ input vector from user-entered 'option'.
	if(option[i] == "mmcar") typet[i] = 1
	else if(option[i] == "mmi") typet[i] = 2
	else if(option[i] == "mmigrp") typet[i] = 3
	else typet[i] = 4 ## "mmdp"
  }

  ## build dummy list elements for C++ function in case of no selection of some options choices
  if( !("mmcar" %in% option) )
  {
	Omegamat 	= matrix(0, 2, 2)
	omegaplus	= matrix(0, 1, 2)
	Omega		= list(Omega = Omegamat)
	omega.plus	= list(omega.plus = omegaplus)
	ngs		= 1
  }

  if( !("mmigrp" %in% option) )
  {
	mmat		= diag(rep(1,2))
	Mmats		= list(Mmat = mmat)
  }

  shape.dp	= 4
  res 		= mmmultPost(y, X, Z, W.case, Mmats, Omega, omega.plus, ngs, subject, typet, n.iter, n.burn, n.thin, shape.dp, rate.dp, strength.mm)
  stopifnot(length(res) == 19)
  stopifnot(ncol(res$Tauu) == nty)
  stopifnot(is.list(res$U))
  Utest = do.call("cbind",res$U)
  stopifnot(ncol(Utest) == sum(numt))

  ##################################################################
  ## summary (short-hand) results
  ##################################################################
  summary.results			<- summary_quantiles(model.output = res, Nfixed = Nfixed, Nrandom = Nrandom, Nsubject = Nsubject)
  summary.results$X			<- X
  summary.results$Z			<- Z
  summary.results$n.fix_degree		<- n.fix_degree
  summary.results$map.subject		<- map.subject
  summary.results$time			<- time ## not used in accessor functions; just reporting back to user to let them know that sorted by subject
  summary.results$map.trt		<- map.trt
  summary.results$model 		<- option
  summary.results$group.i		<- group.i
  stopifnot(length(summary.results$u.summary) == length(ulabs))
  names(summary.results$u.summary)	<- as.character(ulabs) ## list object
  summary.results$numt			<- numt
  summary.results$ulabs			<- ulabs ## for naming elements of U list object

  residuals					= colMeans(res$Residuals)
 
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
			## memo: if have nuisance covariates, need input of Nrandom to construct time-based random effects since Q > Nrandom
   }

 } ## end conditional statement on creating growth curves

   
 if(plot.out == TRUE)
 {
   ##################################################################
   ## plots
   ##################################################################
   ## note: group.i is of length numt[i] (number of sessions)
   plot.results = mcmcPlots(subjecti.u = subjecti.u, bmat.summary = summary.results$bmat.summary,  
				groupi.u = group.i, u.summary = summary.results$u.summary, 
				M = res$M, Tauu = res$Tauu, Taub = res$Taub, Taue = res$Taue, Deviance = res$Deviance, ulabs = ulabs)

 } #end conditional statement on whether to plot 
   
 ##################################################################
 ## function output
 ##################################################################
 if(plot.out == TRUE )
 {
   if( (!is.null(time) & length(unique(time)) > 1) & !is.null(n.fix_degree) ) ## a set of growth curves were generated from time-based covariates
   {
   	plot.results$p.gcall = gc.plot$p.gcall; plot.results$p.gcsel = gc.plot$p.gcsel
   	if( "mmdp" %in% option )
   	{
		
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U,  M = res$M, S = res$optPartition[[3]], 
			Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, bigSu = res$bigSu, Mu = res$Mu, Su = res$Su, summary.results = summary.results, plot.results = plot.results,
			residuals = residuals, dat.growthCurve = gc.plot$plot.dat, dat.gcdata = gc.plot$dat.data) 
   	}else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, 
			M = res$M, S = res$S, Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, plot.results = plot.results,
			residuals = residuals, dat.growthCurve = gc.plot$plot.dat, dat.gcdata = gc.plot$dat.data)
   	}
   }else{ ## is.null(time) = TRUE
     	if( "mmdp" %in% option )
     	{
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U,  M = res$M, S = res$optPartition[[3]], 
			Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, bigSu = res$bigSu, Mu = res$Mu, Su = res$Su, summary.results = summary.results, plot.results = plot.results, residuals = residuals)
     	}else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, 
			M = res$M, S = res$optPartition[[3]], Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, plot.results = plot.results, residuals = residuals)
   	}  ## end conditional statement on choice	
   } ## end conditional statement on whether is.null(time)
 }else{ ## plot.out = FALSE
   if(  "mmdp" %in% option )
   {
    		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U,  M = res$M, S = res$optPartition[[3]], 
			Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, bigSu = res$bigSu, Mu = res$Mu, Su = res$Su, summary.results = summary.results, residuals = residuals)
   }else{
		resot = list(Deviance = res$Deviance, Beta = res$Beta, Alpha = res$Alpha, B = res$B, U = res$U, 
			M = res$M, S = res$optPartition[[3]], Num = res$Num, Residuals = res$Residuals, bigSmin = res$bigSmin, phat = res$optPartition[[1]], ordscore = res$optPartition[[2]],
			Tau.u = res$Tauu, Tau.e = res$Taue, Tau.b = res$Taub, summary.results = summary.results, residuals = residuals)
   }
 } ## end conditional statement on plot.out

 ##
 ## return list output for dpgrowmm.default()
 ##

 resot$call		<- match.call()
 resot$Nrandom     	<- ncol(resot$summary.results$Z)
 resot$Nsubject		<- length(unique(subject))
 resot$subject		<- unique(subjecti.u) ## will employ for labeling B with user input subject labels
 class(resot)		<- c("dpgrowmult")
 return(resot)

} ### end function dpgrowmm.default()


#####################################################
## .Call statements to C++ functions
#####################################################
#' Bayesian mixed effects model with a DP prior on by-subject effects and more than one multiple membership random effects term
#'
#' An internal function to \code{\link{dpgrowmult}}
#'
#' @export mmmultPost
#' @aliases mmmultPost 
#' @param y An \emph{N x 1} response (of subject-measure cases)
#' @param X Fixed effects design matrix
#' @param Z Random effects design matrix.  Assumed grouped by \code{subjects}
#' @param Wcases A list object containing \code{N x 1} multiple membership (MM) weight matrices; one for each MM term.
#' @param Mmats A list object containing \code{S[i] x G[i]} matrices, where \code{S[i]} are the number of effects for MM term \code{i} and 
#'		\code{G[i]} is the number of (unique) groups for the \code{S[i]} effects.  The length of this list is equal to the number of terms under prior option set to \code{"mmigrp"}
#' @param Omegas A list object containing \code{S[i] x S[i]} unnormalized adjacency matrices for MM term \code{i}, each with entries equal to 1 where two effects communicate
#'	and 0, otherwise.  Diagonal elements are zero.  The length of this list is equal to the number terms with prior option set to \code{"mmcar"}
#' @param omegapluses A list object containing \code{S x 1} vector of row sums of \code{Omega[[i]]}.  The length of \code{omegaplus} should equal \code{Omega}.
#' @param ngs A numeric vector containing the number of total groups for each MM term under prior option \code{"mmcar"}.
#' @param subjects An \emph{N x 1} set of subject identifiers.
#' @param typet A numeric vector of length equal to the number of MM terms, where each entry specifies the prior formulation for that effect term (block).
#'		Prior formulation options are \code{1 = "mmcar", 2 = "mmi", 3 = "mmigrp", 4 = "mmdp"}.
#' @param niter The number of MCMC iterations
#' @param nburn The number of MCMC burn-in iterations to discard
#' @param nthin The step increment of MCMC samples to return
#' @param shapealph The shape parameter for the \eqn{\Gamma} prior on the DP concentration parameter.  
#'	The rate parameter is set of \code{1}.
#' @param ratebeta The rate parameter for the \eqn{\Gamma} prior on the DP concentration parameter. Default value is \code{1}.
#' @param strength.mm The shape and rate parameters for the \eqn{\Gamma} prior on the CAR precision parameter, \eqn{\tau_{\gamma}}
#' @return res A list object containing MCMC runs for all model parameters.
#' @seealso \code{\link{dpgrowmm}}
#' @author Terrance Savitsky \email{tds151@@gmail.com}
#' @note Intended as an internal function for \code{\link{dpgrowmult}}
mmmultPost = function (y, X, Z, Wcases, Mmats, Omegas, omegapluses, ngs, subjects, typet, niter, nburn, nthin, shapealph, ratebeta, strength.mm) {
    stopifnot(nrow(X) == nrow(Z))
    stopifnot(nrow(Wcases[[1]]) == nrow(Z))
    stopifnot(length(y) == nrow(X))
    stopifnot(length(omegapluses) == length(Omegas))
    stopifnot(length(Wcases) == length(typet))
    res <- .Call("mmmult", y, X, Z, Wcases, Mmats, Omegas, omegapluses, ngs, subjects, typet, niter, nburn, nthin, shapealph, ratebeta, strength.mm, package = "growcurves")
}

####################################
## accessor methods
####################################
#' S3 functions of dpgrowmult
#'
#' produces quantile summaries for model parameters
#' for an \code{dgprowmult} object.
#'
#' @param object A \code{dpgrowmult} object
#' @param ... Ignored
#' @export 
#' @method summary dpgrowmult
#' @aliases summary.dpgrowmult 
summary.dpgrowmult <- function(object,...)
{
  res 		<- list(call = object$call, summary.results = object$summary.results)
  class(res) 	<- "summary.dpgrowmult"
  return(res)
}

#' Print summary statistics for sampled model parameters
#'
#' provides credible intervals of sampled parameters for 
#' \code{dpgrowmult} object
#'
#' @param x A \code{dpgrowmult} object
#' @param ... Ignored
#' @export 
#' @method print summary.dpgrowmult
#' @noRd
print.summary.dpgrowmult <- function(x,...)
{
  cat("Call:\n")
  print(x$call)
  cat("\nCredible Intervals and Fit Statistics\n")
  print(x$summary.results)
}

#' Produce samples of MCMC output
#'
#' provides posterior sampled values for every model parameter of a
#' \code{dpgrowmult} object
#'
#' @param object A \code{dpgrowmult} object
#' @param ... Ignored
#' @export samples dpgrowmult
#' @aliases samples.dpgrowmult
#' @method samples dpgrowmult
#' @aliases samples.dpgrowmult
samples.dpgrowmult <- function(object,...)
{
  B				<- as.data.frame(object$B)
  names(B)			<- paste(rep(1:object$Nrandom, each = object$Nsubject), rep(object$subject, object$Nrandom), sep=".") ## 1.1, 1.2, ...., 1.299
  Beta				<- as.data.frame(object$Beta)
  names(Beta)			<- colnames(object$summary.results$X)

  names(object$U)		<- as.character(object$summary.results$ulabs)

  if("mmdp" %in% object$summary.results$model )
  {
  	res 		<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, B = B, Gamma = object$U, 
			   	Residuals = object$Residuals, M = object$M, S = object$S, Num.per.cluster = object$Num, bigSmin = object$bigSmin, phat = object$phat, ordscore = object$ordscore,
			   	Tau.gamma = object$Tau.u, Tau.b = object$Tau.b, Tau.e = object$Tau.e, bigSgamma = object$bigSu, Mu = object$Mgamma, Sgamma = object$Su)
  }else{
	res 		<- list(Deviance = object$Deviance, Alpha = object$Alpha, Beta = Beta, B = B, Gamma = object$U, 
				Residuals = object$Residuals, M = object$M, S = object$S, Num.per.cluster = object$Num, bigSmin = object$bigSmin, phat = object$phat, ordscore = object$ordscore,
				Tau.gamma = object$Tau.u, Tau.b = object$Tau.b, Tau.e = object$Tau.e)
  }

  if( !is.null(object$dat.growthCurve) ) ## Add growth curve data set if user chooses growth curve option
  {
	res$dat.growthCurve = object$dat.growthCurve
  }	 
  class(res) 	<- "samples.dpgrowmult"
  return(res)
}

#' Produce model plots
#'
#' Builds model plots, including MCMC trace plots, analysis of session effects and subject growth curves
#'
#' @param x A \code{dpgrowmult} object
#' @param plot.out A \code{boolean} object.  If \code{TRUE}, plots are rendered.  In either case, plots are stored
#' @param ... Ignored
#' @export 
#' @method plot dpgrowmult
#' @return res a list object of class \code{plot.dpgrowmm} of 3 items:
#'	\item{plot.results}{	\code{ggplot2} plot objects; see \code{\link{mcmcPlots}}. }
#'	\item{dat.growcurve}{	A \code{data.frame} containing fields \code{c("fit","time","subject","trt")}
#'		with \code{P*T} rows, where \code{P} is the length of \code{subject} and \code{T = 10} are the number of in-subject
#'		predictions for each subject.  This object may be used to construct additional growth curves using - see \code{\link{growplot}}.} 
#'      \item{dat.gcdata}{	A \code{data.frame} containing fields  \code{c("fit","time","subject","trt")} with \code{N} rows, where \code{N} are the
#'		number of subject-time cases.  This object contains the actual data for all subjects used to co-plot with predicted growth curves.}
#' @aliases plot.dpgrowmult 
plot.dpgrowmult <- function(x, plot.out = TRUE, ...)
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
  class(res)	<- "plot.dpgrowmult"
  return(res)
}