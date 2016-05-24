## Update on Sep 19 2015 for R2
## Add phy.code and phy.cat
## 		- phy.code: N x C.phy matrix, first column ID
## 		- phy.cat: C x 2 matrix, translation rule
##      - phy.unknown: scalar, value of unknown category in phy.code
##      - things changed: pnb
## Update on Sep 30 2015 
## Add probbase overwrite, and new SCI compatibility


#' Implement InSilicoVA methods
#' 
#' This function implements InSilicoVA model. The InSilicoVA model is fitted
#' with MCMC implemented in Java. For more detail, see the paper on
#' \url{http://arxiv.org/abs/1411.3042}.
#' 
#' For Windows user, this function will produce a popup window showing the
#' progress. For Mac and Unix user, this function will print progress messages
#' on the console. Special notice for users using default R GUI for mac, the
#' output will not be printed on console while the function is running, and
#' will only be printed out after it is completed. Thus if you use a Mac, we
#' suggest using either RStudio for mac, or running R from terminal.
#' 
#' The chains could be set to run automatically longer. If set
#' \code{auto.length} to be TRUE, the chain will assess convergence after
#' finishing the length K chain input by user using Heidelberger and Welch's
#' convergence diagnostic. If convergence is not reached, the chain will run
#' another K iterations and use the first K iterations as burn-in. If the chain
#' is still not converged after 2K iterations, it will proceed to another 2K
#' iterations and again use the first 2K iterations as burn-in. If convergence
#' is still not reached by the end, it will not double the length again to
#' avoid heavy memory use. A warning will be given in that case. The extended
#' chains will be thinned in the same way.
#' 
#' For more detail of model specification, see the paper on
#' \url{http://arxiv.org/abs/1411.3042}.
#' 
#' @aliases insilico print.insilico
#' @param data The original data to be used. It is suggested to use similar
#' input as InterVA4, with the first column being death IDs and 245 symptoms. 
#' The only difference in input is InsilicoVA takes three levels: ``present'',
#' ``absent'', and ``missing (no data)''. Similar to InterVA software,
#' ``present'' symptoms takes value ``Y''; ``absent'' symptoms take take value
#' ``NA'' or ``''. For missing symptoms, e.g., questions not asked or answered
#' in the original interview, corrupted data, etc., the input should be coded
#' by ``.'' to distinguish from ``absent'' category. The order of the columns does
#' not matter as long as the column names are correct. It can also include more 
#' unused columns than the standard InterVA4 input. But the first column should be 
#' the death ID. For example input data format, see \code{RandomVA1} and 
#' \code{RandomVA2}.
#' @param isNumeric Indicator if the input is already in numeric form. If the
#' input is coded numerically such that 1 for ``present'', 0 for ``absent'',
#' and -1 for ``missing'', this indicator could be set to True to avoid
#' conversion to standard InterVA format.
#' @param updateCondProb Logical indicator. If FALSE, then fit InSilicoVA model without re-estimating conditional probabilities.
#' @param keepProbbase.level Logical indicator when \code{updateCondProb} is
#' FALSE. If TRUE, then only estimate the InterVA's conditional probability
#' interpretation table; if FALSE, estimate the whole conditional
#' probability matrix. Default to TRUE.
#' @param CondProb Customized conditional probability matrix to use.It should be strict the same configuration as InterVA-4 software. That is, it should be a matrix of 245 rows of symptoms and 60 columns of causes, arranged in the same order as in InterVA-4 specification. The elements in the matrix should be the conditional probability of corresponding symptom given the corresponding cause, represented in alphabetic form indicating levels. For example input, see \code{\link{condprob}}
#' @param CondProbNum Customized conditional probability matrix to use if specified fully by numerical values between 0 and 1. If it is specified, re-estimation of conditional probabilities will not be performed, i.e., \code{updateCondProb} will be set to FALSE.
#' @param datacheck Logical indicator for whether to check the data satisfying
#' InterVA rules. Default set to be TRUE. If \code{warning.write} is set to
#' true, the inconsistent input will be logged in file warnings.txt. It's
#' strongly suggested to be set to TRUE.
#' @param datacheck.missing Logical indicator for whether to perform data check before deleting complete missing symptoms. Default to TRUE.
#' @param warning.write Logical indicator for whether to save the changes made
#' to data input by \code{datacheck}. If set to TRUE, the changes will be
#' logged in file warnings.txt in current working directory.
#' @param external.sep Logical indicator for whether to separate out external
#' causes first. Default set to be TRUE. If set to TRUE, the algorithm will
#' estimate external causes, e.g., traffic accident, accidental fall, suicide,
#' etc., by checking the corresponding indicator only without considering other
#' medical symptoms. It is strongly suggested to set to be TRUE.
#' @param Nsim Number of iterations to run. Default to be 4000.
#' @param thin Proportion of thinning for storing parameters. For example, if
#' thin = k, the output parameters will only be saved every k iterations.
#' Default to be 10
#' @param burnin Number of iterations as burn-in period. Parameters sampled in
#' burn-in period will not be saved.
#' @param auto.length Logical indicator of whether to automatically increase
#' chain length if convergence not reached.
#' @param conv.csmf Minimum CSMF value to check for convergence if auto.length
#' is set to TRUE. For example, under the default value 0.02, all causes with
#' mean CSMF at least 0.02 will be checked for convergence.
#' @param jump.scale The scale of Metropolis proposal in the Normal model.
#' Default to be 0.1.
#' @param levels.prior Vector of prior expectation of conditional probability
#' levels. They do not have to be scaled. The algorithm internally calibrate
#' the scale to the working scale through \code{levels.strength}. If NULL the
#' algorithm will use InterVA table as prior.
#' @param levels.strength Scaling factor for the strength of prior beliefs in
#' the conditional probability levels. Larger value constrain the posterior
#' estimates to be closer to prior expectation. Defult value 1 scales
#' \code{levels.prior} to a suggested scale that works empirically.
#' @param trunc.min Minimum possible value for estimated conditional
#' probability table. Default to be 0.0001
#' @param trunc.max Maximum possible value for estimated conditional
#' probability table. Default to be 0.9999
#' @param subpop This could be the column name of the variable in data that is to
#' be used as sub-population indicator, or a list of column names if more than one 
#' variable are to be used. Or it could be a vector of sub-population assignments 
#' of the same length of death records. It could be numerical indicators or character 
#' vectors of names. 
#' @param java_option Option to initialize java JVM. Default to ``-Xmx1g'',
#' which sets the maximum heap size to be 1GB. If R produces
#' ``java.lang.OutOfMemoryError: Java heap space'' error message, consider
#' increasing heap size using this option, or one of the following: (1)
#' decreasing \code{Nsim}, (2) increasing \code{thin}, or (3) disabling
#' \code{auto.length}.
#' @param seed Seed used for initializing sampler. The algorithm will produce
#' the same outcome with the same seed in each machine.
#' @param phy.code A matrix of physician assigned cause distribution. The
#' physician assigned causes need not be the same as the list of causes used in
#' InSilicoVA and InterVA-4. The cause list used could be a higher level
#' aggregation of the InSilicoVA causes. See \code{phy.cat} for more detail.
#' The first column of \code{phy.code} should be death ID that could be matched
#' to the symptom dataset, the following columns are the probabilities of each
#' cause category used by physicians.
#' @param phy.cat A two column matrix describing the correspondence between
#' InSilicoVA causes and the physician assigned causes. Note each InSilicoVA
#' cause (see \code{causetext}) could only correspond to one physician assigned
#' cause. See \code{SampleCategory} for an example. 'Unknown' category should
#' not be included in this matrix.
#' @param phy.unknown The name of the physician assigned cause that correspond
#' to unknown COD.
#' @param phy.external The name of the physician assigned cause that correspond
#' to external causes. This will only be used if \code{external.sep} is set to
#' TRUE. In that case, all external causes should be grouped together, as they
#' are assigned deterministically by the corresponding symptoms.
#' @param phy.debias Fitted object from physician coding debias function (see
#' \code{\link{physician_debias}}) that overwrites \code{phy.code}.
#' @param exclude.impossible.cause logical indicator to exclude impossible causes based on the age and gender of the death.
#' @param indiv.CI credible interval for individual probabilities. If set to NULL, individual COD distributions will not be calculated to accelerate model fitting time. See \code{\link{get.indiv}} for details of updating the C.I. later after fitting the model.
#' @param ... not used
#' 
#' @return \item{id}{A vector of death ID. Note the order of the ID is in
#' general different from the input file. See \code{report} for organizing the
#' report.}
#' \item{data}{Cleaned numerical data.}
#' \item{indiv.prob}{Matrix of individual mean cause of death distribution.
#' Each row corresponds to one death with the corresponding ID.}
#' 
#' \item{csmf}{Matrix of CSMF vector at each iterations after burn-in and
#' thinning. Each column corresponds to one cause.}
#' 
#' \item{conditional.probs}{If the model is estimated with
#' \code{keepProbbase.level} = TRUE, this value gives a matrix of each
#' conditional probability at each level at each iterations. Each column
#' corresponds to one level of probability. If \code{keepProbbase.level} =
#' FALSE, this value gives a three-dimensional array. If \code{updateCondProb} =
#' FALSE, the value will be set to NULL. See \code{report} for more analysis.}
#' 
#' \item{missing.symptoms}{Vector of symptoms missing from all input data.}
#' 
#' \item{external}{Logical indicator of whether the model is fitted with
#' external causes separated calculated.}
#' \item{impossible.causes}{Impossible cause-symptom pairs, if any.}
#' \item{indiv.CI}{The posterior credible interval to compute for individual COD probability distributions. If set to NULL, only the posterior mean of the individual COD probabilities will be produced. Default to be 0.95.}
#' \item{indiv.prob.median}{median probability of each cause of death for each individual death.}
#' \item{indiv.prob.lower}{lower CI bound for the probability of each cause of death for each individual death.}
#' \item{indiv.prob.upper}{upper CI bound for the probability of each cause of death for each individual death.}
#'  
#' @author Zehang Li, Tyler McCormick, Sam Clark
#' 
#' Maintainer: Zehang Li <lizehang@@uw.edu>
#' @seealso \code{\link{plot.insilico}}, \code{\link{summary.insilico}}, \code{\link{physician_debias}}
#' @references Tyler H. McCormick, Zehang R. Li, Clara Calvert, Amelia C.
#' Crampin, Kathleen Kahn and Samuel J. Clark(2014) \emph{Probabilistic
#' cause-of-death assignment using verbal autopsies},
#' \url{http://arxiv.org/abs/1411.3042} \cr \emph{Working paper no. 147, Center
#' for Statistics and the Social Sciences, University of Washington}
#' @keywords InSilicoVA
#' @examples
#' \dontrun{
#' data(RandomVA1) 
#' fit0<- insilico(RandomVA1, subpop = NULL,  
#'                 Nsim = 20, burnin = 10, thin = 1 , seed = 1,
#' 			 auto.length = FALSE)
#' summary(fit0)
#' summary(fit0, id = "d199")
#' 
#' ##
#' ## Scenario 1: standard input without sub-population specification
#' ##
#' fit1<- insilico(RandomVA1, subpop = NULL,  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#' 		   auto.length = FALSE)
#' summary(fit1)
#' plot(fit1)
#' 
#' ##
#' ## Scenario 2: standard input with sub-population specification
#' ##
#' data(RandomVA2)
#' fit2<- insilico(RandomVA2, subpop = list("sex"),  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#' 		   auto.length = FALSE)
#' summary(fit2)
#' plot(fit2, type = "compare")
#' plot(fit2, which.sub = "Men")
#' 
#' ##
#' ## Scenario 3: standard input with multiple sub-population specification
#' ##
#' fit3<- insilico(RandomVA2, subpop = list("sex", "age"),  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#' 		   auto.length = FALSE)
#' summary(fit3)
#' 
#' ##
#' ## Scenario 3: standard input with multiple sub-population specification
#' ##
#' fit3<- insilico(RandomVA2, subpop = list("sex", "age"),  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#' 		   auto.length = FALSE)
#' summary(fit3)
#' 
#' ##
#' ## Scenario 5 - 7 are special situations rarely needed in practice,
#' ##   but included here for completeness. 
#' ##   The below examples use no sub-population or physician codes, 
#' ##   but specifying sub-population is still possible as in Scenario 2 - 4.
#' ## 
#' 
#' ##
#' ## Scenario 5: skipping re-estimation of conditional probabilities
#' ##
#' # Though in practice the need for this situation is very unlikely, 
#' # use only the default conditional probabilities without re-estimation
#' fit5<- insilico(RandomVA1, subpop = NULL,  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'               updateCondProb = FALSE, 
#' 		   auto.length = FALSE) 
#' summary(fit5)
#' 
#' ##
#' ## Scenario 6: modify default conditional probability matrix
#' ##
#' # Load the default conditional probability matrix 
#' data(condprob)
#' # The conditional probabilities are given in levels such as I, A+, A, A-, etc.
#' condprob[1:5, 1:5]
#' # To modify certain cells 
#' new_cond_prob <- condprob
#' new_cond_prob["elder", "HIV/AIDS related death"] <- "C"
#' # or equivalently
#' new_cond_prob[1, 3] <- "C"
#' 
#' fit6<- insilico(RandomVA1, subpop = NULL,  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'               CondProb = new_cond_prob, 
#' 		   auto.length = FALSE) 
#' # note: compare this with fit1 above to see the change induced 
#' #  by changing Pr(elder | HIV) from "C+" to "C".
#' summary(fit6)
#' 
#' ##
#' ## Scenario 7: modify default numerical values in conditional probabilities directly
#' ##
#' # Load the default conditional probability matrix 
#' data(condprobnum)
#' # The conditional probabilities are given in numerical values in this dataset
#' condprobnum[1:5, 1:5]
#' # To modify certain cells, into any numerical values you want 
#' new_cond_prob_num <- condprobnum
#' new_cond_prob_num["elder", "HIV/AIDS related death"] <- 0.004
#' # or equivalently
#' new_cond_prob_num[1, 3] <- 0.005
#' 
#' fit7<- insilico(RandomVA1, subpop = NULL,  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'               CondProbNum = new_cond_prob_num, 
#' 		   auto.length = FALSE) 
#' # note: compare this with fit1, fit5, and fit6
#' summary(fit7)
#' 
#' ##
#' ## Scenario 8: physician coding
#' ## see also the examples in physician_debias() function section
#' ##
#' # Load sample input for physicians
#' data(RandomPhysician)
#' # The symptom section looks the same as standard input
#' head(RandomPhysician[, 1:5])
#' # At the end of file, including a few more columns of physician id and coded cause
#' head(RandomPhysician[, 245:250])
#' 
#' # load Cause Grouping (if physician-coded causes are in larger categories)
#' data(SampleCategory)
#' head(SampleCategory)
#' 
#' # existing doctor codes in the sample dataset
#' doctors <- paste0("doc", c(1:15))
#' causelist <- c("Communicable", "TB/AIDS", "Maternal",
#'                "NCD", "External", "Unknown")
#' phydebias <- physician_debias(RandomPhysician, 
#'	phy.id = c("rev1", "rev2"), phy.code = c("code1", "code2"), 
#'	phylist = doctors, causelist = causelist, 
#'	tol = 0.0001, max.itr = 100)
#' 
#' fit8 <- insilico(RandomVA1, subpop = NULL,  
#'               Nsim = 1000, burnin = 500, thin = 10 , seed = 1,
#'               phy.debias = phydebias,
#'               phy.cat = SampleCategory, 
#'               phy.external = "External", phy.unknown = "Unknown",
#' 		   auto.length = FALSE) 
#' summary(fit8)
#' 
#' }
#' @export insilico
insilico <- function(data, isNumeric = FALSE, updateCondProb = TRUE, keepProbbase.level = TRUE,  CondProb = NULL, CondProbNum = NULL, datacheck = TRUE, datacheck.missing = TRUE, warning.write = FALSE, external.sep = TRUE, Nsim = 4000, thin = 10, burnin = 2000, auto.length = TRUE, conv.csmf = 0.02, jump.scale = 0.1, levels.prior = NULL, levels.strength = 1, trunc.min = 0.0001, trunc.max = 0.9999, subpop = NULL, java_option = "-Xmx1g", seed = 1, phy.code = NULL, phy.cat = NULL, phy.unknown = NULL, phy.external = NULL, phy.debias = NULL, exclude.impossible.cause = TRUE, indiv.CI = NULL, ...){ 
	
	# handling changes throughout time
	  args <- as.list(match.call())
	  if(!is.null(args$length.sim)){
	  	Nsim <- args$length.sim
	  	cat("length.sim argument is replaced with Nsim argument, will remove in later versions.\n")
	  }

	fit <- insilico.fit(data = data, 
						isNumeric = isNumeric, 
						updateCondProb = updateCondProb, 
						keepProbbase.level = keepProbbase.level, 
						CondProb = CondProb, 
						CondProbNum = CondProbNum, 
						datacheck = datacheck, 
						datacheck.missing = datacheck.missing, 
						warning.write = warning.write, 
						external.sep = external.sep, 
						Nsim = Nsim, 
						thin = thin, 
						burnin = burnin, 
						auto.length = auto.length, 
						conv.csmf = conv.csmf, 
						jump.scale = jump.scale, 
						levels.prior = levels.prior, 
						levels.strength = levels.strength, 
						trunc.min = trunc.min, 
						trunc.max = trunc.max, 
						subpop = subpop, 
						java_option = java_option, 
						seed = seed, 
						phy.code = phy.code, 
						phy.cat = phy.cat, 
						phy.unknown = phy.unknown, 
						phy.external = phy.external, 
						phy.debias = phy.debias, 
						exclude.impossible.cause = exclude.impossible.cause, 
						indiv.CI = indiv.CI)
	return(fit)  	
} 
