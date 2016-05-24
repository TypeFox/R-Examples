#' Simulate \eqn{g2} 
#' 
#' This function can be used to simulate genotype data, draw subsets of loci and calculate the
#' respective \eqn{g2} values. Every subset of markers is drawn independently to give insights
#' into the variation and precision of \eqn{g2} calculated from a given number of markers and individuals. 
#' 
#'
#' @param n_ind number of individuals to sample from the population
#' @param H_nonInb true genome-wide heteorzygosity of a non-inbred individual
#' @param meanF mean realized inbreeding f
#' @param varF variance in realized inbreeding f
#' @param subsets a vector specifying the sizes of marker-subsets to draw. Specifying  \code{subsets = c(2, 5, 10, 15, 20)} 
#'        would draw marker sets of 2 to 20 markers. The minimum number of markers to calculate g2 is 2.
#' @param reps number of resampling repetitions
#' @param type specifies g2 formula. Type "snps" for large datasets and "msats" for smaller datasets.
#' @param CI Confidence intervals to calculate (default to 0.95)
#'
#' @details The \code{simulate_g2} function simulates genotypes from which subsets of loci can be sampled independently. 
#'          These simulations can be used to evaluate the effects of the number of individuals 
#'          and loci on the precision and magnitude of \eqn{g2}. The user specifies the number of simulated individuals (\code{n_ind}), the subsets of 
#'          loci (\code{subsets}) to be drawn, the heterozygosity of non-inbred individuals (\code{H_nonInb}) and the 
#'          distribution of \emph{f} among the simulated individuals. The \emph{f} values of the simulated individuals are sampled
#'          randomly from a beta distribution with mean (\code{meanF}) and variance (\code{varF}) specified by the user 
#'          (e.g. as in wang2011). This enables the simulation to mimic populations with known inbreeding 
#'          characteristics, or to simulate hypothetical scenarios of interest. For computational simplicity, allele 
#'          frequencies are assumed to be constant across all loci and the simulated loci are unlinked. Genotypes
#'          (i.e. the heterozygosity/homozygosity status at each locus) are assigned stochastically based on the \emph{f}
#'          values of the simulated individuals. Specifically, the probability of an individual being heterozygous at
#'          any given locus (\eqn{H}) is expressed as \eqn{H = H0(1-f)} , where \eqn{H0} is the user-specified heterozygosity of a 
#'          non-inbred individual and \emph{f} is an individual's inbreeding coefficient drawn from the beta distribution.
#'          
#' @return
#' \code{simulate_g2} returns an object of class "inbreed".
#' The functions `print` and `plot` are used to print a summary and to plot the g2 values with means and confidence intervals
#' 
#' An `inbreed` object from  \code{simulate_g2} is a list containing the following components:
#' \item{call}{function call.}
#' \item{estMat}{matrix with all r2(h,f) estimates. Each row contains the values for a given subset of markers}
#' \item{true_g2}{"true" g2 value based on the assigned realized inbreeding values}
#' \item{n_ind}{specified number of individuals}
#' \item{subsets}{vector specifying the marker sets}
#' \item{reps}{repetitions per subset}
#' \item{H_nonInb}{true genome-wide heteorzygosity of a non-inbred individual}
#' \item{meanF}{mean realized inbreeding f}
#' \item{varF}{variance in realized inbreeding f}
#' \item{min_val}{minimum g2 value}
#' \item{max_val}{maximum g2 value}
#' \item{all_CI}{confidence intervals for all subsets}
#' \item{all_sd}{standard deviations for all subsets}
#' 
#' @author  Marty Kardos (marty.kardos@@ebc.uu.se) &
#'          Martin A. Stoffel (martin.adam.stoffel@@gmail.com) 
#'          
#' @examples 
#' data(mouse_msats)
#' genotypes <- convert_raw(mouse_msats)
#' sim_g2 <- simulate_g2(n_ind = 10, H_nonInb = 0.5, meanF = 0.2, varF = 0.03,
#'                       subsets = c(4,6,8,10), reps = 100, 
#'                       type = "msats")
#' plot(sim_g2)
#' @export


simulate_g2 <- function(n_ind = NULL, H_nonInb = 0.5, meanF = 0.2, varF = 0.03,
                        subsets = NULL, reps = 100, type = c("msats", "snps"),
                        CI = 0.95) {
################################################################################
# simulate a population with variable inbreeding
# then estimate g2 from independently sampled / non-overlapping subsets of loci 
################################################################################
    if ((H_nonInb > 1) | (H_nonInb < 0)) stop("H_nonInb has to be a value between 0 and 1")
    if ((meanF > 1) | (meanF < 0)) stop("meanF has to be a value between 0 and 1")
    if (varF < 0) stop("meanF has to be a value between 0 and 1")

# number of individuals to sample from the population
if (is.null(n_ind))   stop("Specify the number of individuals to sample with n_ind")  
# subsets of loci to sample   
if (is.null(subsets)) stop("specify the size of loci subsamples in 'subsets', i.e. subsets = c(2,4,6,8) to 
                           calculate g2 from up to 8 loci")
if (any(subsets < 2)) stop("Specify a minimum of 2 markers in subsets")
if (!isTRUE(all(subsets == floor(subsets)))) stop("'subsets' must only contain integer values")
    
    # check g2 function argument
    if (length(type) == 2){
        type <- "msats"
    } else if (!((type == "msats")|(type == "snps"))){
        stop("type argument needs to be msats or snps")
    } 
    
    # define g2 function
    if (type == "msats") {
        g2_fun <- g2_microsats
    } else if (type == "snps") {
        g2_fun <- g2_snps
    }
    
# total number of loci to simulate
n_loc <- subsets[length(subsets)]
allLoci <- reps*n_loc                            

##############################################
# sample F values from a beta distribution
# following previous work from Jinliang Wang
# (2011, Heredity 107, pp. 433-443)
##############################################

alpha <- ((1 - meanF) / varF - (1 / meanF)) * meanF ^ 2    # parameters of the beta distribution given the mean and variance of realized F
beta <- alpha * ((1 / meanF) - 1)

Fs <- stats::rbeta(n_ind, alpha, beta)   # vector of realized F for the simulated individuals
true_g2 <- stats::var(Fs)/((1-mean(Fs))^2)

#-------------------------
# simulate the individual 
# genotypes
#-------------------------

hets <- NULL        # initialize a data frame to store the genotypic information

for (i in 1:n_ind) {
    
	thisHet <- NULL                       # randomly select a TRUE genome-wide MLH (i.e., the proportion of hypothetically infinitely many loci that are heterozygous in the ith individual)
	thisHet <- H_nonInb*(1-Fs[i])

	rands <- NULL                         # randomly generated numbers between 0 and 1 that are used to determine whether the individual is heterozygous at each locus
	rands <- stats::runif(allLoci,min=0,max=1)

	theseHets <- NULL
	theseHets <- as.numeric(rands < thisHet)
	
	hets <- rbind(hets,theseHets)
	}

#------------------------------------------------------------------------
# repetitively subsample the loci independently, each time estimating g2
#------------------------------------------------------------------------

estMat <- NULL    # matrix to store teh g2 estimates
sampCols <- 1:ncol(hets)    # vector of loci that are available for sampling


estMat <- NULL

for (i in 1:length(subsets))    # loop through the differen subsample sizes
	{
	theseEsts <- rep(NA,reps)    # vector to store the estimates from this number of loci
	for (j in 1:reps)
		{
        theseSampCols <- NULL                                       # get a new independent sample of loci
		theseSampCols <- sample(sampCols,subsets[i],replace=FALSE)

		theseGenos <- NULL
		theseGenos <- hets[,theseSampCols]

		theseEsts[j] <- g2_fun(theseGenos)[2][[1]]
		sampCols <- sampCols[-theseSampCols]

		}

	estMat <- rbind(estMat,theseEsts)
	sampCols <- 1:ncol(hets)    # reconstitute the original vector of loci that are available for sampling
	print(paste("done with subsampling of ",subsets[i]," loci",sep=""))

	}
estMat <- unname(estMat)
# get the upper and lower bounds of the y-axis

minG2 <- min(estMat)
maxG2 <- max(estMat)

# calculate CIs and SDs
calc_CI <- function(estMat_subset) {
    stats::quantile(estMat_subset, c((1-CI)/2,1-(1-CI)/2), na.rm=TRUE)
}

all_CI <- t(apply(estMat, 1, calc_CI))
all_sd <- apply(estMat, 1, stats::sd)

res <- list(call=match.call(),
            estMat = estMat,
            true_g2 = true_g2,
            n_ind = n_ind,
            subsets = subsets,
            reps = reps,
            H_nonInb =  H_nonInb,
            meanF = meanF,
            varF = varF,
            min_val = minG2,
            max_val = maxG2,
            all_CI = all_CI,
            all_sd = all_sd
            )

class(res) <- "inbreed"
return(res)
}





