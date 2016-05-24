#' Calculate a species' standardized trait field
#'
#' Calculate the null-model standardized effect size of a species' trait field.
#'
#' @param trait.distance Symmetrical matrix summarizing pairwise trait distances.
#' @param tree Phylo object.
#' @param picante.cdm A picante-style community data matrix with sites as rows, and
#' species as columns.
#' @param metric Phylogenetic metric of choice (see details).
#' @param null Null model of choice (see details).
#' @param randomizations The number of times the input CDM should be randomized and the
#' metrics calculated across it.
#' @param distances.among A symmetric distance matrix, summarizing the distances among all
#' plots from the cdm. For use with the dispersal null.
#' @param abundance.matters Default is TRUE. If FALSE, species are sampled from
#' neighboring grid cells with equal probability. For use with the dispersal null.
#' @param abundance.assigned For use with the dispersal null. See details there. 
#' @param cores This function can run in parallel. In order to do so, the user must
#' specify the desired number of cores to utilize. The default is "seq", which runs the
#' calculations sequentially.
#' 
#' @details The trait distance matrix should be symmetrical and "complete". See example.
#' Currently only non-abundance-weighted mean pairwise and interspecific
#' abundance-weighted mean pairwise phylogenetic distances are implemented. The
#' only null models that are currently implemented are the richness and dispersal nulls.
#' The function could be improved by tapping into any of the metrics and nulls defined
#' in defineMetrics and defineNulls.
#'
#' @return Data frame of standardized effect sizes of species' trait fields. Table
#' includes the observed trait field, the mean and standard deviation of the species'
#' trait field after randomization with the chosen null model, and the resulting
#' species-specific standarized effect size.
#'
#' @export
#'
#' @importFrom foreach foreach %dopar%
#' @importFrom doParallel registerDoParallel
#'
#' @references Miller, Wagner, Harmon & Ricklefs. In review. Radiating despite a lack of
#' character: closely related, morphologically similar, co-occurring honeyeaters have
#' diverged ecologically.
#'
#' @examples
#' #simulate tree with birth-death process
#' tree <- geiger::sim.bdtree(b=0.1, d=0, stop="taxa", n=50)
#'
#' #simulate trait evolution up the tree. Make 2-d trait space and find distances between
#' #species in that space
#' traits <- evolveTraits(tree)
#'
#' #calculate the distances betweeen species
#' dists <- as.matrix(dist(traits[[2]], diag=TRUE, upper=TRUE))
#'
#' #simulate log-normal abundances
#' sim.abundances <- round(rlnorm(5000, meanlog=2, sdlog=1)) + 1
#'
#' #simulate a community data matrix with these inputs
#' cdm <- simulateComm(tree, richness.vector=10:25, abundances=sim.abundances)
#'
#' #example trait field calculations
#' exampleField <- sesTraitField(trait.distance=dists, tree=tree, picante.cdm=cdm, 
#' 	metric="naw.mpd", null="richness", randomizations=10)

sesTraitField <- function(trait.distance, tree, picante.cdm, metric, null, randomizations,
	distances.among=NULL, abundance.matters=TRUE, abundance.assigned="directly",
	cores="seq")
{
	#calculate the observed trait field
	observed <- traitField(trait.distance, picante.cdm, metric)
	
	#previously set up a temporary matrix to save results into, where each column was a 
	#species and each row is a different randomization. still imitate that behavior below
	#with the foreach parallel loop, but matrix makes itself instead with the rbind
	#combine behavior of foreach. 

	#randomize the CDM with the null model, recalculate species' traitFields, then
	#save into the relevant row of the tempMatrix
	if(null == "dispersal")
	{
		#if cores is set to sequential, do not run in parallel
		if(cores == "seq")
		{
			tempMatrix <-
			foreach(i=1:randomizations, .combine='rbind') %do%
			{
				tempCDM <- dispersalNull(picante.cdm=picante.cdm, tree=tree,
					distances.among=distances.among, abundance.matters=abundance.matters,
					abundance.assigned=abundance.assigned)
				traitField(trait.distance, tempCDM, metric)
			}
		}

		#if cores is not set to seq, run in parallel
		if(cores != "seq")
		{
			#register parallel backend
			registerDoParallel(cores)

			tempMatrix <-
			foreach(i=1:randomizations, .combine='rbind') %dopar%
			{
				tempCDM <- dispersalNull(picante.cdm=picante.cdm, tree=tree,
					distances.among=distances.among, abundance.matters=abundance.matters,
					abundance.assigned=abundance.assigned)
				traitField(trait.distance, tempCDM, metric)
			}
			registerDoSEQ()
		}
	}
	else if(null == "richness")
	{
		#if cores is set to sequential, do not run in parallel
		if(cores == "seq")
		{
			tempMatrix <-
			foreach(i=1:randomizations, .combine='rbind') %do%
			{
				tempCDM <- picante::randomizeMatrix(samp=picante.cdm,
					null.model="richness")
				traitField(trait.distance, tempCDM, metric)
			}
		}

		#if cores is not set to seq, run in parallel
		if(cores != "seq")
		{
			#register parallel backend
			registerDoParallel(cores)

			tempMatrix <-
			foreach(i=1:randomizations, .combine='rbind') %dopar%
			{
				tempCDM <- picante::randomizeMatrix(samp=picante.cdm,
					null.model="richness")
				traitField(trait.distance, tempCDM, metric)
			}
			registerDoSEQ()
		}
	}
	else
	{
		stop("null must currently be set to 'dispersal' or 'richness'")
	}

	#find the mean and standard deviation of each column (species)
	metric.mean <- apply(tempMatrix, 2, mean, na.rm=TRUE)
	metric.sd <- apply(tempMatrix, 2, sd, na.rm=TRUE)
	
	#bind these to the observed values, calculate SES values and return a data frame
	results <- data.frame(observed, metric.mean, metric.sd)
	results$SES <- (results$observed-results$metric.mean)/results$metric.sd

	results
}
