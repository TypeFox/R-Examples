################################################################################
# TODO LIST
# TODO: ...

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 16.03.2015: Implemented use of data.table for more efficient calculations.
# 11.09.2014: First version.
#
#' @title Calibrate Inter-locus Balance
#'
#' @description
#' Estimate values for the PCR efficiency parameter per locus that satisfy
#' the target inter-locus balances.
#'
#' @details
#' The inter-locus balance for a kit should be characterised during the
#' internal validation of the kit. The function search for PCR efficiency
#' values per locus that upon simulation are similar to the target
#' inter-locus balances. Use the PCR efficiency value obtained from the
#' \code{calibratePCRsim} function as \code{seed} value.
#' 
#' @param sim integer the number of simulations per calibration cycle.
#' @param target numeric vector with target interlocus balances.
#' @param amount numeric for amount of DNA in ng.
#' @param cell.dna numeric the DNA content of a diploid cell in nanograms (default is 0.006 ng).
#' @param pcr.cyc integer the number of PCR cycles.
#' @param acc.dev numeric, accepted deviation from target.
#' @param step.size numeric, the probability of PCR is changed by this value.
#' @param seed numeric, start value for optimisation of the PCR probability.
#' @param max.eff numeric, maximal value for estimated PCR efficiency.
#' @param progress logical, print progress to console.
#' @param debug logical to print debug information.
#' 
#' @return vector with estimated PCR efficiencies for each locus.
#' 
#' @importFrom data.table data.table :=
#' @importFrom stats var
#' @importFrom utils head
#' 
#' @export
#' 
#' @examples
#' # Experimental inter-locus balances for the STR kit to be simulated (sums to 1).
#' target <- c(0.20, 0.10, 0.15, 0.25, 0.30)
#'  
#' # Find PCR efficiency values that upon simulation
#' # satisfy the experimental data for 0.5 ng of input DNA.
#' set.seed(10) # For reproducibility.
#' calibrateLb(sim=10, target=target, amount=0.5, seed=0.85, progress=FALSE)
#' 
#' # Locus specific PCR efficency parameters can now be used as parameters.
#' # [1] 0.858 0.816 0.841 0.871 0.883 

calibrateLb <- function(sim=100, target, amount=0.5, cell.dna=0.006, pcr.cyc=30,
                        acc.dev=0.001, step.size=0.001, seed=0.85, max.eff=0.98,
                        progress=TRUE, debug=FALSE){

	# Declare variables.
	optimised <- rep(FALSE, length(target))
	lap <- 0
	simProp <- matrix(NA, nrow=sim, ncol=length(target))
  
  # Constants.
  volFixed <- 10    # Extraction volume and aliquote extract for PCR amplification.

  # Check seed.  
  if(seed > max.eff){
    seed <- max.eff
    message("'seed' can't be higher than 'max.eff'! Setting seed = max.eff = ", max.eff)
  }

	# PREPARE --------------------------------------------------------------------
	
	# Create a data frame with a DNA profile.
	markers <- rep(paste("Marker", seq(1:length(target))), each=2)
	alleles <- rep(c("Allele1","Allele2"), length(target))
	df <- data.frame(Marker=markers, Allele=alleles)

  # Simulate profile
	res <- suppressMessages(simProfile(data=df, sim=sim, name="Sample"))
  
  # Simulate sample
	res <- suppressMessages(simSample(data=res, cells=amount/cell.dna, sd.cells=0,
                                    conc=NULL, sd.conc=0, vol=NULL, sd.vol=0,
	                                  cell.dna=NULL, haploid=FALSE))
	# Simulate extraction.
  res <- suppressMessages(simExtraction(data=res, vol.ex=volFixed, sd.vol=0,
                                        prob.ex=1, sd.prob=0))

	# Start with seeded PCR probability for all loci.
	prob <- rep(seed, length(target))

	# Create and initiate kit parameter data frame.
	kitparam <- data.frame(Marker=unique(markers), PCR.Prob=prob,
                         PCR.Prob.Stutter=0,
	                       stringsAsFactors=FALSE)  
	
  # OPTIMISE ------------------------------------------------------------------
  

	# Repeat until all PCR probabilities have been optimised.
	while(!all(optimised)){
	
    # Update PCR efficiency.
 	  kitparam$PCR.Prob <- prob
    
		# Simulate using current PCR probability.
		simData <- suppressMessages(simPCR(data=res, kit=kitparam, pcr.prob=rep(prob, each=2),
                         stutter.prob=0, sd.stutter.prob=0,
                         pcr.cyc=pcr.cyc, vol.aliq=volFixed, sd.vol.aliq=0, 
                   vol.pcr=volFixed, sd.vol.pcr=0, stutter=FALSE, debug=FALSE))


    # Convert to data.table.
		dt <- data.table::data.table(simData)

    # Calculate sum peak height per locus/marker (LPH).
		simPh <- dt[,list(LPH=sum(PCR.Amplicon, na.rm=TRUE)), by=list(Sample.Name, Marker)]

		# Calculate total peak height sum per sample (TPH).
		simPh <- simPh[,TPH:=sum(LPH, na.rm=TRUE), by=list(Sample.Name)]

		# Calculate proportion for each loci (LB).
    simPh$LB <- simPh$LPH/simPh$TPH

    # Calculate mean locus balance per marker (MLB) and variance per marker (VLB) across all simulations.
		simLb <- simPh[,list(MLB=mean(LB, na.rm=TRUE), VLB=var(LB, na.rm=TRUE)), by=list(Marker)]
		
    # Compare the means to the target locus balance.
    diffLb <- simLb$MLB - target

		# Check if any simulated mean is within accepted range.
		acceptedLoci <- abs(diffLb) < acc.dev
    
    if(debug){
      print("lap")
      print(lap)
      print("simPh")
      print(head(simPh))
      print("simProp")
      print(head(simProp))
      print("simLb")
      print(simLb)
      print("diffLb")
      print(diffLb)
      print("acceptedLoci")
      print(acceptedLoci)
    }
    
		if(any(acceptedLoci)) {
			
			# Set flag for accepted loci.
			optimised[acceptedLoci] <- TRUE
			
		}

		# Optimise as long as not all is accepted.
		if(!all(acceptedLoci)) {

			# Increase counter.
			lap <- lap + 1

			# Indicate direction of change for PCR probability.
			sign <- diffLb 
			sign[diffLb <0] <- 1
			sign[diffLb >=0] <- -1
			sign[acceptedLoci] <- 0

			# Print progress.
			if(progress){
				message("Iteration#",lap)
				message("PCR probability:")
				print(prob)
				message("Simulated mean inter locus balance (", sim, " simulations):")
				print(simLb$MLB)
				message("Change PCR probability:")
				print(sign)
			}

			# Change PCR probability in the desired direction.
			prob[!acceptedLoci] <- prob[!acceptedLoci] + sign[!acceptedLoci] * step.size

			# Decrease all if any PCR probability is more than 1.
			if (any(prob > max.eff)){
				prob <- prob - step.size
			}
			# Increase all if any PCR probability is less than 0.
			if (any(prob < 0)){
				prob <- prob + step.size
			}

		} # End If optimise.

	} # End While loop.
		
	return(prob)

}
