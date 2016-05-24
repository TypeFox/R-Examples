################################################################################
# TODO LIST
# TODO: ...

# NOTE:
# NB! pcr.prob fixed to 1 for the PCR simulation. Only the estimated number of
# molecules is currently affected in the calibration process as the PCR efficiency
# is reduced by step.size.
# Reducing both the estimate and the pcr.prob result in constant simulated peak height...
# Although this approach will not give the exact pcr.prob it give close to expected peak heights.
# IF the number of molecules known (i.e. picked cells) then it should be possible to reduce
# pcr.prob until data fits.

################################################################################
# CHANGE LOG (10 last changes)
# 14.04.2016: Version 1.0.0 released.
# 14.04.2016: Added detection threshold range.
# 09.09.2014: First version.

#' @title Calibrate PCRsim
#'
#' @description
#' Estimates the detection threshold, peak height scaling factor,
#' and PCR efficiency needed to calibrate PCRsim.
#'
#' @details
#' To calibrate the PCR simulator perform the following experiments:
#' 1) Prepare single source samples with optimal amount of DNA (different or replicates)
#' Calculate the average total peak height (sum of peak heights) across all samples
#' to set the \code{target} parameter.
#' 2) Prepare a serial dilution (preferrably from intact cells, or from single
#' source crime scene samples). It is suitable to go from approximately optimal amount
#' down to low concentrations with drop-outs or completely blank profiles.
#' Quantify each dilution as accurately as possible. Amplify using normal 
#' procedure and analyse the PCR product.
#' Require a dataset with sample names ('Sample.Name'), and (average) concentration or
#' (average) amount in columns named 'Concentration' and 'Amount' respectively. Also
#' requires the average peak height ('H').
#' 
#' NB! Samples with either zero average peak height or zero number of molecules is removed automatically.
#' In addition, if replicate quants, samples where one replicate was negative can be removed manually.

#' @param data data.frame with observed DNA result. Required columns are sample names ('Sample.Name'),
#' average peak height ('H'), and (mean) DNA concentration ('Mean').
#' @param target integer average peak height ('H') as observed average across replicate
#'  analyses of samples with \code{dna.amount} of DNA.
#'  If NULL it will be estimated from the linear regression at \code{dna.amount} of DNA.
#' @param ref data.frame with known profiles for the samples in \code{data}. 
#' Required columns are 'Sample.Name', 'Marker', and ('Allele'), and (mean) DNA concentration ('Mean').
#' @param quant data.frame with (average) 'Concentration' or 'Amount' for the samples in \code{data}. 
#' @param ignore.case logical \code{TRUE} to ignore case in sample name matching. 
#' @param kit character string to specify the STR DNA typing kit to simulate.
#' If NULL all markers in \code{db} will be used.
#' @param db data.frame with the allele frequency database (if random profiles are simulated).
#' @param fixed.profile data.frame with columns 'Marker' and 'Allele' (if fixed profiles are simulated).
#' @param sim integer the number of simulations per calibration cycle.
#' @param pcr.cyc integer the number of PCR cycles.
#' @param pcr.aliq integer the aliquot DNA extract transferred to the PCR reaction.
#' @param pcr.vol integer the total PCR reaction volume.
#' @param ce.aliq integer the aliquot PCR product used for capillary electrophoresis.
#' @param minimize logical \code{TRUE} stops when the squared difference is minimized,
#'  \code{FALSE} continues until the PCR efficiency is 0.
#' @param step.size numeric size of PCR efficiency reduction for each calibration cycle.
#' @param cell.dna numeric the DNA content of a diploid cell in nanograms (default is 0.006 ng).
#' @param filter logical \code{TRUE} to retrieve known alleles defined in \code{ref} from \code{data}.
#' @param dna.amount numeric the amount of DNA in nanograms (ng) to be used in simulation.
#' NB! must correspond to the amount in samples used to calculate \code{target}.
#' @param plot.data logical to show linear regression data plot with marked target.
#' @param decimals interger for number of decimal places in plot title.
#' @param debug logical to print debug information.
#' @param ext.debug logical to print extended debug information.
#' 
#' @importFrom ggplot2 ggplot aes_string aes ggtitle labs geom_point geom_smooth geom_segment
#'  geom_rect annotate scale_x_continuous scale_y_continuous
#' @importFrom stats lm coefficients quantile
#' @importFrom utils head tail str
#' 
#' @export

calibratePCRsim <- function(data, target=NULL, ref=NULL, quant=NULL, ignore.case=TRUE,
                      kit=NULL, db=NULL, fixed.profile=NULL, sim=1,
                      pcr.cyc, ce.aliq=1, pcr.aliq=17.5, pcr.vol=25,
                      minimize=TRUE, step.size=0.001, cell.dna=0.006, dna.amount=0.5, filter=TRUE,
                      plot.data=TRUE, decimals=4, debug=FALSE, ext.debug=FALSE){

  # Initialise the PCR efficiency variable to 1 (100%).
  pcreff <- 1.0
  
  # Set metric to correct column name.
  metric <- "H"
  
  # PREPARE ###################################################################

  # CALCULATE PEAK HEIGHT METRICS ---------------------------------------------
  
  # Check if peak height metrics are available.
  if(!metric %in% names(data)){

    if(!is.null(ref)){

      message(paste("Metric", metric, "is not available in data."))
      
      if(filter){

        message(paste("Filter known alleles from 'data'."))
        
        # Filter the data.
        data <- strvalidator::filterProfile(data=data, ref=ref,
                                            add.missing.loci=TRUE,
                                            keep.na=TRUE,
                                            ignore.case=ignore.case,
                                            invert=FALSE, debug=ext.debug)
        
      }
      
      # Check if data has zygosity indicator.
      if(!'Heterozygous' %in% names(data)){

        message(paste("Adding Heterozygous indicator."))
        
        if(!'Heterozygous' %in% names(ref)){
          
          # Add indicator to reference dataset.
          ref <- strvalidator::calculateHeterozygous(data=ref, debug=debug)
          
          # Add indicator to dataset.
          data <- strvalidator::addData(data=data, new.data=ref, by.col="Sample.Name",
                                        then.by.col="Marker", exact=FALSE,
                                        ignore.case=ignore.case, debug=ext.debug)
          
        } else {
          
          # Add indicator to dataset.
          data <- strvalidator::addData(data=data, new.data=ref, by.col="Sample.Name",
                                        then.by.col="Marker", exact=FALSE,
                                        ignore.case=ignore.case, debug=ext.debug)
          
        }
        
      }

      message(paste("Calculating metric ", metric, ".", sep=""))
      
      # Calculate average peak height.
      data <- strvalidator::calculateHeight(data=data, na=0, add=FALSE, exclude=NULL,
                                            debug=ext.debug)
      
      # Add amount to dataset.
      data <- strvalidator::addData(data=data, new.data=quant, by.col="Sample.Name",
                                    then.by.col=NULL, exact=FALSE,
                                    ignore.case=ignore.case, debug=ext.debug)
      
      
    } else {
      
      stop("'ref' is required to calculate peak height metrics.")
      
    }
    
    
  } else {
    
    message(paste(metric, "is available in 'data'."))
    
  }
  
  # Make new column with metric to simplify code.
  data$Metric <- data[ , metric]
  
  # CALCULATE NUMBER OF CELLS -------------------------------------------------

  if("Amount" %in% names(data)){

    if(!is.numeric(data$Amount)){
      data$Amount <- as.numeric(data$Amount)
    }

    message(paste("Using 'Amount' /", cell.dna, "to estimate number of cells."))
    
    # Estimate the number of cells.
    data$Cells <- data$Amount / cell.dna
    
  } else if("Concentration" %in% names(data)){

    if(!is.numeric(data$Concentration)){
      data$Concentration <- as.numeric(data$Concentration)
    }
    
    message(paste("Using ( 'Concentration' *", pcr.aliq, ")/", cell.dna,
                  "to estimate number of cells."))
    
    # Estimate number of cells.
    data$Cells <- (data$Concentration * pcr.aliq) / cell.dna
    
  } else {
    
    stop("'data' must contain a column 'Amount' or 'Concentration'.")
    
  }
  
  if(debug){
    print("data:")
    print(head(data))
    print(tail(data))
    print(str(data))
  }

  # ESTIMATE EMPIRICAL THRESHOLD ----------------------------------------------
  
  # EMPIRICAL DETECTION THRESHOLD
  # First get the minimum number of cells resulting in a peak height > 0 (any peak > LDT).
  # Remove samples with zero peak height.
  dfmin <- data[data$Metric > 0, ]
  # Get minimum number of cells.
  minM <- round(suppressWarnings(min(dfmin$Cells)), 2)
  
  # Second get the maximum number of cells resulting in a blank profile (all peaks < LDT).
  # Remove samples with non-zero peak height.
  dfmax <- data[data$Metric == 0, ]
  # Get minimum number of cells.
  maxM <- round(suppressWarnings(max(dfmax$Cells)), 2)
  
  # Print reulst.
  message(paste("Estimated empirical detection threshold (@ LDT) is between:",
              minM, "and", maxM, "cells."))

  # DISCARD NULL RESULTS ------------------------------------------------------
  
  # Remove samples with zero peak height.
  if(any(data$Metric == 0)){
    tmp1 <- nrow(data)
    data <- data[data$Metric != 0, ]
    tmp2 <- nrow(data)
    message(paste("Removed", (tmp1 - tmp2), "rows with", metric, "= 0"))
  }
  
  # Remove samples with zero cells.
  if(any(data$Cells == 0)){
    tmp1 <- nrow(data)
    data <- data[data$Cells != 0, ]
    tmp2 <- nrow(data)
    message(paste("Removed", tmp1 - tmp2, "rows with Cells = 0"))
  }

  # CALCULATE PARAMETERS ------------------------------------------------------
  
  # Calculate the number of cells to use in simulation.
  targetCells <- dna.amount / cell.dna

  # Estimate target from regression.
  if(is.null(target)){
    
    # Regression: log peak height by log amount.
    a_model <- lm(formula=log(Metric) ~ log(Amount), data=data)
    a_intercept <- coefficients(a_model)[1]
    a_slope <- coefficients(a_model)[2]
    # a_fit <- summary(a_model)
    # a_sigma <- a_fit$sigma
    target <- round(exp(as.numeric((a_intercept + a_slope * log(dna.amount)))),0)
    message(paste("Parameter estimated from linear regression: target =", target, "RFU"))

  }
  
  # CALCULATE -----------------------------------------------------------------

  # Pre-allocate result vectors.
  res_pcreff <- rep(NA, sim)
  res_t_intercept <- rep(NA, sim)
  res_t_slope <- rep(NA, sim)
  res_t_sigma <- rep(NA, sim)
  res_intercept <- rep(NA, sim)
  res_slope <- rep(NA, sim)
  res_sigma <- rep(NA, sim)
  res_obsH <- rep(NA, sim)
  res_simH <- rep(NA, sim)
  res_sqd <- rep(NA, sim)

  # Initiate variables.
  element <- 1          # Vector index.
  previous.sqd <- NULL  # Previous squared difference.

  repeat{ # Alt. put repeat in function, f(pcreff), with output diff^2
    # Also think about confidence interval and variation compared to observed. m(pcreff) = m of pcreff = average TPH.
    # M=E[m(pcreff)]?  Calculate the number of simulations needed for confidence in the result.
    
    # Progress.
    message(paste("Simulating with PCR efficiency:", pcreff))

    # Calculate number of molecules after 'pcr.cyc' cycles PCR Nc = N0 * E^c,
    # where Nc is the number of molecules after c cycles, N0 is the number of starting molecules,
    # and E is the PCR efficiency (a number {1-2} where 2=100% efficiency), hence the 1 + pcreff.
    data$Molecules <- data$Cells * (1 + pcreff)^pcr.cyc
    
    # Calculate number of molecules loaded to CE (e.g. 1 uL from the 25uL PCR reaction).
    data$MoleculesCE <- data$Molecules * (ce.aliq / pcr.vol)
    
    # Logistic regression -----------------------------------------------------
    #    To calculate (threshold) number of molecules from peak height.
    
    # Regression: log peak height by log number of molecules.
    t_model <- lm(formula=log(MoleculesCE) ~ log(Metric), data=data)
    t_intercept <- coefficients(t_model)[1]
    t_slope <- coefficients(t_model)[2]
    t_fit <- summary(t_model)
    t_sigma <- t_fit$sigma
    
    if(debug){
      print(paste("T intercept=", t_intercept,
                  ", T slope=", t_slope,
                  ", T sigma=", t_sigma))
      print(paste("Threshold (number of molecules) @ LDT=1 RFU",
                  round(exp(as.numeric((t_intercept + t_slope * log(1)))),0)))
    }
    
    # Calculate the 5th and 95th percentile for the threshold.
    t_range_mol <- quantile(floor(exp(rnorm(100000, t_intercept, t_sigma))), c(0.05, 0.95))
    t_range_cell <- c(NA,NA)
    t_range_cell[1] <- (t_range_mol[1] / (ce.aliq / pcr.vol)) / (1 + pcreff)^pcr.cyc
    t_range_cell[2] <- (t_range_mol[2] / (ce.aliq / pcr.vol)) / (1 + pcreff)^pcr.cyc
    t_range_amount <- c(NA,NA)
    t_range_amount[1] <- t_range_cell[1] * cell.dna
    t_range_amount[2] <- t_range_cell[2] * cell.dna
    t_range_df <- data.frame(xmin = t_range_amount[1],
                             xmax = t_range_amount[2],
                             ymin = -Inf, ymax = Inf,
                             color = "red",
                             Amount=0.5, Metric=50)
    
    if(debug){
      print("Threshold range (5th and 95th percentile)")
      print(paste("Molecules in CE:", paste(floor(t_range_mol), collapse=", ")))
      print(paste("Input cells:", paste(t_range_cell, collapse=", ")))
      print(paste("Input Amount:", paste(t_range_amount, collapse=", ")))
    }

    # Logistic regression -----------------------------------------------------
    # To calculate peak height from number of molecules.
    
    # Regression: log average peak height by log number of molecules .
    model <- lm(log(Metric) ~ log(MoleculesCE), data=data)
    intercept <- coefficients(model)[1]
    slope <- coefficients(model)[2]
    fit <- summary(model)
    sigma <- fit$sigma

    if(debug){
      print(paste("Scaling intercept=", intercept, ", slope=", slope, ", sigma=", sigma))
    }
    
    # Simulate ----------------------------------------------------------------
    
    # Simulate profiles.
    dfsim <- suppressMessages(simProfile(data=fixed.profile, kit=kit, sim=sim, db=db))
    
    if(debug){
      print("Simulated profiles")
      print(head(dfsim))
      print(tail(dfsim))
    }
    
    # Simulate sample.
    dfsim <- suppressMessages(simSample(data=dfsim, cells=targetCells, sd.cells=0,
                                        conc=NULL, sd.conc=0,
                                        vol=NULL, sd.vol=0,
                                        cell.dna=cell.dna, haploid=FALSE,
                                        debug=ext.debug))
    
    # Simulate extraction.
    dfsim <- suppressMessages(simExtraction(data=dfsim, vol.ex=pcr.aliq, sd.vol=0,
                                            prob.ex=1, sd.prob=0,
                                            cell.dna=cell.dna, debug=ext.debug))
    
    # Simulate PCR.
    dfsim <- suppressMessages(simPCR(data=dfsim, kit=NULL, pcr.cyc=pcr.cyc,
                                     pcr.prob=1, sd.pcr.prob=0,
                                     stutter.prob=0, sd.stutter.prob=0,
                                     vol.aliq=pcr.aliq, sd.vol.aliq=0,
                                     vol.pcr=pcr.vol, sd.vol.pcr=0,
                                     stutter=FALSE, debug=ext.debug))
    
    if(debug){
      print("Simulated data after PCR:")
      print(head(dfsim))
      print(tail(dfsim))
    }
    
    # Combine homozygote alleles.
    dfsimph <- suppressMessages(compact(data=dfsim, per.sample=TRUE, col="PCR.Amplicon", sim=TRUE))
    
    
    # Simulate capillary electrophoresis.
    dfsimph <- suppressMessages(simCE(data=dfsimph, vol=ce.aliq, sd.vol=0,
                                      intercept=intercept, slope=slope, sigma=sigma,
                                      t.intercept=t_intercept, t.slope=t_slope, t.sigma=t_sigma,
                                      debug=ext.debug))

    if(debug){
      print("Simulated data after CE:")
      print(head(dfsimph))
      print(tail(dfsimph))
    }
    
    # Analyse -----------------------------------------------------------------
    
    # Add heterozygote flags.
    dfsimph <- strvalidator::calculateHeterozygous(data=dfsimph, debug=ext.debug)

    if(debug){
      print("Simulated data after adding heterozygous indicator:")
      print(head(dfsimph))
      print(tail(dfsimph))
    }
    
    # Calculate peak height metrics.
    dfsimH <- strvalidator::calculateHeight(data=dfsimph, na=NULL, add=FALSE, debug=ext.debug)

    # Calculate total mean.
    dfSimMean <- mean(dfsimH[ , metric], na.rm=TRUE)
      
    # Calculate squared differences.
    sqd <- (dfSimMean - target)^2
    
    if(debug){
      print("Simulated peak height metrics:")
      print(head(dfsimH))
      print(tail(dfsimH))
      print(str(dfsimH))
      print(paste("Simulated mean:", dfSimMean))
      print(paste("Target mean:", target))
      print(paste("Squared Difference:", sqd))
    }
    
    # Create result -------------------------------------------------------------
    
    # Add in vectors.
    res_pcreff[element] <- pcreff
    res_t_intercept[element] <- t_intercept
    res_t_slope[element] <- t_slope
    res_t_sigma[element] <- t_sigma
    res_intercept[element] <- intercept
    res_slope[element] <- slope
    res_sigma[element] <- sigma
    res_obsH[element] <- target  
    res_simH[element] <- dfSimMean
    res_sqd[element] <- sqd
    
    # Prepare next loop.
    pcreff <- pcreff - step.size
    element <- element + 1

    # Initiate of NULL.
    if(is.null(previous.sqd)){
      previous.sqd <- sqd
    }
    
    # Stop when PCR efficiency is below 0 (0%).
    if(pcreff < 0){
      break
    }
    
    # If stop when sqd is minimized.
    if(minimize){
      if(sqd > previous.sqd){
        break
      }
    }
    
    # Store previous value.
    previous.sqd <- sqd
    
  } # End repeat.
  
  # Create result dataframe.
  res <- data.frame("PCR.Efficiency"=res_pcreff,
                    "T.intercept"=res_t_intercept,
                    "T.slope"=res_t_slope,
                    "T.sigma"=res_t_sigma,
                    "Intercept"=res_intercept,
                    "Slope"=res_slope,
                    "Sigma"=res_sigma,
                    "Target"=res_obsH,
                    "Simulated"=res_simH,
                    "Sqd"=res_sqd,
                    stringsAsFactors=FALSE)
  
  # Remove na rows.
  res <- res[!is.na(res$PCR.Efficiency),]
  
  if(plot.data){
    
    a_title <- paste("Calibration data with target peak height (RFU) and amount (ng)",
                     "\n[PCR.Efficiency=", res_pcreff[which.min(res_sqd)], "]",
                     "\n[Detection threshold: T.intercept=", round(unique(res_t_intercept), decimals),
                     "T.slope=", round(unique(res_t_slope), decimals),
                     "T.sigma=", round(unique(res_t_sigma), decimals), "]",
                     "\n[Peak height scaling: intercept=", round(unique(res_intercept), decimals),
                     "slope=", round(unique(res_slope), decimals),
                     "sigma=", round(unique(res_sigma), decimals), "]")
    # Plot data.
    # NB! A workaround for 'object not found' bug is to use environment = environment():
    # http://stackoverflow.com/questions/5106782/use-of-ggplot-within-another-function-in-r
    a_plot <- ggplot(data=data, aes_string(x="Amount", y="Metric"), environment = environment())
    a_plot <- a_plot + labs(title = a_title)
    a_plot <- a_plot + labs(y = "Average peak height")
    a_plot <- a_plot + geom_point(shape=1)
    a_plot <- a_plot + geom_smooth(method='lm', fullrange = TRUE)
    a_plot <- a_plot + geom_segment(aes(x = 0, y = target, xend = dna.amount, yend = target), linetype="dotted")
    a_plot <- a_plot + geom_segment(aes(x = dna.amount, y = 0, xend = dna.amount, yend = target), linetype="dotted")
    a_plot <- a_plot + annotate("text", label = target, x = dna.amount, y = target,
                                colour = "red", hjust = 0, vjust = 0)
    a_plot <- a_plot + annotate("text", label = dna.amount, x = dna.amount, y = 0,
                                colour = "red", hjust = 0, vjust = 0)
    a_plot <- a_plot + geom_rect(data=t_range_df, 
                                 aes(xmin=t_range_df$xmin, xmax=t_range_df$xmax, ymin=0, ymax=Inf),
                                 fill="red", alpha=0.2)
    a_plot <- a_plot + scale_x_continuous(trans='log2')
    a_plot <- a_plot + scale_y_continuous(trans='log2')
    print(a_plot)
    
  }
  
  # Print info about threshold:
  message("Threshold range (5th and 95th percentile)")
  message(paste("Molecules in CE:", paste(floor(t_range_mol), collapse=", ")))
  message(paste("Input cells:", paste(t_range_cell, collapse=", ")))
  message(paste("Input Amount:", paste(t_range_amount, collapse=", ")))
  
  return(res)
  
}