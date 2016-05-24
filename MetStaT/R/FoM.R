# 		FoM (Figures of Merit) tool, part of the 'MetStaT' package  
#		Copyright (C) 2012 Marcel van Batenburg and Tim Dorscheidt
#		
#		This program is free software: you can redistribute it and/or modify
#		it under the terms of the GNU General Public License as published by
#		the Free Software Foundation, either version 3 of the License, or
#		(at your option) any later version.
#		
#		This program is distributed in the hope that it will be useful,
#		but WITHOUT ANY WARRANTY; without even the implied warranty of
#		MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#		GNU General Public License for more details.
#		
#		You should have received a copy of the GNU General Public License
#		along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# 		Email: g.zwanenburg@uva.nl ('MetStaT' contact person) or tdorscheidt@gmail.com
#		Article with more details:
#		New Figures of Merit for Comprehensive Functional Genomics Data: The Metabolomics Case 
#		Van Batenburg MF, Coulier L, van Eeuwijk F, Smilde AK, Westerhuis JA 
#		Analytical Chemistry, Volume 83(9) (2011), pages 3267-3274
###############################################################################

# main method to calculate and plot the figure of merit
FoM.Calculate <- function(data, cols.id = 1, col.value, ids.per.bin, quiet = FALSE, repeats.per.id = -1, fit.type = 1, alpha.steps = 100) {
	# this is a default run of the FoM method, only performing the method on a single value column
	# if you wish to perform the analysis many times, or with many different columns, it is advised to manually perform the FoM.GetIdOrdererMatrix method for all value columns, and then perform the next two methods
	id.ordered.matrix <- FoM.GetIdOrderedMatrix(data, cols.id, col.value)
	binned.ids <- FoM.OrderAndBinIdByMeans(id.ordered.matrix, 1, ids.per.bin, quiet, repeats.per.id)
	result <- FoM.FitBinnedSampleRepeatErrors(binned.ids, fit.type, alpha.steps)
	result
}

# group repeats of the same sample and find the average per sample
FoM.GetIdOrderedMatrix <- function(data, cols.id = 1, cols.values = -1) {
	# first check for common errors in input
	if (sum(cols.id %in% cols.values)>0) {stop("Non-empty intersection between sample-value columns and ID columns. Please exclude the violating column(s) from either set (ID-column is set to 1 by default).")}
	if (!is.matrix(data)) {stop("Input data is not correctly formatted. Make sure it is a numerical matrix.")}
	if (sum(dim(data[2]>cols.values))>0 || sum(dim(data[2]>cols.id))>0) {stop("One of the requested sample-value columns or ID columns is beyond the scope of the input data.")}
	
	# handle conversion of default parameter-input to workable parameter-settings
	if (length(cols.values)==1 && cols.values==-1) {
		cols.values <- 1:dim(data)[2]
		cols.values <- cols.values[-cols.id]
	}
	if (length(cols.values)==0) {stop("No sample-value column(s) present. Perhaps too many ID columns were selected.")}
	if (is.null(rownames(data))) {
		rownames(data) <- paste(sep="","r",1:dim(data)[1])
	}
	if (any(is.na(data[,cols.id]))) {stop("At least one NA present in an ID column of your data. Cannot continue.")}
		
	# on the basis of the ID columns, find out which rows are repeats of the same sample-measurement and assign these a shared new unique ID
	result <- list()
	result$ids <- ASCA.GetRowRepeats(data[,cols.id,drop=FALSE])
	# exclude all ID columns and non-required sample-value columns from input data
	result$data <- data[,cols.values,drop=FALSE]
	
	no.ids <- dim(result$ids$row.patterns)[1]
	no.values <- length(cols.values)
	
	# calculate mean per ID per value-column
	result$means.per.id <- matrix(0, nrow = no.ids, ncol = no.values)
	colnames(result$means.per.id) <- colnames(result$data)
	rownames(result$means.per.id) <- 1:no.ids
	
	for (id.index in 1:no.ids) {
		for (col.no in 1:no.values) {
			values.to.average <- result$data[ result$ids$indices.per.pattern[[id.index]] , col.no ]
			values.to.average <- values.to.average[!is.na(values.to.average)] # ignore NA's when calculating 
			if (length(values.to.average)==0) { result$means.per.id[id.index,col.no] <- NA
			} else { result$means.per.id[id.index,col.no] <- mean(values.to.average) }
		}
	}
	
	result
}

# binning step of FoM method
FoM.OrderAndBinIdByMeans <- function(id.ordered.matrix, value.col, ids.per.bin, quiet = FALSE, repeats.per.id = -1) {	
	# first, find out in which order the IDs (rows containing repeats of the same sample-measurement) should be binned together
	id.order.by.mean <- sort(id.ordered.matrix$means.per.id[ ,value.col], index.return = TRUE)$ix
	
	unordered.result <- list()
	bin.counter <- 0
	current.bin <- 0
	samples.discarded <- 0
	na.encountered <- 0
	for (id in id.order.by.mean) { # go through each ID in order of increasing mean-value
		if (bin.counter == 0 || bin.counter > ids.per.bin) { # either first bin or previous bin is full: start a new one
			current.bin <- current.bin + 1
			unordered.result[[current.bin]] <- list()
			bin.counter <- 1
		}
		current.id.values <- id.ordered.matrix$data[id.ordered.matrix$ids$indices.per.pattern[[id]] ,value.col, drop = FALSE] # get all values for the current id
		na.encountered <- na.encountered + ( length(current.id.values) - length(current.id.values[!is.na(current.id.values)])) # N
		current.id.values <- current.id.values[!is.na(current.id.values)] # filter out NAs
		
		if (length(current.id.values)>1) { # only include IDs in a bin when there are at least two repeats
			if (repeats.per.id != -1 && length(current.id.values)>repeats.per.id) { # only include a max amount of repeats per ID (or multiple IDs, if enough repeats are available)
				# too many repeats for this sample, split the sample up into multiple ranges/subsets of values and add each to a bin
				ranges <- seq(0,length(current.id.values),repeats.per.id)
				for (r in 2:length(ranges)) {
					range.start <- ranges[r-1]+1
					range.end <- ranges[r]
					unordered.result[[current.bin]][[bin.counter]] <- current.id.values[range.start:range.end] 
					bin.counter <- bin.counter + 1          
					if (bin.counter > ids.per.bin) { # either first bin or previous bin is full: start a new one
						current.bin <- current.bin + 1
						unordered.result[[current.bin]] <- list()
						bin.counter <- 1
					}     
				}
			}
			else { # include all repeats in single ID
				unordered.result[[current.bin]][[bin.counter]] <- current.id.values
				bin.counter <- bin.counter + 1
			}
		}
		else {
			samples.discarded <- samples.discarded + 1
		}
	}
	if (length(unordered.result[[current.bin]])==0) { # check whether the last bin created actually received any valid entries (otherwise: delete)
		unordered.result[[current.bin]] <- NULL
		current.bin <- current.bin - 1
	}
	if (current.bin<2) {print("Too few repeating sample were found to create 2 bins or more. Please adjust bin-size or select a column with more sample values."); return(NULL)}
	if (!quiet && samples.discarded > 0) {
		print(paste("During binning,",samples.discarded,"samples had no further repeats, and were therefore discarded."))
	}
	if (!quiet && na.encountered > 0) {
		print(paste("During binning,",na.encountered,"values were NA, and were ignored."))
	}
	# per bin, calculate mean, and calculate the mean variance over all IDs in that bin
	unordered.result$mean.per.bin <- array(0,dim=c(current.bin))
	unordered.result$error.per.bin <- array(0,dim=c(current.bin))
	for (bin in 1:current.bin) {
		variances.per.id <- c()
		for (samples.per.id in unordered.result[[bin]]) {
			mean.this.id <- mean(samples.per.id)
			variances.per.id <- c(variances.per.id,
					sum( (samples.per.id-mean.this.id)^2 ) / (length(samples.per.id)-1)  
			)
		}
		unordered.result$mean.per.bin[bin] <- mean(unlist(unordered.result[[bin]]))
		unordered.result$error.per.bin[bin] <- mean(variances.per.id)
	}
	ordered.result <- list()
	# check whether the order of mean.per.bin is still correct (value removal and NA handling can affect mean intensity)
	if (sum(order(unordered.result$mean.per.bin)!=1:length(unordered.result$mean.per.bin)) > 0) {
		correct.order <- order(unordered.result$mean.per.bin)
		for (b in 1:current.bin) {
			ordered.result[[b]] <- unordered.result[[correct.order[b]]] 
		}
		ordered.result$mean.per.bin <- unordered.result$mean.per.bin[correct.order] 
		ordered.result$error.per.bin <- unordered.result$error.per.bin[correct.order]
	} else {
		ordered.result <- unordered.result 
	}
	
	ordered.result$name.of.value.col <- colnames(id.ordered.matrix$data)[value.col]
	ordered.result
}

# method to iterativaly find the best fit per alpha value and plot/return the overall best fit
FoM.FitBinnedSampleRepeatErrors <- function(ordered.bins.by.mean, fit.type = 1, alpha.steps = 100) {
	result <- ordered.bins.by.mean
	# define requested fit function
	if (fit.type==1) {
		local.fit.function <- function(...) {
			lm(...)
		} 
	} else if (fit.type==2) { 

		local.fit.function <- function(...) {
			prevWarnOption <- options(warn=-1)
			output <- rlm(...,maxit=800)
			options(prevWarnOption)
			# if the robust regression method does not converge, use classic regression instead 
			if (is.na(output) || is.na(output$coefficients) || length(output$coefficients)==0) {
				print("The regression method chosen did not converge. Using classical regression instead.")
				output <- lm(...)	
			}
			output
		}	
	} else { stop("Selected fit type is not valid.")}
	# alpha is the mean-intensity (x-axis) at which the additive fit stops and the multiplicative fit starts
	
	alphas <- seq(from=min(ordered.bins.by.mean$mean.per.bin), to=max(ordered.bins.by.mean$mean.per.bin), length.out=alpha.steps)
	result$alphas <- alphas
	ssq.res.per.alpha <- c() # used to store the sum of squared residuals per alpha
	# no matter the precise alpha; there are only as many additive regression possibilities as there are bins (minus 1), so find these first:
	no.bins <- length(ordered.bins.by.mean$mean.per.bin)
	ad.reg.coeff.per.bin <- array(0,no.bins) # each of these additive regressions possibilites has a unique coefficient 
	ad.reg.res.ssq.per.bin <- array(0,no.bins) # orderedBinsByMean$errorPerBin^2 # each of these additive regressions possibilites has a unique sum of squared residuals (
	for (ad in 2:no.bins) {
		additive.fit <- local.fit.function(ordered.bins.by.mean$error.per.bin[1:ad] ~ 1)
		ad.reg.coeff.per.bin[ad] <- additive.fit$coefficients[1]
		ad.reg.res.ssq.per.bin[ad] <- sum((additive.fit$residuals)^2)
	}
	# Now for the multiplicative regression. Although the origin is different per alpha, the number of bin-sets to use for this regression is limited to the number of bins (minus 1)
	ad.coeff.per.alpha <- array(NA, length(alphas))
	mu.coeff.per.alpha <- array(NA, length(alphas))
	tot.res.ssq.per.alpha <- array(NA, length(alphas))
	
	for (mu in 1:(no.bins-1)) {
		# each of these bin-sets to be used for the multiplicative regressions has a new origin. The mean value of the origin is equal to alpha, but the error value is equal to the regression-coefficient of the additive regression
		current.alphas.indices.to.use <- which(alphas<=ordered.bins.by.mean$mean.per.bin[mu]) # which values of alpha should be used for the current multiplicative regression
		
		if (length(current.alphas.indices.to.use) > 0) {			
			error.value.of.origin <- 0 # by default, the origin is set to zero if no additive regression is done (when the regression was skipped due to not having enough bins (0 or 1 bin))
			if (mu > 2) {  # an additive regression was done, so find its error-value
				error.value.of.origin <- ad.reg.coeff.per.bin[mu-1]
				ad.coeff.per.alpha[current.alphas.indices.to.use] <- ad.reg.coeff.per.bin[mu-1]
				#print(current.alphas.indices.to.use)
				#print(ad.reg.coeff.per.bin[mu-1])
			}
			
			error.per.bin.for.new.origin <- ordered.bins.by.mean$error.per.bin[mu:no.bins] - error.value.of.origin
			
			for (alpha.index in current.alphas.indices.to.use) {
				mean.per.bin.for.new.origin <- ordered.bins.by.mean$mean.per.bin[mu:no.bins] - alphas[alpha.index]
				multiplicative.fit <- local.fit.function(error.per.bin.for.new.origin ~ mean.per.bin.for.new.origin + 0)
				mu.coeff.per.alpha[alpha.index] <- multiplicative.fit$coefficients[1]
				tot.res.ssq.per.alpha[alpha.index] <- sum((multiplicative.fit$residuals)^2)
				if (mu > 2) {tot.res.ssq.per.alpha[alpha.index] <- tot.res.ssq.per.alpha[alpha.index] + ad.reg.res.ssq.per.bin[mu-1]}
			}
			alphas[current.alphas.indices.to.use] <- NA # exclude the currently used alphas from any future multiplicative regressions
		}
	}
	# the only remaining fit is the exclusively additive fit
	current.alphas.indices.to.use <- which(!is.na(alphas))
	ad.coeff.per.alpha[current.alphas.indices.to.use] <- ad.reg.coeff.per.bin[no.bins]
	tot.res.ssq.per.alpha[current.alphas.indices.to.use] <- ad.reg.res.ssq.per.bin[no.bins]
	
	result$tot.res.ssq.per.alpha <- tot.res.ssq.per.alpha
	result$ad.coeff.per.alpha <- ad.coeff.per.alpha
	result$mu.coeff.per.alpha <- mu.coeff.per.alpha
	
	# plot mean intensity vs. (mean) variance per bin, and plot the best fit
	if (is.null(ordered.bins.by.mean$name.of.value.col) || is.na(ordered.bins.by.mean$name.of.value.col)) {plot.name <- "Best fit for FoM"} else {plot.name <- paste(sep="","Best fit for FoM (",ordered.bins.by.mean$name.of.value.col,")")}
	plot(result$mean.per.bin,result$error.per.bin,xlab="Intensity (mean of binned sample-repeats)", ylab = "Error (mean of sd per sample in bin)", main = plot.name)
	best.alpha.index <- order(tot.res.ssq.per.alpha)[1]
	if (!is.na(ad.coeff.per.alpha[best.alpha.index])) { # an additive fit was part of the best fit, so plot this 
		lines(c(0,result$alphas[best.alpha.index]),c(result$ad.coeff.per.alpha[best.alpha.index],result$ad.coeff.per.alpha[best.alpha.index]))
		if (!is.na(mu.coeff.per.alpha[best.alpha.index])) {
			lines(c(result$alphas[best.alpha.index],result$alphas[length(alphas)]),
					c(result$ad.coeff.per.alpha[best.alpha.index],
							result$ad.coeff.per.alpha[best.alpha.index]+(result$alphas[length(alphas)]-result$alphas[best.alpha.index])*mu.coeff.per.alpha[best.alpha.index]))
		} else {  # if only an additive fit was done, then stretch the line to the end
			lines(c(0,result$alphas[length(alphas)]),c(result$ad.coeff.per.alpha[best.alpha.index],result$ad.coeff.per.alpha[best.alpha.index]))
		}
	} else {
		if (!is.na(mu.coeff.per.alpha[best.alpha.index])) {
			lines(c(result$alphas[best.alpha.index],result$alphas[length(alphas)]),c(0,
							0+(result$alphas[length(alphas)]-result$alphas[best.alpha.index])*mu.coeff.per.alpha[best.alpha.index]))
		}  
	}
	
	result$best.fit <- matrix(c(result$alphas[best.alpha.index], result$ad.coeff.per.alpha[best.alpha.index], result$mu.coeff.per.alpha[best.alpha.index], result$tot.res.ssq.per.alpha[best.alpha.index]))
	rownames(result$best.fit) <- c("alpha","additiveCoefficient","multiplicativeCoefficient","ssQresiduals")
	
	# prune results and change order of results to most significant first
	newResult <- list()
	if (!is.null(ordered.bins.by.mean$name.of.value.col) && !is.na(ordered.bins.by.mean$name.of.value.col)) { newResult$name.of.value.col <- ordered.bins.by.mean$name.of.value.col } 
	newResult$best.fit <- result$best.fit
	newResult$alphas <- result$alphas
	newResult$tot.res.ssq.per.alpha <- result$tot.res.ssq.per.alpha
	newResult$ad.coeff.per.alpha <- result$ad.coeff.per.alpha
	newResult$mu.coeff.per.alpha <- result$mu.coeff.per.alpha
	newResult
}
