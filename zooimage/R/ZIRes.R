## Copyright (c) 2004-2012, Ph. Grosjean <phgrosjean@sciviews.org>
##
## This file is part of ZooImage
## 
## ZooImage is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 2 of the License, or
## (at your option) any later version.
## 
## ZooImage is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
## 
## You should have received a copy of the GNU General Public License
## along with ZooImage.  If not, see <http://www.gnu.org/licenses/>.

print.ZIRes <- function (x, ...)
{
	X <- x
	class(X) <- "data.frame"
	print(X)
	## Are there size spectra?
	spectrum <- attr(x, "spectrum")
	if (length(spectrum)) {
		cat("\nWith size spectrum:\n")
		print(spectrum)
	}
	invisible(x)
}

## TODO... with inspirations from histSpectrum() and plotAbdBio()
#plot.ZIRes <- function (x, y, ...)
#{
#	
#}

#histSpectrum <- function (spect, class = 1:18 * 0.3 / 3 + 0.2, lag = 0.25,
#log.scale = TRUE, width = 0.1, xlab = "classes (mm)",
#ylab = if (log.scale) "log(abundance + 1)/m^3" else "Abundance (ind./m^3",
#main = "", ylim = c(0, 2), plot.exp = FALSE)
#{	
#	## Plot of histograms and optionally line for exponential decrease
#	## for size spectra
#	plot.exp <- isTRUE(as.logical(plot.exp))
#	log.scale <- isTRUE(as.logical(log.scale))
#	if (plot.exp) {
#		spect.lm <- lm(spect ~ class)
#		print(summary(spect.lm))
#		slope <- format(coef(spect.lm)[2], digits = 3)
#		main <- paste(main, " (slope = ", slope, ")", sep = "")
#		class2 <- class - lag
#		spect.lm2 <- lm(spect ~ class2)
#		if (log.scale) {
#			spect <- 10^spect - 1
#			expdat <- 10^predict(spect.lm2) - 1
#		}
#	}
#	barplot(spect, width = 0.1, space = 0, xlab = xlab, ylab = ylab,
#		main = main, ylim = ylim)
#	if (plot.exp) {
#		if (log.scale) {
#			abline(coef = coef(spect.lm2), col = 2, lwd = 2)
#		} else {
#			lines(class2, expdat, col = 2, lwd = 2)
#		}
#		return(invisible(spect.lm2))
#	}
#}
#
#plotAbdBio <- function (t, y1, y2, y3, ylim = c(0,3), xlab = "Date",
#ylab = "log(abundance + 1)", main = "", cols = c("green", "blue", "red"),
#pchs = 1:3, hgrid = 1:3, vgrid = t, vline = NULL, xleg = min(vgrid),
#yleg = ylim[2], legend = c("series 1", "series 2", "series 3"), type = "o")
#{	
#	## Custom plot for abundance and biomass
#	plot(t, y1, type = type, ylim = ylim, xlim = range(vgrid), ylab = ylab,
#		xlab = xlab, main = main, col = cols[1], xaxt = "n", pch = pchs[1])
#	axis(1, at = vgrid, labels = format(vgrid, "%b"))
#	lines(t, y2, type = type, col = cols[2], pch = pchs[2])
#	lines(t, y3, type = type, col = cols[3], pch = pchs[3])
#	
#	## Grid
#	abline(h = hgrid, col = "gray", lty = 2)
#	abline(v = vgrid, col = "gray", lty = 2)
#	
#	## Vertical line(s) to spot particular time events
#	if (!is.null(vline))
#		abline(v = as.Date(vline), lty = 2, lwd = 2, col = 2)
#	if (!is.null(xleg))
#		legend(xleg, yleg, legend, col = cols, lwd = 1, pch = pchs,
#			bg = "white")
#}

rbind.ZIRes <- function (..., deparse.level = 1)
{
	## Same as rbind.data.frame, but take care also to combine spectrum attributes
	res <- rbind.data.frame(..., deparse.level = deparse.level)
	
	attr(res, "spectrum") <- do.call("c", lapply(list(...), attr,
		which = "spectrum"))
	res
}

## Calculate abundances, biomasses and size spectra per class in a sample
processSample <- function (x, sample, keep = NULL, detail = NULL, classes = "both",
header = c("Abd", "Bio"), biomass = NULL, breaks = NULL)
{
	## Check arguments
	if (missing(sample)) {
		sample <- unique(sampleInfo(x$Label, type = "sample", ext = ""))
		if (length(sample) != 1) {
			warning("'sample' not provided, or 'x' does not contain a single sample 'Label'")
			return(NULL)
		}
	}
	if (!is.character(sample) || length(sample) != 1)
		stop("'sample' must be a single character string")
	
	header <- as.character(header)
	if (length(biomass)) {
		if (length(header) < 2) 
			stop("you must provide headers for abundances and biomasses")
		header <- header[1:2]
	} else {
		if (length(header) < 1)
			stop("You must provide a header for abundances")
		header <- header[1]
	}
	
	if (!length(x$Dil) || !is.numeric(x$Dil)) {
		warning("'Dil' column is missing or not numeric in 'x'")
		return(NULL)
	}
	
	## Extract only data for a given sample
	allSamples <- unique(sampleInfo(x$Label, type = "sample", ext = ""))
	if (!sample %in% allSamples){
		warning("Sample not found in 'x'")
		return(NULL)
	}
	x <- x[allSamples == sample, ]
	
	## Retrieving classes
	classes <- as.character(classes)[1]
	Cl <- switch(classes,
		Class = x$Class,
		Predicted = x$Predicted,
		both = { # Use Class where it is defined, otherwise, use predicted
			res <- x$Class
			if (is.null(res)) x$Predicted else {
				isMissing <- is.na(res)
				res[isMissing] <- x$Predicted[isMissing]
				res
			}
		},
		x[, classes])
	if (length(Cl) != NROW(x) || !is.factor(Cl)) { # There is a problem retrieving classes!
		warning("problem while retrieving classes (are they defined?)")
		return(NULL)
	} else x$Cl <- Cl
	
	## Subsample, depending on which classes we keep
	if (length(keep)) {
		keep <- as.character(keep)
		if (!all(keep %in% levels(x$Cl))) {
			warning("one or more 'keep' levels are not found in the classes")
			return(NULL)
		}
		x <- x[x$Cl %in% keep, ] # Select keep levels
	}
	if (NROW(x) == 0) {
		warning("no data left for this sample in 'x' when 'keep' is applied")
		return(NULL)
	}
		
	## Data for biomass calculation
	if (length(biomass)) {
		if (inherits(biomass, "data.frame")) { # Parameters vary by class
			## We need Class, P1, P2 and P3, and among groups, we need [other] 
			if (NCOL(biomass) != 4 && !is.factor(biomass[, 1]))
				stop("you must provide a data frame with four columns and first one as factor")
			if (!"[other]" %in% levels(biomass[, 1]))
				stop("you must include '[other]' in the levels of factor for biomass conversion")
			## Make sure the three other columns are numeric
			biomass[, 2] <- as.numeric(biomass[, 2])
			biomass[, 3] <- as.numeric(biomass[, 3])
			biomass[, 4] <- as.numeric(biomass[, 4])
			## Place P1, P2 and P3 according to class in x
			isother <- biomass[, 1] == "[other]"
			defbio <- biomass[isother, 2:4]
			nms <- as.character(biomass[!isother, 1])
			P1bio <- structure(biomass[!isother, 2], names = nms)
			P1 <- P1bio[Cl]
			P1[is.na(P1)] <- defbio[1]
			x$P1 <- as.numeric(P1)
			P2bio <- structure(biomass[!isother, 3], names = nms)
			P2 <- P2bio[Cl]
			P2[is.na(P2)] <- defbio[2]
			x$P2 <- as.numeric(P2)
			P3bio <- structure(biomass[!isother, 4], names = nms)
			P3 <- P3bio[Cl]
			P3[is.na(P3)] <- defbio[3]
			x$P3 <- as.numeric(P3)
			
		} else if (length(biomass) == 3 && is.numeric(biomass)) { # Same parameters for all classes
			x$P1 <- biomass[1]
			x$P2 <- biomass[2]
			x$P3 <- biomass[3]
		} else stop("wrong 'biomass', must be NULL, a vector of 3 values or a data frame with Class, P1, P2 and P3")
		if(!is.numeric(x$ECD)) stop("'ECD' required for biomasses")
		x$BioWeight <- (x$P1 * x$ECD^x$P3 + x$P2) * x$Dil
	}
	
	## Split among detail, if provided
	Cl <- as.character(x$Cl)
	if (length(detail)) {
		# We want more details for one ore more groups...
		detail <- as.character(detail)
		## 'total' and 'others' calculated differently!
		detail <- detail[detail != "[total]" & detail != "[other]"]
		
		Cl[!Cl %in% detail] <- "[other]"
		x$Cl <- Cl
		res <- tapply(x$Dil, Cl, sum, na.rm = TRUE)
		res <- res[c(detail, "[other]")]
		res <- c(res, '[total]' = sum(x$Dil, na.rm = TRUE))
		names(res) <- paste(header[1], names(res))
		
		if (!missing(biomass)) {
			resbio <- tapply(x$BioWeight, Cl, sum, na.rm = TRUE)
			resbio <- resbio[c(detail, "[other]")]
			resbio <- c(resbio, '[total]' = sum(x$BioWeight, na.rm = TRUE))
			names(resbio) <- paste(header[2], names(resbio))
			res <- c(res, resbio)
		}
		
	} else { # Total abundance (and biomass) only
		res <- sum(x$Dil, na.rm = TRUE)
		if (!missing(biomass))
			res <- c(res, sum(x$BioWeight, na.rm = TRUE))
		names(res) <- paste(header, "[total]")
	}
	
	## Make the result a data frame with first column being Id, and make it
	## a ZIRes object inheriting from data frame
	res <- structure(data.frame(Id = sample, t(res), check.names = FALSE),
		class = c("ZI3Res", "ZIRes", "data.frame"))
	
	## Do we calculate size spectra?
	if (length(breaks)) {
		if (!is.numeric(breaks) || length(breaks) < 2)
			stop("'breaks' must be a vector of two or more numerics or NULL")
		
		if(!is.numeric(x$ECD)) {
			warning("'ECD' required for size spectra")
			return(NULL)
		}
		
		## For each image, calculate size spectra per classes
		tcut <- function (items, data, breaks) {
			data <- data[items, ]
			x <- data$ECD 
			## Cut by class
			res <- tapply(x, data$Cl, function (x, breaks)
				table(cut(x, breaks = breaks)), breaks = breaks)
			## For empty classes, make sure to get zero
			res <- lapply(res, function (x, breaks)
				if (is.null(x)) table(cut(-1, breaks = breaks)) else x,
				breaks = breaks)
			## Turn this into a matrix and multiply by Dil for this image
			do.call("rbind", res) * data$Dil[1]
		}
		
		## Get abundance breaks (ind/vol) by image and by class
		## and sum over all images
		if (length(detail)) {
			x$Cl <- factor(x$Cl, levels = c(detail, "[other]"))
		} else x$Cl <- as.factor(x$Cl)
		spectrum <- Reduce("+", tapply(1:NROW(x), x$Label, tcut,
			data = x, breaks = breaks))
		
		## Place [other] at the end and add [total]
		isother <- rownames(spectrum) == "[other]"
		if (any(isother)) {
			spectrum <- rbind(spectrum[!isother, , drop = FALSE],
				spectrum[isother, , drop = FALSE],
				'[total]' = apply(spectrum, 2, sum))
		} else {
			spectrum <- rbind(spectrum,
			'[total]' = apply(spectrum, 2, sum))	
		}
		
		## Eliminate all ine except total if detail is not provided
		if (!length(detail))
			spectrum <- spectrum[NROW(spectrum), , drop = FALSE]

		## Put this in a 'spectrum' attribute (named list)
		spectrum <- list(spectrum)
		names(spectrum) <- sample
		attr(res, "spectrum") <- spectrum
	}
	
	res
}

processSampleAll <- function (path = ".", zidbfiles, ZIClass = NULL, keep = NULL,
detail = NULL, classes = "both", header = c("Abd", "Bio"), biomass = NULL,
breaks = NULL)
{
	## First, switch to that directory                                       
	if (!checkDirExists(path)) return(invisible(FALSE))
	initdir <- setwd(path)
	on.exit(setwd(initdir))
	path <- "."	# Indicate we are now in the right path
	
	## Get the list of ZID(B) files to process
	if (missing(zidbfiles) || !length(zidbfiles)) {	# Compute them from path
		zidbfiles <- zidbList(".")
		## If no .zidb files, try .zid files instead
		if (!length(zidbfiles)) zidbfiles <- zidList(".")
	}
	
	## If there is no files to process, exit now
	if (!length(zidbfiles)) {
		warning("there are no ZID(B) files to process in ", getwd())
		return(invisible(FALSE))
	}
			
	## Process samples in the .zidb files
	message("Processing sample statistics for ZIDB files...")
	flush.console()
	nfiles <- length(zidbfiles)
	res <- NULL
	for (i in 1:nfiles) {
		progress(i, nfiles)
		zidbfile <- zidbfiles[i]
		if (hasExtension(zidbfile, "zidb")) {
			dat <- zidbDatRead(zidbfile)
		} else dat <- zidDatRead(zidbfile)
		## Do we predict the classes in the sample?
		if (length(ZIClass)) dat <- predict(ZIClass, dat, class.only = FALSE)

		## Process that one sample and merge with the rest
		res <- rbind(res, processSample(dat, keep = keep, detail = detail,
			classes = classes, header = header, biomass = biomass,
			breaks = breaks))
	}
	progress(101) # Clear progression indicator
	message(" -- Done! --")
	
	res
}
