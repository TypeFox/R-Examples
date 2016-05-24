################################################################################
# File: SpatialDesign-methods
# Purpose: Define S3 and S4 methods for class SpatialDesign
# Programmer: Tom Kincaid
# Date: June 5, 2015
################################################################################

summary.SpatialDesign <- function(object, ..., auxvar = NULL, spframe = NULL,
	tess_ind = TRUE, sbc_ind = FALSE, nrows = 5, dxdy = TRUE) {

# Ensure that object belongs to class SpatialDesign
	if(class(object) != "SpatialDesign")
			stop("\nThe object argument must be a member of class SpatialDesign.")

# Create the survey design summary
	dsum <- dsgnsum(object, auxvar = auxvar)

# Calculate spatial balance metrics for the survey design
	if(is.null(spframe)) {
		spbal <- NULL
	} else {
		if(!(class(spframe) %in% c("SpatialPointsDataFrame",
		 "SpatialLinesDataFrame", "SpatialPolygonsDataFrame")))
			stop("\nThe spframe argument must be a member of class SpatialPointsDataFrame, \nSpatialLinesDataFrame, or SpatialPolygonsDataFrame.")
		spbal <- spbalance(object, spframe, tess_ind, sbc_ind, nrows, dxdy)
	}

# Return the results list
	invisible(list("design summary" = dsum, "spatial balance statistics" = spbal))
}
setMethod("summary", signature(object = "SpatialDesign"),
	summary.SpatialDesign)

plot.SpatialDesign <- function(x, y, ..., spframe = NULL, stratum = NULL,
	mdcaty = NULL, auxvar = NULL, pdffile = NULL, width = 8, height = 10) {

# Ensure that x belongs to class SpatialDesign
	if(class(x) != "SpatialDesign")
		stop("\nThe x argument must be a member of class SpatialDesign.")

# Assign the sites (attributes) data frame and the design list from the
# SpatialDesign object

	sites <- x@data
	design <- x@design

# Determine whether an appropriate survey design frame object was supplied

	if(is.null(spframe))
		stop("\nAn object containing the survey design frame must be supplied as the spframe \nargument.")
	if(!(class(spframe) %in% c("SpatialPointsDataFrame", "SpatialLinesDataFrame", "SpatialPolygonsDataFrame")))
		stop("\nThe spframe argument must be a member of class SpatialPointsDataFrame, \nSpatialLinesDataFrame, or SpatialPolygonsDataFrame.")

# Assign strata names and number of strata from the design list

	snames <- names(design)
	nstrata <- length(snames)

# If stratum equals NULL, ensure that the design list specifies a single stratum
# Otherwise, ensure that the name provided for stratum identifies a column in
# the attributes data frame of the spframe object

	if(is.null(stratum)) {
		if(length(snames) > 1)
			stop("\nThe column from the attributes data frame of the spframe object that identifies \nstratum membership was not provided and the design list specifies more than one \nstratum.")
	} else {
		temp <- match(stratum, names(spframe@data), nomatch=0)
		if(temp == 0)
			stop(paste("\nThe value provided for the column from the attributes data frame of the \nspframe object that identifies stratum membership for each element in the \nframe, \"", stratum, "\", does not occur among the columns in the data slot \nof the spframe object.", sep=""))
}

# Ensure that the stratum variable in the sites data frame is a factor

	if(!is.factor(sites$stratum))
		sites$stratum <- as.factor(sites$stratum)


# Ensure that unequal probability category is a factor, and assign category
# names and number of categories from the attributes data frame

	sites$mdcaty <- factor(sites$mdcaty)
	cnames <- levels(sites$mdcaty)
	ncaty <- length(cnames)

# Determine the selection type for each strata from the design list

	seltype <- character(nstrata)
	for(i in 1:nstrata) seltype[i] <- design[[i]]$seltype

# If mdcaty equals NULL, ensure that the attributes data frame of the SpatialDesign object does not include unequal probability categories
# Otherwise, ensure that the name provided for mdcaty identifies a column in
# the attributes data frame of the spframe object

	if(is.null(mdcaty)) {
		if(length(cnames) > 1)
			stop("\nThe column from the attributes data frame of the spframe object that identifies \nunequal probability categories was not provided and the attributes data frame \nof the SpatialDesign object includes a column containing unequal probability \ncategories.")
	} else {
		temp <- match(mdcaty, names(spframe@data), nomatch=0)
		if(temp == 0)
			stop(paste("\nThe value provided for the column from the attributes data frame of the spframe \nobject that identifies unequal probability categories for each element in the \nframe, \"", mdcaty, "\", does not occur among the columns in the data slot \nof the SpatialDesign object.", sep=""))
	}

# Ensure that names provided by the auxvar argument occur in the attributes
# table of the SpatialDesign object

	if(!is.null(auxvar)) {
		temp <- match(auxvar, names(sites), nomatch=0)
		if(any(temp == 0)) {
			temp.str <- vecprint(auxvar[temp == 0])
			stop(paste("\nThe following values in the vector of auxiliary variable names do not occur \namong the columns in the data slot of the SpatialDesign object:\n", temp.str, sep=""))
		}
	}

# If requested, open the PDF file

	if(!is.null(pdffile))
		pdf(file = pdffile, width = width, height = height)

#
# This section handles finite survey designs
# 

	if(class(spframe) == "SpatialPointsDataFrame") {

# Plot the complete set of sample sites

		plot(spframe@coords, pch=20, xlab="x-coordinate", ylab="y-coordinate",
			main="Sample Sites", ...)
		points(sites$xcoord, sites$ycoord, pch=20, col="red")

# Plot sample sites by stratum

		if(nstrata > 1) {
			for(i in 1:nstrata) {
				ind <- spframe@data[, stratum] == snames[i]
				plot(spframe@coords[ind,], pch=20, xlab="x-coordinate",
					ylab="y-coordinate", main=paste("Sample Sites for Stratum:", snames[i]), ...)
				ind <- sites$stratum == snames[i]
				points(sites$xcoord[ind], sites$ycoord[ind], pch=20, col="red")
			}
		}

# Plot sample sites by design category

		if(all(seltype != "Continuous") & ncaty > 1) {
			for(i in 1:ncaty) {
				ind <- spframe@data[, mdcaty] == cnames[i]
				plot(spframe@coords[ind,], pch=20, xlab="x-coordinate",
					ylab="y-coordinate", main=paste("Sample Sites for Design Category:", cnames[i]), ...)
				ind <- sites$mdcaty == cnames[i]
				points(sites$xcoord[ind], sites$ycoord[ind], pch=20, col="red")
			}
		}

# Plot sample sites color-coded by design category for each stratum

		if(nstrata > 1 & all(seltype != "Continuous") & ncaty > 1) {
			cols <- rainbow(ncaty, s=0.75)
			for(i in 1:nstrata) {
				ind <- spframe@data[, stratum] == snames[i]
				plot(spframe@coords[ind,], pch=20, xlab="x-coordinate",
					ylab="y-coordinate", main=paste("Sample Sites Color-Coded by Design Category \nStratum:", snames[i]), ...)
				legend(x="bottomright", inset=0.05, legend=cnames, pch=20, col=cols)
				for(j in 1:ncaty) {
					ind <- sites$stratum == snames[i] & sites$mdcaty == cnames[j]
					points(sites$xcoord[ind], sites$ycoord[ind], pch=20, col=cols[j])
				}
			}
		}

# For each auxiliary variable, plot sample sites color-coded by category

		nvar <- length(auxvar)
		if(nvar > 0) {
			for(i in 1:nvar) {
				sites[,auxvar[i]] <- factor(sites[,auxvar[i]])
				cnames <- levels(sites[,auxvar[i]])
				ncaty <- length(cnames)
				cols <- rainbow(ncaty, s=0.75)
				plot(spframe@coords, pch=20, xlab="x-coordinate",
					ylab="y-coordinate", main=paste("Sample Sites Color-Coded by", auxvar[i], "Category"), ...)
				legend(x="bottomright", inset=0.05, legend=cnames, pch=20, col=cols)
				for(j in 1:ncaty) {
					ind <- sites[,auxvar[i]] == cnames[j]
					points(sites$xcoord[ind], sites$ycoord[ind], pch=20, col=cols[j])
				}
			}
		}

#
# This section handles linear and area survey designs
# 

	} else {

# Plot the complete set of sample sites

		plot(spframe, axes=TRUE, pch=20, xlab="x-coordinate", ylab="y-coordinate",
			main="Sample Sites", ...)
		points(sites$xcoord, sites$ycoord, pch=20, col="red")

# Plot sample sites by stratum

		if(nstrata > 1) {
			for(i in 1:nstrata) {
				ind <- spframe@data[, stratum] == snames[i]
				plot(spframe[ind,], axes=TRUE, pch=20, xlab="x-coordinate",
					ylab="y-coordinate", main=paste("Sample Sites for Stratum:", snames[i]), ...)
				ind <- sites$stratum == snames[i]
				points(sites$xcoord[ind], sites$ycoord[ind], pch=20, col="red")
			}
		}

# Plot sample sites by design category

		if(all(seltype != "Continuous") & ncaty > 1) {
			for(i in 1:ncaty) {
				ind <- spframe@data[, mdcaty] == cnames[i]
				plot(spframe[ind,], axes=TRUE, pch=20, xlab="x-coordinate",
					ylab="y-coordinate", main=paste("Sample Sites for Design Category:", cnames[i]), ...)
				ind <- sites$mdcaty == cnames[i]
				points(sites$xcoord[ind], sites$ycoord[ind], pch=20, col="red")
			}
		}

# Plot sample sites color-coded by design category for each stratum

		if(nstrata > 1 & all(seltype != "Continuous") & ncaty > 1) {
			cols <- rainbow(ncaty, s=0.75)
			for(i in 1:nstrata) {
				ind <- spframe@data[, stratum] == snames[i]
				plot(spframe[ind,], axes=TRUE, pch=20, xlab="x-coordinate",
					ylab="y-coordinate", main=paste("Sample Sites Color-Coded by Design Category \nStratum:", snames[i]), ...)
				legend(x="bottomright", inset=0.05, legend=cnames, pch=20, col=cols)
				for(j in 1:ncaty) {
					ind <- sites$stratum == snames[i] & sites$mdcaty == cnames[j]
					points(sites$xcoord[ind], sites$ycoord[ind], pch=20, col=cols[j])
				}
			}
		}

# For each auxiliary variable, plot sample sites color-coded by category

		nvar <- length(auxvar)
		if(nvar > 0) {
			for(i in 1:nvar) {
				sites[,auxvar[i]] <- factor(sites[,auxvar[i]])
				cnames <- levels(sites[,auxvar[i]])
				ncaty <- length(cnames)
				cols <- rainbow(ncaty, s=0.75)
				plot(spframe, axes=TRUE, pch=20, xlab="x-coordinate",
					ylab="y-coordinate", main=paste("Sample Sites Color-Coded by", auxvar[i], "Category"), ...)
				legend(x="bottomright", inset=0.05, legend=cnames, pch=20, col=cols)
				for(j in 1:ncaty) {
					ind <- sites[,auxvar[i]] == cnames[j]
					points(sites$xcoord[ind], sites$ycoord[ind], pch=20, col=cols[j])
				}
			}
		}
	}

# If necesssary, close the PDF file

	if(!is.null(pdffile))
		graphics.off()

}
setMethod("plot", signature(x = "SpatialDesign", y = "missing"),
	plot.SpatialDesign)
