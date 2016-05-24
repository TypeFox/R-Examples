pno <- function(path_bioclim, path_model, subset = NULL, 
                bin_width = 1, 	bin_number = NULL){
	
	# read and reorganize BIOCLIM layer: slow
	# ---------------------------------------
	hdr <- read.table(path_bioclim, row.names = 1, nrows = 6)
	bc.cells.dim <- hdr[grep("nrows|ncols", rownames(hdr)), ]
	nodata <- hdr[rownames(hdr) == "NODATA_value", ]
	bc <- scan(path_bioclim, what = "i", skip = 6, quiet = TRUE)
	bc <- as.numeric(bc)
	bc.na <- which(bc == nodata)
	bc[bc.na] <- NA
	
	# create bins for climate values:
	# -------------------------------
	bc_min <- min(bc, na.rm = TRUE)
	bc_max <- max(bc, na.rm = TRUE)
	if (is.null(bin_number)) {
		by <- bin_width
	}		
	else {
		by <- ((bc_max - bc_min)/(bin_number - 1))
	}
	cats  <- seq(from = bc_min, to = bc_max, by = by)
	cat("\nBioclimatic data ranges from", bc_min, "to", bc_max, 		"and will be binned into", length(cats), 		
	    "categories with a bin width of", by, "\n")
	
	# get species names
	specs <- list.files(path = path_model, pattern = ".asc",
	    full.names = TRUE)
	if (length(specs) == 0) 
	    stop("\nNo models found. Please check your path_model ", 
	        "argument!")
	# clamping maps:
	clamping <- grep("clamping", specs)
	if (length(clamping) > 0)
		specs <- specs[-clamping]
	if (!is.null(subset))
		specs <- specs[grep(paste(subset, collapse = "|"),
		    specs)]

	# get taxon names
	label <- gsub(path_model, "", specs)
	label <- gsub("/|.asc", "", label)
		
	# initialize matrix to store results
	# ----------------------------------
	x <- matrix(nrow = length(cats), ncol = length(specs) + 1)
	colnames(x) <- c("variable", label)
	x[, 1] <- cats
	
	# DEFINE: calculate cumulative probabilities for each value
	# on the environmental gradient
	# -----------------------------
	cum.prob.bin <- function(x, bc, enm, normalizer, bin){
	    id <- which(bc > x - bin/2 & bc <= x + bin/2)			
		sum(enm[id], na.rm = TRUE) / normalizer
	}
	
	# loop over MAXENT models (MM)
	# ----------------------------
	for (h in seq(along = specs)) {
		
		# reading data
		# ------------
		cat("\n\nReading ENM for", specs[h], "...")
		hdr <- read.table(specs[h], row.names = 1, nrows = 6)
		enm.cells.dim <- hdr[grep("nrows|ncols", 
		    rownames(hdr)), ]
	    nodata <- hdr[rownames(hdr) == "NODATA_value", ]
	    enm <- scan(specs[h], what = "i", skip = 6, quiet = TRUE)
	    enm <- as.numeric(enm)
	    enm.na <- which(enm == nodata)
	    enm[enm.na] <- NA
		cat("finished.")
		
		if (!identical(bc.cells.dim, enm.cells.dim))
			stop("Resolution and/or extent of ENM and", 				            " bioclimatic layer do not match!")
			
		# check differing coastlines
		# --------------------------
	    px <- abs(length(bc.na) - length(enm.na)) 
	    if (px > 0)
	        cat("\n\nWARNING: Raster maps differ by", 
	            px, "cells",
			    "(perhaps due to differing coastlines.)",
		        "\nThese cells are treated as having zero", 
		        "probability in the MAXENT distribution.")
		    
	    normalizer <- sum(enm, na.rm = TRUE)
		
	    cat("\n\nBinning probabilities for", specs[h], "...")
	    prob <- sapply(cats, cum.prob.bin, bc = bc, 
	        enm = enm, normalizer = normalizer, bin = by)
	    cat("finished.")
		
	    # remove ENM
	    # ----------------------------
	    remove(enm)
		
	    # fill data into matrix
	    # ---------------------
	    x[, h + 1] <- prob
			
	} # end of loop over MAXENT models
	    
    # original scale of Temperature values:
	# -------------------------------------
	if (length(grep("Temperature|Diurnal", path_bioclim)) == 1)
		x[, 1] <- x[, 1] / 10
	if (length(grep("Temperature_Seasonality|Isothermality",
	    path_bioclim)) == 1)
		x[, 1] <- x[, 1] / 100
	x
}
