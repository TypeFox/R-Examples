
read.DICOM.RT <- function(path, exclude=NULL, recursive=TRUE, verbose=TRUE, limit=NULL, DVH=TRUE, zDVH=FALSE, ...) {
	if (length(list.files(path)) == 0 && file.exists(path)) {
    	filenames <- path
	}
	else {
      filenames <- list.files(path, full.names=TRUE, recursive=recursive)
    }
	if (! is.null(exclude)) {
    	filenames <- grep(exclude, filenames, ignore.case=TRUE, value=TRUE, invert=TRUE)
  	}
  	if (length(filenames) < 1) {
  		warning("No files to read from path '", path, "'", sep="")
  		return()
  	}
	if (verbose) {
		cat("Reading ", length(filenames), " DICOM files from path: '", path, "' ... ", sep="")
	}
	DICOMs <- readDICOM(path, verbose=FALSE, exclude=exclude, recursive=recursive, ...)
	
	if (verbose) {
		cat("FINISHED\nExtracting CT data ... ", sep="")
	}
	modalities <- as.character(unlist(lapply(DICOMs$hdr, function(x) {x[which(x[,"name"]=="Modality"), "value"]})))

#############################
## IMPORT CT IMAGE FILE(S) ##
#############################	

	CT <- as.numeric(which(modalities == "CT"))
	frame.ref.CT <- as.character(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "FrameOfReferenceUID"), "value"])
	## ASSUMES CONSTANT VOXEL SIZE AND SLICE THICKNESS FOR ALL DICOM FILES IN CT!!!!
	voxel.size <- as.numeric(unlist(strsplit(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "PixelSpacing"), "value"], " ")))
	image.position <- as.numeric(unlist(strsplit(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "ImagePositionPatient"), "value"], " ")))
	z.slices <- unlist(lapply(DICOMs$hdr[CT], function(x) { as.numeric(unlist(strsplit(x[which(x[,"name"] == "ImagePositionPatient"), "value"], " "))[3]) }))
	voxel.size <- c(voxel.size, as.numeric(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "SliceThickness"), "value"]))
	patient.name <- as.character(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "PatientsName"), "value"])
	patient.ID <- as.character(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "PatientID"), "value"])
	slope <- as.numeric(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "RescaleSlope"), "value"])
	intercept <- as.numeric(DICOMs$hdr[[CT[1]]][which(DICOMs$hdr[[CT[1]]][,"name"] == "RescaleIntercept"), "value"])
	CT <- create3D(list(hdr=DICOMs$hdr[CT], img=DICOMs$img[CT]))
	CT <- CT*slope + intercept
	dimnames(CT) <- list((1:dim(CT)[1]-1)*voxel.size[1]+image.position[1], (1:dim(CT)[2]-1)*voxel.size[2]+image.position[2], z.slices) 
	if (verbose) {
		cat("FINISHED [", length(z.slices), " slices, ", paste(sprintf("%.*f", 1, voxel.size), collapse="x", sep=""), "mm res]\n", sep="")
	}
	
	if (length(which(modalities == "RTPLAN")) < 1) {
		dose.rx <- NA		
	}
#################################
## IMPORT RT PLAN PARAMETER(S) ##
#################################	

	for (i in as.numeric(which(modalities == "RTPLAN"))) {
		if (verbose) {
			cat("Reading RT plan information from file: '", filenames[i], "' ... ", sep="")
		}
		if (inherits(try(DICOM.i <- readDICOMFile(filenames[i], skipSequence=FALSE), silent=TRUE), "try-error")) {
			if (verbose) {
				warning("Unable to read DICOM file: ", filenames[i])
				cat("ERROR\n")
			}
			next
		}
		N.fractions <- as.numeric(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "NumberOfFractionsPlanned"), "value"])
		structureset.exists <- toupper(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "RTPlanGeometry"), "value"]) == "PATIENT"
		if (structureset.exists) {
			structureset.ID <- as.character(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "ReferencedSOPInstanceUID"), "value"])
		}
		else {
			structureset.ID <- ""
		}
		dose.ref.type <- unique(toupper(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "DoseReferenceType"), "value"]))
		if ((length(which(DICOM.i$hdr[,"name"] == "DoseReferenceSequence")) == 1) & (length(dose.ref.type) == 1)) {
			switch(dose.ref.type,
				TARGET = {
					dose.rx <- max(as.numeric(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "TargetPrescriptionDose"), "value"]), na.rm=TRUE)			
				},
				ORGAN_AT_RISK = {
					dose.rx <- max(as.numeric(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "DeliveryMaximumDose"), "value"]), na.rm=TRUE)			
				},
				dose.rx <- NA
			)
		}
		else {
			dose.rx <- NA
		}
		if (verbose) {
			cat("FINISHED\n")
		}
	}
	
#########################
## IMPORT DOSE GRID(S) ##
#########################	

	doses <- c()
	for (i in as.numeric(which(modalities == "RTDOSE"))) {
		if (verbose) {
			cat("Reading dose grid from file: '", filenames[i], "' ... ", sep="")
		}
		if (length(doses) > 1) {
			warning("Multiple dose grid files specified (may only use one dose grid at a time)")
			cat("ERROR\n")
			break
		}
		if (inherits(try(DICOM.i <- readDICOMFile(filenames[i], skipSequence=FALSE), silent=TRUE), "try-error")) {
			if (verbose) {
				warning("Unable to read DICOM file: ", filenames[i])
				cat("ERROR\n")
			}
			next
		}
		DICOMs$hdr[[i]] <- DICOM.i$hdr
		frame.ref.CT.i <- as.character(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "FrameOfReferenceUID"), "value"])
		
		if (length(frame.ref.CT.i) < 1) {
			frame.ref.CT.i <- frame.ref.CT
			warning("No reference frame in dose grid file")
		}
		if (frame.ref.CT != frame.ref.CT.i) {
			if (verbose) {
				warning("Reference frame mismatch")
				cat("ERROR\n")
			}			
			next
		}
		pixel.size <- as.numeric(unlist(strsplit(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "PixelSpacing"), "value"], " ")))
		z.size <- mean(diff(as.numeric(unlist(strsplit(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "GridFrameOffsetVector"), "value"], " ")))))
		image.origin <- as.numeric(unlist(strsplit(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "ImagePositionPatient"), "value"], " ")))
		grid.scale <- as.numeric(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "DoseGridScaling"), "value"])
		dose.units <- toupper(DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "DoseUnits"), "value"])
		plan.type <- DICOM.i$hdr[which(DICOM.i$hdr[,"name"] == "DoseSummationType"), "value"]
		# Type of Dose Summation... Defined Terms: PLAN = dose calculated for entire RT Plan MULTI_PLAN = dose calculated for 2 or more RT Plans FRACTION = dose calculated for a single Fraction Group within RT Plan BEAM = dose calculated for one or more Beams within RT Plan BRACHY = dose calculated for one or more Brachy Application Setups within RT Plan CONTROL_POINT = dose calculated for one or more Control Points within a Beam
		DICOMs$img[[i]] <- DICOM.i$img * grid.scale
		temp <- array(data=NA, dim=dim(DICOM.i$img)[c(2:1,3)])
		for (j in 1:dim(DICOM.i$img)[3]) {
			temp[,,j] <- t(DICOMs$img[[i]][,,j])
		}
		dimnames(temp) <- list(
			(1:dim(temp)[1]-1)*pixel.size[1]+image.origin[1],
			rev(2*mean(range(as.numeric(dimnames(CT)[[2]]))) - (((1:dim(temp)[2])-1)*pixel.size[2]+image.origin[2])),
			(1:dim(temp)[3]-1)*z.size+image.origin[3]
		)
		DICOMs$img[[i]] <- temp
		if (length(unique(dose.units)) != 1) {
			warning("Disagreement of dose specification within dose grid file")
			if (verbose) {
				cat("ERROR\n")
			}
			doses <- c(doses, NA)
			next
		}
		dose.units <- dose.units[1]
		switch(dose.units,
			GY = dose.units <- "Gy",
			CGY = dose.units <- "cGy",
			{
				warning("Dose not specified as 'Gy' or 'cGy'")
				if (verbose) {
					cat("ERROR\n")
				}
				doses <- c(doses, NA)
				next
			}
		)
		if (verbose) {
			cat("[Dose type = ", plan.type, ", Units = ", dose.units, "] ", sep="")
		}
		attr(DICOMs$img[[i]], "dose.units") <- dose.units
		doses <- c(doses, list(DICOMs$img[[i]]))
#		doses.hdr <- DICOM.i$hdr
		if (verbose) {
			cat("FINISHED\n")
		}
		## EXTRACT DVH DATA
#		if (FALSE) { ### THE FOLLOWING COMMANDS BYPASSED DUE TO UNRELIABILITY OF DVHs STORED WITHIN DICOM-RT DATA
#		if (DVH) {
#			if (verbose) {
#				cat("Extracting existent DVHs ... ", sep="")
#			}
#			dvh.start <- which(DICOM.i$hdr[,"name"] == "DVHType")
#			dvh.end <- which(DICOM.i$hdr[,"name"] == "DVHMeanDose")
#			if (length(dvh.start) < 1) {
#				warning(paste("Dose file '", file, "' contained no recognizable DVH structure(s)", sep=""))
#				if (verbose) {
#					cat("ERROR\n")
#				}
#				next
#			}
#			else if (length(dvh.start) == 1) {
 #  			 	dvhs <- list(DICOM.i$hdr[dvh.start:dvh.end,])	
#			}
#			else {
#				dvhs <- mapply(function(start, end) list(DICOM.i$hdr[start:end,]), dvh.start, dvh.end)
#			}
#			DVH.list <- lapply(dvhs,
#				function(dvh) {
#					ID <- dvh[which(dvh[, "name"] == "ReferencedROINumber"), "value"]
#					type <- toupper(dvh[which(dvh[, "name"] == "DVHType"), "value"])
#					switch(type,
#						CUMULATIVE = type <- "cumulative",
#						DIFFERENTIAL = type <- "differential"
#					)
#					dose.units <- toupper(dvh[which(dvh[, "name"] == "DoseUnits"), "value"])
#					switch(dose.units,
#						GY = dose.units <- "Gy",
#						CGY = dose.units <- "cGy",
#						{
#							warning("Dose not specified as 'Gy' or 'cGy'")
#							return(new("DVH", patient=patient.name, ID=patient.ID, structure.name=ID))
#						}
#					)
#					vol.units <- toupper(dvh[which(dvh[, "name"] == "DVHVolumeUnits"), "value"])
#					switch(vol.units,
#						CM3 = vol.type <- "absolute",
#						PERCENT = vol.type <- "relative",
#						vol.type <- "relative"
#					)
#					dvh.length <- as.numeric(dvh[which(dvh[, "name"] == "DVHNumberOfBins"), "value"])
#					data <- as.numeric(unlist(strsplit(dvh[which(dvh[, "name"] == "DVHData"), "value"], " ")))
#					vols <- data[1:dvh.length*2]
#					scale <- as.numeric(dvh[which(dvh[, "name"] == "DVHDoseScaling"), "value"])
#					doses <- cumsum(data[1:dvh.length*2-1]*scale)
#					min <- dose.rx*as.numeric(dvh[which(dvh[, "name"] == "DVHMinimumDose"), "value"])/100
#					mean <- dose.rx*as.numeric(dvh[which(dvh[, "name"] == "DVHMeanDose"), "value"])/100
#					max <- dose.rx*as.numeric(dvh[which(dvh[, "name"] == "DVHMaximumDose"), "value"])/100
#					return(new("DVH", patient=patient.name, ID=patient.ID, structure.name=ID, type=type, dose.units=dose.units, volume.type=vol.type, dose.type="absolute",doses=doses,volumes=vols,dose.min=min,dose.mean=mean,dose.max=max,dose.fx=N.fractions,dose.rx=dose.rx))
#				}
#			)
#			DVH.list.names <- unlist(lapply(DVH.list, function(dvh) {return(dvh$structure.name)}))
#			if (verbose) {
#				cat("FINISHED\n")
#			}		
#		}
#		else {
#			DVH.list.names <- DVH.list <- NULL
#		}
	}
	if (length(which(modalities == "RTDOSE")) < 1) {
		warning("Unable to extract DVH data from DICOM-RT (no dose grid available)")
		DVH <- FALSE
	}

##################################
## IMPORT STRUCTURE SET FILE(S) ##
##################################	

	first <- TRUE
	data.old <- data <- list(set=NULL, name=NULL, points=NULL, DVH=NULL)
	use.dose.grid <- c()
	for (i in as.numeric(which(modalities == "RTSTRUCT"))) {
		if (verbose) {
			cat("Reading structure set from file: '", filenames[i], "' ... ", sep="")
		}
		DICOMs$hdr[[i]] <- DICOM.i <- readDICOMFile(filenames[i], skipSequence=FALSE)$hdr
			
		structures <- DICOM.i[which(DICOM.i[,"name"] %in% c("ROIName", "ROINumber")),]
		N <- dim(structures)[1]/2
		frame.ref.CT.i <- as.character(DICOM.i[which(DICOM.i[,"name"] == "FrameOfReferenceUID"), "value"])
		structureset.i.ID <- as.character(DICOM.i[which(DICOM.i[,"name"] == "SOPInstanceUID"), "value"])
#		use.dose.grid <- c(use.dose.grid, !structureset.i.ID %in% structureset.ID)
		use.dose.grid <- c(use.dose.grid, TRUE) ## CALCULATE DVH FROM DOSE GRID RATHER THAN USING EXISTENT DVHs IN DICOM DOSE FILE
		if (length(frame.ref.CT.i) < 1) {
			frame.ref.CT.i <- frame.ref.CT
			warning("No reference frame in structure set file")
		}
		if (frame.ref.CT != frame.ref.CT.i) {
			if (verbose) {
				warning("Reference frame mismatch")
				cat("ERROR\n")
			}			
			next
		}
		structureset <- as.character(DICOM.i[which(DICOM.i[,"name"] == "StructureSetName"), "value"])
		if (length(structureset) == 0) {
			structureset <- as.character(DICOM.i[which(DICOM.i[,"name"] == "StructureSetLabel"), "value"])
		}
		if (N < 1) {
			if (verbose) {
				warning("Empty structure set")
				cat("ERROR\n")
			}			
			next
		}
		if (verbose) {
			cat("(", N, " structures identified) ", sep="")
		}
		structure.IDs <- as.numeric(structures[1:N*2-1, "value"])
		names(structure.IDs) <- structures[1:N*2, "value"]
		colors <- as.numeric(which(DICOM.i[,"name"] %in% c("ROIDisplayColor")))
		col <- c()
		structures <- as.numeric(which(DICOM.i[,"name"] == "ReferencedROINumber"))
		contour.seq <- as.numeric(which(DICOM.i[,"name"] == "ContourSequence"))
		contours <- as.numeric(which(DICOM.i[,"name"] == "ContourData"))
		if (length(contour.seq) < 1) {
			warning(paste("Structure set from file '", filenames[i], "' is empty", sep=""))
			if (verbose) {
				cat("ERROR\n")
			}			
			next			
		}
		if (!first) {
			data.old$set <- c(data.old$set, data$set)
			data.old$name <- c(data.old$name, list(data$name))
			data.old$points <- c(data.old$points, list(data$points))
			if (DVH) {
				data.old$DVH <- c(data.old$DVH, list(data$DVH))
			}
		}
		else {
			data.old <- list(set=NULL, name=NULL, points=NULL, DVH=NULL)
			first <- FALSE
		}
		data <- list(set=structureset, name=names(structure.IDs), points=vector("list", N), DVH=vector("list", N))
		N.ROIs <- length(structures)
		used <- c()
		for (j in 1:length(contour.seq)) {
			structures.j <- structures[which(structures > contour.seq[j])]
			structure.j <- structures.j[which.min(structures.j)]
			data.j <- strsplit(DICOM.i[intersect(contour.seq[j]:structure.j, contours), "value"], " ")
			struct.ID.j <- which(structure.IDs == as.numeric(DICOM.i[structure.j, "value"]))
#			if (DVH) {
#				DVH.j <- which(DVH.list.names == as.numeric(DICOM.i[structure.j, "value"]))
#			}
			if (length(struct.ID.j) < 1) {
				warning(paste("Expected structure not matched in DICOM file", sep=""))
				next
			}
			used <- c(used, struct.ID.j)
			if (length(data.j) < 1) {
				warning(paste("Structure '", names(structure.IDs)[struct.ID.j], "' is empty", sep=""))
				data$points[[struct.ID.j]] <- NA
				if (DVH) {
#					if (length(DVH.j) < 1) {
						data$DVH[[struct.ID.j]] <- new("DVH", patient=patient.name, ID=patient.ID, structure.name=names(structure.IDs)[struct.ID.j])
#					}
#					else {
#						DVH.j <- DVH.list[[DVH.j]]
#						DVH.j$structure.name <- names(structure.IDs)[struct.ID.j]
#						data$DVH[[struct.ID.j]] <- DVH.j
#					}
					next
				}
			}
			data.j <- lapply(data.j,
				function(x) {
					x <- as.numeric(x)
					if (length(x) < 3) {
						warning(paste("Structure '", names(structure.IDs)[struct.ID.j], "' is missing slices", sep=""))
						return(NA)
					}
					x <- cbind(x[1:(length(x)/3)*3-2], x[1:(length(x)/3)*3-1], x[1:(length(x)/3)*3])
					x[,2] <- sum(range(as.numeric(dimnames(CT)[[2]]))) - x[,2]
					return(x)
				}
			)			
			data$points[[struct.ID.j]] <- data.j
			if (DVH) {
#				if (length(DVH.j) < 1) {
					data$DVH[[struct.ID.j]] <- new("DVH", patient=patient.name, ID=patient.ID, structure.name=names(structure.IDs)[struct.ID.j])
#				}
#				else {
#					DVH.j <- DVH.list[[DVH.j]]
#					DVH.j$structure.name <- names(structure.IDs)[struct.ID.j]
#					data$DVH[[struct.ID.j]] <- DVH.j
#				}
			}
		}
		if (length(setdiff(1:N, used)) > 0) {
			warning(paste("Structure(s) ", paste("'", names(structure.IDs)[setdiff(1:N, used)], "'", collapse=", ", sep=""), " are empty", sep=""))
			for (k in setdiff(1:N, used)) {
				data$points[[k]] <- NA
				if (DVH) {
					data$DVH[[k]] <- new("DVH", patient=patient.name, ID=patient.ID, structure.name=names(structure.IDs)[k])
				}
			}
		}
		if (verbose) {
			cat("FINISHED\n")
		}
	}
	data$set <- c(data.old$set, data$set)
	data$name <- c(data.old$name, list(data$name))
	data$points <- c(data.old$points, list(data$points))
	if (DVH) {
		data$DVH <- c(data.old$DVH, list(data$DVH))
	}

	if (length(unlist(data$name)) >= 1) {
		if (verbose) {
			cat("Processing (", length(unlist(data$name)), ") structures:\n", sep="")
		}		
	}
	else {
		warning("No structure set(s) available for import")
		if (length(doses) > 0) {
			return(new("RTdata", name=path, CT=CT, dose=doses[[1]]))
		}
		else {
			return(new("RTdata", name=path, CT=CT))
		}
	}
	N <- length(data$name)
	struct.list <- new("structure.list")
	if (is.null(limit)) {
		limit <- Inf
	}
	for (i in 1:N) {
		for (j in 1:length(data$name[[i]])) {
			struct.i <- data$points[[i]][[j]]
			if (length(unlist(struct.i, recursive=FALSE)) > limit) {
				if (verbose) {
					cat("  ", data$set[[i]], ": ", data$name[[i]][j], " [", length(struct.i), " axial slice(s), ", length(unlist(struct.i, recursive=FALSE))/3, " point(s)] ... skipped\n", sep="")
				}
				if (DVH) {
					struct.list <- c(struct.list, new("structure3D", name=paste(data$name[[i]][j], data$set[[i]]), DVH=data$DVH[[i]][[j]]))
				}
				else {
					struct.list <- c(struct.list, new("structure3D", name=paste(data$name[[i]][j], data$set[[i]])))
				}
				next
			}
			else if (length(unlist(struct.i, recursive=FALSE)) > 1) {
				if (verbose) {
					cat("  ", data$set[[i]], ": ", data$name[[i]][j], " [", length(struct.i), " axial slice(s), ", length(unlist(struct.i, recursive=FALSE)), " point(s)] ... ", sep="")
				}
				if (identical(struct.i, NA)) {
					pts.i <- matrix(nrow=0, ncol=3)
				}
				else if (is.null(dim(struct.i))) {
					pts.i <- matrix(NA, nrow=0, ncol=3)
					for (k in 1:length(struct.i)) {
						pts.i <- rbind(pts.i, struct.i[[k]])
					}
				}						
				struct.i <- new("structure3D", name=paste(data$name[[i]][j], data$set[[i]]), vertices=pts.i, closed.polys=struct.i)
				if (DVH & use.dose.grid[i]) {
					if (zDVH) {
						if (verbose) {
							cat("calculating zDVH from dose grid ... ")
						}
						dvh.i <- calculate.DVH(struct.i, doses[[1]], resolution.xyz=c(pmin(voxel.size[1:2]/4, pixel.size/8, apply(range(struct.i),2,diff)[1:2]/100, na.rm=TRUE), voxel.size[3]), method="axial", dose.units=dose.units)
					}
					else {
						if (verbose) {
							cat("calculating DVH from dose grid ... ")
						}
						dvh.i <- calculate.DVH(struct.i, doses[[1]], resolution.xyz=c(pmin(voxel.size[1:2]/4, pixel.size/8, apply(range(struct.i),2,diff)[1:2]/100, na.rm=TRUE), voxel.size[3]), method="ATC", dose.units=dose.units)
					}
					if (is.null(dvh.i)) {
						warning(paste("Unable to calculate DVH for structure '", data$name[[i]][j], "_", data$set[[i]], "'", sep=""))
						if (verbose) {
							cat("ERROR\n")
						}			
						struct.list <- c(struct.list, struct.i)
						next						
					}
					struct.i$DVH <- dvh.i
					struct.i$DVH$dose.rx <- as.numeric(dose.rx)
					struct.i$DVH$dose.fx <- N.fractions
				}
				else if (DVH) {
					struct.i$DVH <- data$DVH[[i]][[j]]
				}
				struct.list <- c(struct.list, struct.i)
				if (verbose) {
					cat("FINISHED\n")
				}
			}
			else if ((length(unlist(struct.i, recursive=FALSE)) == 1) & (!is.na(struct.i))) {
				if (verbose) {
					cat("  ", data$set[[i]], ": ", data$name[[i]][j], " [", length(struct.i), " axial slice(s), ", length(unlist(struct.i, recursive=FALSE))/3, " point(s)] ... ", sep="")
				}
				if (identical(struct.i, NA)) {
					pts.i <- matrix(nrow=0, ncol=3)
				}
				else if (is.null(dim(struct.i))) {
					pts.i <- matrix(NA, nrow=0, ncol=3)
					for (k in 1:length(struct.i)) {
						pts.i <- rbind(pts.i, struct.i[[k]])
					}
				}					
				struct.i <- new("structure3D", name=paste(data$name[[i]][j], data$set[[i]]), vertices=pts.i, closed.polys=struct.i)
				if (DVH & use.dose.grid[i]) {
					if (zDVH) {
						if (verbose) {
							cat("calculating zDVH from dose grid ... ")
						}
						dvh.i <- calculate.DVH(struct.i, doses[[1]], resolution.xyz=c(pmin(voxel.size[1:2]/4, pixel.size/8, apply(range(struct.i),2,diff)[1:2]/100, na.rm=TRUE), voxel.size[3]), method="axial", dose.units=dose.units)
					}
					else {
						if (verbose) {
							cat("calculating DVH from dose grid ... ")
						}
						dvh.i <- calculate.DVH(struct.i, doses[[1]], resolution.xyz=c(pmin(voxel.size[1:2]/4, pixel.size/8, apply(range(struct.i),2,diff)[1:2]/100, na.rm=TRUE), voxel.size[3]), method="ATC", dose.units=dose.units)
					}	
					if (is.null(dvh.i)) {
						warning(paste("Unable to calculate DVH for structure '", data$name[[i]][j], "_", data$set[[i]], "'", sep=""))
						if (verbose) {
							cat("ERROR\n")
						}			
						struct.list <- c(struct.list, struct.i)
						next						
					}
					struct.i$DVH <- dvh.i
					struct.i$DVH$dose.rx <- as.numeric(dose.rx)
					struct.i$DVH$dose.fx <- N.fractions
				}
				else if (DVH) {
					struct.i$DVH <- data$DVH[[i]][[j]]
				}
				struct.list <- c(struct.list, struct.i)
				if (verbose) {
					cat("FINISHED\n")
				}
			}
			else {
				if (verbose) {
					cat("  ", data$set[[i]], ": ", data$name[[i]][j], " [EMPTY] ... FINISHED\n", sep="")
				}
				if (DVH) {
					struct.list <- c(struct.list, new("structure3D", name=paste(data$name[[i]][j], data$set[[i]]), DVH=data$DVH[[i]][[j]]))				
				}
				else {
					struct.list <- c(struct.list, new("structure3D", name=paste(data$name[[i]][j], data$set[[i]])))				
				}
			}
		}
	}

	if (length(doses) > 0) {
		return(new("RTdata", name=path, CT=CT, dose=doses[[1]], structures=struct.list))
	}
	else {
		return(new("RTdata", name=path, CT=CT, structures=struct.list))		
	}
	## FOR OTHER FILES LOAD SPECIFIC DICOM files with skipSequence=FALSE and re-store hdr info in DICOM list!

}

