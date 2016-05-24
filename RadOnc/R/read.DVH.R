read.DVH <- function (file, type=NA, verbose=TRUE, collapse=TRUE) {
	type <- match.arg(tolower(type), choices=c(NA, "aria10", "aria11", "aria8", "dicom", "cadplan", "tomo", "monaco", "raystation"), several.ok=TRUE)
	if (length(file) < 1) {
		warning("argument 'file' is missing, with no default")
		return()
	}
	
	read.DVH.file <- function (file, type, verbose=TRUE) {
		switch(type, 
			aria10 = return(read.DVH.Aria10(file=file, verbose=verbose)),
			aria8 = return(read.DVH.Aria8(file=file, verbose=verbose)),
			aria11 = return(read.DVH.Aria11(file=file, verbose=verbose)),
			dicom = return(read.DVH.DICOM(path=file, verbose=verbose)),
			cadplan = return(read.DVH.CadPlan(file=file, verbose=verbose)),
			tomo = return(read.DVH.TomoTherapy(file=file, verbose=verbose)),
			monaco = return(read.DVH.Monaco(file=file, verbose=verbose)),
			raystation = return(read.DVH.RayStation(file=file, verbose=verbose)),
			warning("DVH file format not specified for file '", file, "'")
		)
		return()
	}		

	
	if (length(file) == 1) {
		if (length(type) > 1) {
			warning("length of 'file' and 'type' do not match")
			type <- type[1]
		}
		return(read.DVH.file(file, type, verbose))
	}
	if (length(type) < length(file)) {
		if (length(type) > 1) {
			warning("length of 'file' and 'type' do not match")
		}
		type <- rep(type[1], length(file))
	}
	else if (length(type) > length(file)) {
		warning("length of 'file' and 'type' do not match")
		type <- type[1:length(file)]	
	}
	if (all(is.na(type))) {
		warning("'type' not specified")
		return()
	}
	DVH.list <- mapply(function(x, y, z) {list(read.DVH.file(x, y, z))}, file, type, verbose)
    names(DVH.list) <- sapply(DVH.list, function(x) {if (is.null(x)) {return("")} else {paste(x[[1]]@ID, x[[1]]@patient, sep="_")}})
    if (collapse) {
		return(as(DVH.list, "DVH.list"))    	
    }
    else {
		return(DVH.list)
    }
}


read.DVH.Aria10 <- function (file, verbose=TRUE) {
	if (!(fid <- file(file, open="r"))) {
		warning(paste("Could not open file '", file, "'", sep=""))		
		return()
	}
	if (verbose) {
		cat("Reading DVH file ('", file, "')... ", sep="")
	}
	data <- readLines(fid)
	close(fid)
	
    ## IDENTIFY STRUCTURES
    struct.start <- grep("^Structure: ", data, perl=TRUE)
    struct.end <- struct.start + diff(c(struct.start, length(data)+1)) - 1
    if (length(struct.start) < 1) {
		warning(paste("File '", file, "' contained no recognizable DVH structure(s)", sep=""))
		if (verbose) {
			cat("ERROR\n")
		}
		return()
    }
    else if (length(struct.start) == 1) {
   	 	structures <- list(data[struct.start:struct.end])
    }
    else {
	    structures <- mapply(function(start, end) list(data[start:end]), struct.start, struct.end)
	}
	# EXTRACT HEADER INFO
    header <- data[1:(struct.start[1]-1)]
    patient <- sub("^.*: (.*?)\\s*([(].*[)]|$).*", "\\1", header[grep("^Patient Name.*: ", header, ignore.case=TRUE, perl=TRUE)], perl=TRUE)
    ID <- sub("^.*: (.+$)", "\\1", header[grep("^Patient ID.*: ", header, ignore.case=TRUE, perl=TRUE)])
    plan.sum <- (grepl("Summed", header[grep("^Comment.*: ", header, ignore.case=TRUE, perl=TRUE)]))
    date <- sub("^.*: (.+$)", "\\1", header[grep("^Date.*: ", header, ignore.case=TRUE, perl=TRUE)])
	plan <- sub("^.*: (.+$)", "\\1", header[grep("^Plan.*: ", header, ignore.case=TRUE, perl=TRUE)])
	DVH.type <- header[grep("^Type.*: ", header, ignore.case=TRUE, perl=TRUE)]
	if (grepl("Cumulative", DVH.type, ignore.case=TRUE, perl=TRUE)) {
		DVH.type <- "cumulative"
	}
	else {
		DVH.type <- "differential"
	}
	# EXTRACT PRESCRIPTION DOSE AND DOSE UNITS
	dose.rx <- header[grep("^Prescribed dose.*: ", header, ignore.case=TRUE, perl=TRUE)]
    dose.units <- toupper(sub("^.*[[](.*)[]].*", "\\1", dose.rx, perl=TRUE))

	if (dose.units == "GY") {
		dose.units <- "Gy"
	}
	else if (dose.units == "CGY") {
		dose.units <- "cGy"
	}
	dose.rx <- suppressWarnings(as.numeric(sub("^Prescribed dose.*: ", "", dose.rx, ignore.case=TRUE, perl=TRUE)))
    rx.isodose <- suppressWarnings(as.numeric(sub(".*: ", "", header[grep("^[%] for dose.*: ", header, ignore.case=TRUE, perl=TRUE)])))
    
	if (verbose) {
		cat("[exported on ", date, "]\n", sep="")
		cat("  Patient: ", patient, " (", ID, ")\n", sep="")
		cat("  Plan: ", plan, "\n  Dose: ", if (is.na(dose.rx)) {"NOT SPECIFIED"} else {paste(dose.rx, dose.units, " (at ", rx.isodose, "% isodose line)", sep="")}, "\n", sep="")		
	}

	# EXTRACT DVH DATA FOR EACH STRUCTURE
	DVH.list <- lapply(structures,
		function (data) {
		    name <- sub("^.*: (.+$)", "\\1", data[grep("^Structure.*: ", data, ignore.case=TRUE, perl=TRUE)])
		    if (length(name) < 1) {
				warning("Invalid DVH file format, could not extract structure")
				return(new("DVH"))
		    }
		    volume <- suppressWarnings(as.numeric(sub("^.*: (.+$)", "\\1", data[grep("^Volume.*: ", data, ignore.case=TRUE, perl=TRUE)], perl=TRUE)))
		    if (length(volume) < 1) {
				warning("Invalid DVH file format, could not extract structure volume information")
				return(new("DVH"))
		    }
        	header <- grep("Dose [[](%|Gy|cGy)[]].*Volume", data, ignore.case=TRUE, perl=TRUE)
			if (grepl("^\\s*Dose [[](cGy|Gy)[]].*Volume", data[header], ignore.case=TRUE, perl=TRUE)) {
				dose.type <- "absolute"
			}
			else {
				dose.type <- "relative"
			}
			if (grepl(".*Volume [[][%][]]", data[header], ignore.case=TRUE, perl=TRUE)) {
				volume.type <- "relative"
			}
			else {
				volume.type <- "absolute"
			}
			getDose <- function(dose) {
				if (grepl("[[]%[]]", dose)) {
					if (dose.type == "absolute") {
						dose <- suppressWarnings(as.numeric(sub(".*: ", "", dose)))*dose.rx/rx.isodose		
					}
					else {
						dose <- suppressWarnings(as.numeric(sub(".*: ", "", dose)))
					}
				}
				else {
					if (dose.type == "absolute") {
						dose <- suppressWarnings(as.numeric(sub(".*: ", "", dose)))
					}
					else {
						dose <- suppressWarnings(as.numeric(sub(".*: ", "", dose)))*rx.isodose/dose.rx		
					}
				}
				return(dose)				
			}

			dose.min <- getDose(data[grep("^Min Dose.*: ", data, ignore.case=TRUE, perl=TRUE)])
			dose.max <- getDose(data[grep("^Max Dose.*: ", data, ignore.case=TRUE, perl=TRUE)])
			dose.mean <- max(0, getDose(data[grep("^Mean Dose.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			if (verbose) {
				cat("  ..Importing structure: ", name, "  [volume: ", volume, "cc, dose: ", dose.min, " - ", dose.max, dose.units, "]\n", sep="")
			}

			dose.mode <- max(0, getDose(data[grep("^Modal Dose.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.median <- max(0, getDose(data[grep("^Median Dose.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.STD <- max(0, getDose(data[grep("^STD.*: ", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)

		    equiv.sphere <- suppressWarnings(as.numeric(sub("^.*: (.+$)", "\\1", data[grep("^Equiv. Sphere Diam.*: ", data, ignore.case=TRUE, perl=TRUE)])))
		    conf.ind <- suppressWarnings(as.numeric(sub("^.*: (.+$)", "\\1", data[grep("^Conformity Index.*: ", data, ignore.case=TRUE, perl=TRUE)])))
		    gradient <- suppressWarnings(as.numeric(sub("^.*: (.+$)", "\\1", data[grep("^Gradient Measure.*: ", data, ignore.case=TRUE, perl=TRUE)])))

			con <- textConnection(data[(header+1):length(data)])
			dvh <- read.table(con, header=FALSE, stringsAsFactors=FALSE)
			close(con)
			data.dose <- dvh[, 1]
			if (plan.sum) {
				data <- dvh[, 2]
			}
			else {
				data <- dvh[, 3]
			}
			if (DVH.type == "differential") {
				data <- data * diff(c(-data.dose[1], data.dose))
				temp.doses <- data.dose - diff(c(-data.dose[1], data.dose))/2
				data.dose <- c(temp.doses, (2*data.dose - temp.doses)[length(temp.doses)])
				data <- diffinv(-data, xi=sum(data))
			}
			return(new("DVH", patient=patient, ID=ID, dose.min=dose.min, dose.max=dose.max, dose.mean=dose.mean, dose.mode=dose.mode, dose.median=dose.median, dose.STD=dose.STD, equiv.sphere=equiv.sphere, conf.index=conf.ind, gradient=gradient, plan.sum=plan.sum, dose.rx=dose.rx, rx.isodose=rx.isodose, structure.name=name, structure.volume=volume, doses=data.dose, volumes=data, type="cumulative", dose.type=dose.type, dose.units=dose.units, volume.type=volume.type))	
		}
	)
	
	# RETURN DVH LIST
	names(DVH.list) <- unlist(lapply(DVH.list, names))
	return(new("DVH.list", DVH.list))
}

read.DVH.Aria11 <- function (file, verbose=TRUE) {
	return(read.DVH.Aria10(file, verbose))
}


read.DVH.Aria8 <- function (file, verbose=TRUE) {
	warning("Aria 8 format not currently supported")
	return()
}

read.DVH.DICOM <- function(path, verbose=TRUE) {
	dicom <- read.DICOM.RT(path, verbose=verbose, DVH=TRUE)
	if (is.null(dicom)) return()
	return(as(dicom$structures, "DVH.list"))
}

read.DVH.CadPlan <- function(file, verbose=TRUE) {
	if (!(fid <- file(file, open="r"))) {
		warning(paste("Could not open file '", file, "'", sep=""))		
		return()
	}
	if (verbose) {
		cat("Reading DVH file ('", file, "')... ", sep="")
	}
	data <- readLines(fid)
	close(fid)
	
    ## IDENTIFY STRUCTURES
    struct.start <- grep("^Histogram.*:\\s*", data, perl=TRUE)
    struct.end <- struct.start + diff(c(struct.start, length(data)+1)) - 1
    if (length(struct.start) < 1) {
		warning(paste("File '", file, "' contained no recognizable DVH structure(s)", sep=""))
		if (verbose) {
			cat("ERROR\n")
		}
		return()
    }
    else if (length(struct.start) == 1) {
    	structures <- list(data[struct.start:struct.end])
    }
    else {
	    structures <- mapply(function(start, end) list(data[start:end]), struct.start, struct.end)
	}
	# EXTRACT HEADER INFO
    header <- data[1:(struct.start[1]-1)]
    patient <- sub("^.*: (.*)", "\\1", header[grep("^Patient Name.*:\\s*", header, ignore.case=TRUE, perl=TRUE)], perl=TRUE)
    ID <- sub("^.*: (.+$)", "\\1", header[grep("^Patient ID.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
    date <- sub("^.*: (.+$)", "\\1", header[grep("^Date.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
	plan <- sub("^PLAN\\s*(.+$)", "\\1", header[grep("^PLAN\\s+", header, ignore.case=TRUE, perl=TRUE)])
	DVH.type <- header[grep("Dose Volume Histogram", header, ignore.case=TRUE, perl=TRUE)]
	if (grepl("Cumulative", DVH.type, ignore.case=TRUE, perl=TRUE)) {
		DVH.type <- "cumulative"
	}
	else {
		DVH.type <- "differential"
	}

	if (verbose) {
		cat("[exported on ", date, "]\n", sep="")
		cat("  Patient: ", patient, " (", ID, ")\n", sep="")
		cat("  Plan: ", plan, "\n", sep="")		
	}
	
	# EXTRACT DVH DATA FOR EACH STRUCTURE
	DVH.list <- lapply(structures,
		function (data) {
			# EXTRACT PRESCRIPTION DOSE AND DOSE UNITS
			dose.rx <- data[grep("^Prescr[.] dose.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]
			rx.isodose <- data[grep("^[%] for dose.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]

    		dose.units <- toupper(sub("^.*[(](.*)[)].*", "\\1", dose.rx, perl=TRUE))
			if (dose.units == "GY") {
				dose.units <- "Gy"
			}
			else if (dose.units == "CGY") {
				dose.units <- "cGy"
			}
			dose.rx <- suppressWarnings(as.numeric(sub("^Prescr[.] dose.*:\\s*", "", dose.rx, ignore.case=TRUE, perl=TRUE)))
		    rx.isodose <- suppressWarnings(as.numeric(sub("^[%] for dose.*:\\s*", "", rx.isodose, ignore.case=TRUE)))
			if (is.na(dose.rx)) {
				warning("Prescription dose not specified")
			}

		    name <- sub("^.*:\\s*(.+$)", "\\1", data[grep("^Histogram.*:\\s*", data, ignore.case=TRUE, perl=TRUE)])
		    if (length(name) < 1) {
				warning("Invalid DVH file format, could not extract structure")
				return(new("DVH"))
		    }
		    volume <- suppressWarnings(as.numeric(sub("^.*:\\s*(.+$)", "\\1", data[grep("^Volume.*:\\s*", data, ignore.case=TRUE, perl=TRUE)], perl=TRUE)))
		    if (length(volume) < 1) {
				warning("Invalid DVH file format, could not extract structure volume information")
				return(new("DVH"))
		    }
        	header <- grep("Dose\\s*[(]\\s*(%|Gy|cGy)\\s*[)].*Volume", data, ignore.case=TRUE, perl=TRUE)
			if (grepl("^\\s*Dose\\s*[(]\\s*(Gy|cGy)\\s*[)].*Volume", data[header], ignore.case=TRUE, perl=TRUE)) {
				dose.type <- "absolute"
			}
			else {
				dose.type <- "relative"
			}
			if (grepl(".*Volume\\s*[(]\\s*[%]\\s*[)]", data[header], ignore.case=TRUE, perl=TRUE)) {
				volume.type <- "relative"
			}
			else {
				volume.type <- "absolute"
			}
			getDose <- function(dose) {
				if (grepl("[(]\\s*[%]\\s*[)]", dose)) {
					if (dose.type == "absolute") {
						dose <- suppressWarnings(as.numeric(sub(".*:\\s*", "", dose)))*dose.rx/rx.isodose	
					}
					else {
						dose <- suppressWarnings(as.numeric(sub(".*:\\s*", "", dose)))
					}
				}
				else {
					if (dose.type == "absolute") {
						dose <- suppressWarnings(as.numeric(sub(".*:\\s*", "", dose)))
					}
					else {
						dose <- suppressWarnings(as.numeric(sub(".*:\\s*", "", dose)))*rx.isodose/dose.rx
					}
				}
				return(dose)				
			}

			dose.min <- getDose(data[grep("^Dose minimum.*:\\s*", data, ignore.case=TRUE, perl=TRUE)])
			dose.max <- getDose(data[grep("^Dose maximum.*:\\s*", data, ignore.case=TRUE, perl=TRUE)])
			dose.mean <- max(0, getDose(data[grep("^Dose mean.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			if (verbose) {
				cat("  ..Importing structure: ", name, "  [volume: ", volume, "cc, dose: ", dose.min, " - ", dose.max, dose.units, "]\n", sep="")
			}

			dose.mode <- max(0, getDose(data[grep("^Dose modal.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.median <- max(0, getDose(data[grep("^Dose median.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)
			dose.STD <- max(0, getDose(data[grep("^Standard dev.*:\\s*", data, ignore.case=TRUE, perl=TRUE)]), na.rm=TRUE)

			con <- textConnection(data[(header+1):length(data)])
			dvh <- read.table(con, header=FALSE, stringsAsFactors=FALSE)
			close(con)
			data.dose <- dvh[, 1]
			data <- dvh[, 2]
			if (DVH.type == "differential") {
				data <- data * diff(c(-data.dose[1], data.dose))
				temp.doses <- data.dose - diff(c(-data.dose[1], data.dose))/2
				data.dose <- c(temp.doses, (2*data.dose - temp.doses)[length(temp.doses)])
				data <- diffinv(-data, xi=sum(data))
			}
			return(new("DVH", dose.min=dose.min, dose.max=dose.max, dose.mean=dose.mean, dose.mode=dose.mode, dose.median=dose.median, dose.STD=dose.STD, dose.rx=dose.rx, rx.isodose=rx.isodose, structure.name=name, structure.volume=volume, doses=data.dose, volumes=data, type="cumulative", dose.type=dose.type, dose.units=dose.units, volume.type=volume.type))	
		}
	)
	
	# RETURN DVH LIST
	names(DVH.list) <- unlist(lapply(DVH.list, names))
	return(new("DVH.list", DVH.list))
}


read.DVH.TomoTherapy <- function (file, verbose=TRUE) {
	if (!(fid <- file(file, open="r"))) {
		warning(paste("Could not open file '", file, "'", sep=""))		
		return()
	}
	if (verbose) {
		cat("Reading DVH file ('", file, "')... \n", sep="")
	}
	close(fid)
	data <- read.table(file, header=TRUE, sep=",", quote="\"")
	dvh <- list()
	for (i in 1:(dim(data)[2]/3)) {
		name <- sub("(.*)[.]STANDARD[.](.*)", "\\1\\2", colnames(data)[i*3-2], ignore.case=FALSE, perl=TRUE)
		dose.units <- toupper(sub("^Dose[.]*(c?Gy).*", "\\1", colnames(data)[i*3-1], ignore.case=TRUE, perl=TRUE))
		if (dose.units == "GY") {
			dose.units <- "Gy"
		}
		else if (dose.units == "CGY") {
			dose.units <- "cGy"
		}
		if (grepl("^Relative[.]", colnames(data)[i*3], ignore.case=TRUE, perl=TRUE)) {
			volume.type <- "relative"			
		}
		else {
			volume.type <- "absolute"
		}
		if (verbose) {
			cat("  ..Importing structure: ", name, "\n", sep="")
		}
		data.dose <- data[,i*3-1]
		data.volume <- data[,i*3]
		dvh <- c(dvh, new("DVH", structure.name=name, doses=data.dose, volumes=data.volume, dose.units=dose.units, volume.type=volume.type, type="cumulative"))
	}
	return(new("DVH.list", dvh))
}


read.DVH.Monaco <- function (file, verbose=TRUE) {
	if (!(fid <- file(file, open="r"))) {
		warning(paste("Could not open file '", file, "'", sep=""))		
		return()
	}
	if (verbose) {
		cat("Reading DVH file ('", file, "')... ", sep="")
	}
	header <- unlist(strsplit(readLines(fid, n=3),"[ ]*[|][ ]*", perl=TRUE))
	data <- readLines(fid)
	close(fid)
    patient <- sub("^.*: (.*)[~].*", "\\1", header[grep("^Patient ID.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
    ID <- sub("^.*: .*[~](.*)", "\\1", header[grep("^Patient ID.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
    date <- data[length(data)]
	plan <- sub("^Plan Name:\\s*(.+$)", "\\1", header[grep("^Plan Name.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
	dose.units <- sub("^Dose Units:\\s*(.+$)", "\\1", header[grep("^Dose Units.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
	if (dose.units == "%") {
		dose.type <- "relative"
		dose.units <- sub("^Bin Width:.*[(](.*)[)]", "\\1", header[grep("^Bin Width.*:", header, ignore.case=TRUE, perl=TRUE)])
	}
	else {
		dose.type <- "absolute"
	}
	volume.units <- sub("^Volume Units:\\s*(.+$)", "\\1", header[grep("^Volume Units.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
	if (volume.units == "%") {
		volume.type <- "relative"
	}
	else {
		volume.type <- "absolute"
	}

	if (verbose) {
		cat("[exported on ", date, "]\n", sep="")
		cat("  Patient: ", patient, " (", ID, ")\n", sep="")
		cat("  Plan: ", plan, "\n", sep="")		
	}
	data <- strsplit(data[1:(length(data)-3)],"[ ]+")
	data.structures <- unlist(lapply(data, function(x) {x[1]}))
	data.dose <- as.numeric(unlist(lapply(data, function(x) {x[2]})))
	data.volume <- as.numeric(unlist(lapply(data, function(x) {x[3]})))
	structures <- unique(data.structures)
	# EXTRACT DVH DATA FOR EACH STRUCTURE
	DVH.list <- lapply(structures,
		function (structure) {
			which.data <- which(data.structures == structure)
			which.dose <- data.dose[which.data]
			which.volume <- data.volume[which.data]	
			if (identical(which.volume,sort(which.volume,decreasing=TRUE))) {
				DVH.type <- "cumulative"
			}
			else {
				DVH.type <- "differential"
			}
			return(new("DVH", patient=patient, ID=ID, structure.name=structure, doses=which.dose, dose.units=dose.units, volumes=which.volume, type=DVH.type, dose.type=dose.type, volume.type=volume.type))	
		})

	return(DVH.list)
}


read.DVH.RayStation <- function (file, verbose=TRUE) {
	if (!(fid <- file(file, open="r"))) {
		warning(paste("Could not open file '", file, "'", sep=""))		
		return()
	}
	if (verbose) {
		cat("Reading DVH file ('", file, "')... ", sep="")
	}
	header <- readLines(fid, n=3)
	data <- readLines(fid)
	close(fid)

	# EXTRACT HEADER INFORMATION
    patient <- sub("^[#]PatientName.*:\\s*(.*)", "\\1", header[grep("^[#]PatientName.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
    ID <- sub("^[#]PatientId.*:\\s*(.*)", "\\1", header[grep("^[#]PatientId.*:\\s*", header, ignore.case=TRUE, perl=TRUE)])
	plan <- sub("^[#]Dosename:\\s*(.+$)", "\\1", header[grep("^[#]Dosename:\\s*", header, ignore.case=TRUE, perl=TRUE)])
	# EXTRACT PLAN DOSE (IN CGY)
	if (grepl("Plan dose:", plan, ignore.case=TRUE)) {
		dose.rx <- as.numeric(sub("^Plan dose:.+\\s([.0-9]+)c?Gy.*", "\\1", plan, perl=TRUE, ignore.case=TRUE))
		if (toupper(sub("^Plan dose:.+\\s[.0-9]+(c?Gy).*", "\\1", plan, perl=TRUE, ignore.case=TRUE)) == "GY") {
			dose.rx <- dose.rx * 100		
		}
	}
	else {
		dose.rx <- NA
	}
		
    # IDENTIFY STRUCTURES
    struct.start <- grep("^[#]RoiName:", data, perl=TRUE)
    struct.end <- struct.start + diff(c(struct.start, length(data)+1)) - 1
    if (length(struct.start) < 1) {
		warning(paste("File '", file, "' contained no recognizable DVH structure(s)", sep=""))
		if (verbose) {
			cat("ERROR\n")
		}
		return()
    }
    else if (length(struct.start) == 1) {
   	 	structures <- list(data[struct.start:struct.end])
    }
    else {
	    structures <- mapply(function(start, end) list(data[start:end]), struct.start, struct.end)
	}

	if (verbose) {
		cat("\n  Patient: ", patient, " (", ID, ")\n", sep="")
		cat("  Plan: ", plan, "\n", sep="")		
	}
	DVH.type <- "cumulative"
	# EXTRACT DVH DATA FOR EACH STRUCTURE
	DVH.list <- lapply(structures,
		function (data) {
			# EXTRACT STRUCTURE NAME
		    name <- sub("^[#]RoiName:\\s*(.+$)", "\\1", data[grep("^[#]RoiName:\\s*", data, ignore.case=TRUE, perl=TRUE)])
		    if (length(name) < 1) {
				warning("Invalid DVH file format, could not extract structure")
				return(new("DVH"))
		    }
			# EXTRACT DOSE UNITS
			dose.units <- toupper(sub("^[#]Dose unit:\\s*(.+$)", "\\1", data[grep("^[#]Dose unit:\\s*", data, ignore.case=TRUE, perl=TRUE)], perl=TRUE))
		    if (length(dose.units) < 1) {
				warning("Invalid DVH file format, could not extract dose units")
				return(new("DVH"))
		    }
			if (dose.units == "GY") {
				dose.units <- "Gy"
			}
			else if (dose.units == "CGY") {
				dose.units <- "cGy"
			}
			
		    volume <- suppressWarnings(as.numeric(sub("^.*:\\s*(.+)[%]$", "\\1", data[grep("^[#]Roi volume.*:\\s*", data, ignore.case=TRUE, perl=TRUE)], perl=TRUE)))
		    if ((length(volume) == 1) & (volume > 0))  {
				warning(paste(volume, "% volume of structure (", name, ") exists outside measureable dose grid, use DVH data with caution", sep=""))
		    }
			if (verbose) {
				cat("  ..Importing structure: ", name, "  [units: ", dose.units, "]\n", sep="")
			}
			con <- textConnection(data[4:length(data)])
			dvh <- read.table(con, header=FALSE, stringsAsFactors=FALSE)
			close(con)
			data.dose <- dvh[, 1]
			data <- dvh[, 2]
			return(new("DVH", structure.name=name, structure.volume=volume, doses=data.dose, volumes=data, type="cumulative", dose.type="absolute", dose.rx=if (dose.units == "Gy") {dose.rx*100} else {dose.rx}, dose.units=dose.units, volume.type="relative"))	
		}
	)

	# RETURN DVH LIST
	names(DVH.list) <- unlist(lapply(DVH.list, names))
	return(new("DVH.list", DVH.list))

}
