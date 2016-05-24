CQmodel <-
function(p.est = NULL, show = NULL, p.type = NULL, equation = NULL) {

	############Helper functions############
	

	breakup <- function(data, starts, titles) {
		ends = c(starts[-1], length(data))
		mapply(breakdown, titles, starts, ends, list(data), SIMPLIFY = FALSE)
	}


	breakdown <- function(title, m_sec1, m_sec2, shw) {
		shw[m_sec1:m_sec2]
	}


	RMP <- function(table, parts) {

		RMP.lengths = c(10, 9, 8, 6, 6, 6, 8, 6, 5, 6)
		RMP.titles = c("est", "error", "U.fit", "U.Low", "U.High", "U.T", "W.fit", "W.Low", "W.High", "W.T")

		titles <- as.list(as.data.frame(rbind(paste("n_", parts, sep = ""), parts), stringsAsFactors = FALSE))
		titles[parts == "step"] <- "step"
		titles <- unlist(titles)


		left.side.titles <- strsplit(table[5], "ESTIMATE")[[1]][1]

		parts.search <- parts

		parts.search[parts == 'step'] <- '(step|category)'

		line.seps <- numeric()
		
		for (i in 1:length(parts)) {

			temp.col.seps <- numeric()

			if (gregexpr(parts.search[i],left.side.titles)[[1]][1] - 2 > 0){

				temp.col.seps[1] <- gregexpr(parts.search[i],left.side.titles)[[1]][1] - 2
				#print(temp.col.seps)
				
				

				if (i > 1){

					# temp.col.seps[1] <- temp.col.seps[1] 
					# print(temp.col.seps)
	
					if(parts[i] != 'step' & parts[i] != 'category'){
						
						
	
						first.line <- substr(table[7], 1, temp.col.seps[1]+2)
						#print(first.line)
	
						number.col <- (nchar(first.line) - gregexpr("[0-9]+\\s+?$",first.line)[[1]][1] ) 
						#print(number.col)
	
						temp.col.seps[2] <- number.col
						temp.col.seps[1] <- temp.col.seps[1] - temp.col.seps[2] - sum(line.seps)
						
						#print(temp.col.seps)
						
						
	
					}else{

						temp.col.seps[1] <- temp.col.seps[1] - sum(line.seps)

					}



				}



			}



			if (i == 1){
	
				line.seps <- temp.col.seps
	
			}else{
	
				line.seps <- c(line.seps, temp.col.seps)
	
			}
	

		}

		if (length(line.seps)==0){

			line.seps <- nchar(left.side.titles)

		} else{

			line.seps <- c(line.seps,nchar(left.side.titles) - max(line.seps))

		}

		if (all(parts == 'step') & length(parts) == 1) {line.seps <- nchar(left.side.titles)}

		out = grep("^ *(Parameter )?[0-9]", table, value = TRUE)		
		out = gsub("[\\(\\)\\*,]", " ", out)
		out = trim(out)
		
		if(nchar(out[1]) < sum(RMP.lengths)) {
			out1 <- split.right(out, 10)
			multcols <- grep("[^ ] +[^ ]",out1[[2]])
			if(length(multcols) > 0) {
				RMP.lengths <- c(10,7)
				RMP.titles <- c("est","error")
			} else {
				RMP.lengths <- 10
				RMP.titles <- "est"
			}
		}
		out <- split.right(out, sum(RMP.lengths))
		
		if ( (all(parts == 'step') & length(parts) == 1) | imported) {

			left.table <- matrix(apply(out[1],1,function (x) gsub("^\\s+|\\s+$", "", x)), ncol=1)			
			colnames(left.table) <- parts

		}else{
			left.table <- read.fwf(tempify(out[1]), line.seps, col.names = titles, stringsAsFactors = FALSE)
			left.table[sapply(left.table, is.character)] <- sapply(left.table[sapply(left.table, is.character)],function (x) gsub("^\\s+|\\s+$", "", x))

		}


		right.table <- read.fwf(tempify(out[2]), RMP.lengths, col.names = RMP.titles, stringsAsFactors = FALSE)

		cbind(left.table, right.table)

	}


	split.right <- function(table, right) {
		maxchar <- max(sapply(table,nchar))
		left = maxchar - right
		tf <- tempfile()
		write(table, tf)
		out <- read.fwf(tf, c(left, right), stringsAsFactors=FALSE)
		out
	}

	tempify <- function(list) {
		tf <- tempfile()
		write(as.character(as.matrix(list)), tf)
		tf
	}
	
		make.GIN <- function(table) {
		if (length(table) == 5) {
			items <- as.vector(unlist(unique(table[5])))
			if(class(items) != "character")
				items <- paste("Item",items,sep="_")
			#print(items)
			#print(length(items))
			out <- mapply(by.item,items,c(1:length(items)), list(table[4]), list(table[2]), list(max(table[1])),SIMPLIFY=FALSE)
			#print(as.data.frame(out))
			return(t(as.data.frame(out)))
		} else {
			pieces <- as.character(unlist(unique(table[5])))
			return(mapply(break.GIN, pieces, list(table),SIMPLIFY=FALSE))
		}
	}

	break.GIN <- function(piece, table) {
		table <- table[table[5] == piece,]
		#print(make.GIN(table[-c(3, 4, 5)]))
		return(make.GIN(table[-c(3, 4, 5)]))
	}

	by.item <- function(name,item, ids, values, how.long) {

		vals <- values[ids == item]
		if (length(vals) < how.long) {
			extras <- c((length(vals) + 1):how.long)
			vals[extras] <- NA
		}
		vals

	}

	get.names <- function(type, table) {
		return(unique(table[type]))
	}

	split.by <- function(table, split) {
		out <- read.table(tempify(table), sep = split, stringsAsFactors = FALSE, strip.white = TRUE)
		row.names(out) <- t(rename(t(out[1])))
		out <- out[out$V2 != "", ]
		out <- as.list(as.data.frame(t(out[2]), stringsAsFactors = FALSE))
		out
	}

	trim <- function(x) gsub("^\\s+|\\s+$", "", x)

	safe.remove <- function(table, line.at) {
		if (length(line.at) > 0) {
			table <- table[-line.at]
		}
		table
	}

	make.label <- function(lab, dim) {
		return(paste(lab, " (", dim, ")", sep = ""))
	}

	make.labels <- function(labs, dims) {
		return(unlist(as.list(t(sapply(labs, make.label, dims)))))
	}

	reliabilities <- function(table) {
		table <- table[-grep("^Dimension:", table)]
		colons.at <- grep(":", table)
		colons <- table[colons.at]
		out <- split.by(colons, ":")
		out
	}

	numify <- function(table) {
		out <- table
		nums <- suppressWarnings(!is.na(as.numeric(table)))
		out[nums] <- as.numeric(table[nums])
		return(out)

	}

	rename <- function(titles) {
		titles[titles == "SUMMARY OF THE ESTIMATION"] <- "SOE"
		titles[titles == "TABLES OF RESPONSE MODEL PARAMETER ESTIMATES"] <- "RMP"
		titles[grepl("^IMPORTED MODEL", titles)] <- "RMP"
		titles[titles == "TABLES OF POPULATION MODEL PARAMETER ESTIMATES"] <- "PMP"
		titles[grepl("MAP OF .+ AND RESPONSE MODEL PARAMETER ESTIMATES", titles)] <- "MRM"
		titles[grepl("MAP OF .+ AND THRESHOLDS", titles)] <- "MTH"
		titles[titles == "TABLES OF GIN Thresholds"] <- "GIN"
		titles[titles == "TABLES OF GIN Item Parameters"] <- "GIN.deltas"
		titles[titles == "Estimation method was"] <- "method"
		titles[titles == "Assumed population distribution was"] <- "distribution"
		titles[titles == "Constraint was"] <- "constraint"
		titles[titles == "The format"] <- "format"
		titles[titles == "The item model"] <- "equation"
		titles[titles == "Sample size"] <- "participants"
		titles[titles == "Final Deviance"] <- "deviance"
		titles[titles == "Total number of estimated parameters"] <- "parameters"
		titles[titles == "The number of iterations"] <- "iterations"
		titles[titles == "Termination criteria"] <- "criteria"
		titles[titles == "Random number generation seed"] <- "seed"
		titles[titles == "Number of nodes used when drawing PVs"] <- "PV.nodes"
		titles[titles == "Number of nodes used when computing fit"] <- "fit.nodes"
		titles[titles == "Number of plausible values to draw"] <- "n.plausible.values"
		titles[titles == "Maximum number of iterations without a deviance improvement"] <- "max.iterations.no.improvement"
		titles[titles == "Maximum number of Newton steps in M-step"] <- "max.steps"
		titles[titles == "Value for obtaining finite MLEs for zero/perfects"] <- "zero.perfect.value"
		titles[titles == "The regression model"] <- "regression"
		titles[titles == "key 1 scored as 1"] <- "key"
		titles[titles == "The Data File"] <- "data.file"

		titles[titles == "Max iterations"] <- "max.iterations"
		titles[titles == "Parameter Change"] <- "parameter.change"
		titles[titles == "Deviance Change"] <- "deviance.change"

		titles[titles == "REGRESSION COEFFICIENTS"] <- "reg.coef"
		titles[grepl("COVARIANCE/CORRELATION MATRIX", titles)] <- "cov.cor"
		titles[titles == "RELIABILITY COEFFICIENTS"] <- "rel.coef"

		titles
	}


	model <- list()
	if (!(is.null(show))) {
		#ptm <- proc.time()
		shw <- readLines(show)
		shw.starts = grep("^\f==", shw)
		if (length(shw.starts) == 0) {
			shw.starts = grep("^==", shw)
			shw.starts = shw.starts[(shw.starts + 3) %in% shw.starts]
			CQV = 3
		} else {
			CQV = 2
		}
		titles = shw[shw.starts + 2]
		if (length(grep("^IMPORTED MODEL", titles)) == 0) {
			imported = FALSE
		} else {
			imported = TRUE
		}

		shw.titles = rename(titles)


		model <- breakup(shw, shw.starts, shw.titles)
		model$imported <- imported
		
		#print(model)

		
		#######SOE####################
		
		#ptm <- proc.time()
		date.pattern <- "(?:(Sun|Mon|Tue|Wed|Thu|Fri|Sat)\\s+)?(Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\\s+(0[1-9]|[1-2]?[0-9]|3[01])\\s+(2[0-3]|[0-1][0-9]):([0-5][0-9])(?::(60|[0-5][0-9]))?\\s+(19[0-9]{2}|[2-9][0-9]{3})+$"
		date.at <- grep(date.pattern, model[[1]])
		title.date <- model[[1]][date.at]
		model$SOE <- safe.remove(model$SOE, date.at)
		m <- regexpr(date.pattern, title.date)
		model$run.details <- list()
		#class(model$run.details) <- "details"
		date.string <- regmatches(title.date, m)
		model$run.details$date <- strptime(date.string, format = "%a %b %d %H:%M %Y")
		model$title <- trim(paste(unlist(regmatches(title.date, m, invert = TRUE)), collapse = ""))


		file.at <- grep("The Data File: ", model$SOE)
		file.line <- model$SOE[file.at]
		model$SOE <- safe.remove(model$SOE, file.at)
		m <- regexpr("The Data File: ", file.line)
		model$run.details$data.file <- trim(paste(unlist(regmatches(file.line, m, invert = TRUE)), collapse = ""))

		colons.at <- grep(":", model$SOE)
		other.lines <- model$SOE[-colons.at]
		colons <- model$SOE[colons.at]
		two.colons.at <- grep(":.+:", colons)
		other.lines <- append(other.lines, colons[two.colons.at])
		colons <- safe.remove(colons, two.colons.at)
		
		if(!is.null(colons)) {

			SOE <- split.by(colons, ":")
	
			other.lines <- other.lines[-grep("===", other.lines)]
			other.lines <- other.lines[other.lines != ""]
			other.lines <- other.lines[other.lines != "SUMMARY OF THE ESTIMATION"]
			other.lines <- trim(other.lines)
	
			reason.at <- grep("Iterations terminated because", other.lines)
			SOE$termination.reason <- other.lines[reason.at]
			other.lines <- safe.remove(other.lines, reason.at)
	
			deviance.line <- grep("Deviance Change=", other.lines, value = TRUE)
			other.lines <- other.lines[-grep("Deviance Change=", other.lines)]
			termination.criteria <- c(unlist(strsplit(SOE$criteria, ",")), deviance.line)
			SOE <- c(SOE, split.by(termination.criteria, "="))
			SOE$criteria <- NULL
	
			model$run.details$format <- SOE$format
			model$run.details$key <- SOE$key
			model$run.details <- c(model$run.details, other.lines)
			SOE <- numify(SOE)
	
			model <- c(SOE["equation"], SOE["participants"], SOE["deviance"], SOE["parameters"], model)
			class(SOE) <- "SOE"
	
	
			model$SOE <- SOE
		
		}
		
		if(is.null(model$equation)) {
			if(is.null(equation)) {
				stop("Please specify the model equation.")
			} else {
				model$equation <- equation
			}
		}
		
		
		
		##############RMP######################
		


		additive.parts = unlist(strsplit(model$equation, "[+|-]"))
		parts = strsplit(additive.parts, "\\*")
		model$additive.parts <- additive.parts
		model$parts <- parts
		
		if (imported) {
			params <- RMP(model$RMP, "Parameters")
			model$RMP <- list()
			model$RMP$item <- params
		} else {
			RMP.tables <- breakup(model$RMP, grep("  VARIABLES", model$RMP) - 2, additive.parts)
			model$RMP = mapply(RMP, RMP.tables, parts, SIMPLIFY = FALSE)
			model$run.details$names <- mapply(get.names, parts[parts == additive.parts], model$RMP[parts == additive.parts])
		}
		
		
		
		#print(model$RMP)
		

		##########PMP###########
		if(!is.null(model$PMP)) {
		PMP.starts <- grep("^====+", model$PMP)
		PMP.heads <- grep(date.string, model$PMP)
		PMP.starts <- PMP.starts[!(PMP.starts + 1) %in% PMP.heads]
		PMP.titles <- rename(trim(paste(model$PMP[PMP.starts + 1], model$PMP[PMP.starts + 2], sep = "")))
		PMP <- breakup(model$PMP, PMP.starts, PMP.titles)
		variance.line = grep("Variance", PMP$cov.cor, value = TRUE)
		m <- gregexpr("[0-9]+\\.[0-9]+(?![0-9])(?![\\s]*\\))", variance.line, perl = TRUE)
		PMP$variances <- as.numeric(unlist(regmatches(variance.line, m)))
		m <- gregexpr("[0-9\\.\\-#IO]+(?=[\\s]*\\))", variance.line, perl = TRUE)
		errors <- as.numeric(unlist(regmatches(variance.line, m)))
		PMP$nDim = length(PMP$variances)
		length(errors) <- length(PMP$variances)
		PMP$dimensions = "Main dimension"
		if (PMP$nDim > 1) {

			if (CQV == 2){
				
				PMP$cov.cor <- read.fwf(tempify(PMP$cov.cor[8:(8 + PMP$nDim - 1)]), widths = c(25, rep(9, PMP$nDim)), row.names = 1, 
				strip.white = TRUE, stringsAsFactors = FALSE)
				
			}
			
			if (CQV == 3){
				
				PMP$cov.cor <- read.fwf(tempify(PMP$cov.cor[8:(8 + PMP$nDim - 1)]), widths = c(25, rep(18, PMP$nDim)), row.names = 1, 
				strip.white = TRUE, stringsAsFactors = FALSE)
				
			}

			PMP$dimensions <- row.names(PMP$cov.cor)
			names(PMP$cov.cor) <- PMP$dimensions
			PMP$cor.matrix <- PMP$cov.cor
			PMP$cor.matrix[upper.tri(PMP$cor.matrix)] = t(PMP$cor.matrix)[upper.tri(t(PMP$cor.matrix))]
			diag(PMP$cor.matrix) = 1
			PMP$cov.matrix <- PMP$cov.cor
			PMP$cov.matrix[lower.tri(PMP$cov.matrix)] = t(PMP$cov.matrix)[lower.tri(t(PMP$cov.matrix))]
			diag(PMP$cov.matrix) <- PMP$variances
		}
		if (length(errors) > 0) {
			PMP$variances <- cbind(PMP$variances, errors)
		}
		PMP$cov.cor <- NULL

		PMP$reg.coef = gsub("[\\(\\),]", " ", PMP$reg.coef)
		start <- grep("CONSTANT", PMP$reg.coef)
		end <- grep("^--+$", PMP$reg.coef) - 1
		PMP$reg.coef <- read.table(tempify(PMP$reg.coef[start:end]), stringsAsFactors = FALSE, strip.white = TRUE, row.names = 1)
		if (ncol(PMP$reg.coef) == PMP$nDim) {
			names(PMP$reg.coef) <- PMP$dimensions
		} else if (ncol(PMP$reg.coef) == 2 * PMP$nDim) {
			names(PMP$reg.coef) <- unlist(as.list(as.data.frame(rbind(PMP$dimensions, "S. errors"), stringsAsFactors = FALSE)))
		}
		PMP$reg.coef <- t(PMP$reg.coef)

		rel.starts <- grep("Dimension:", PMP$rel.coef)
		PMP$rel.coef <- breakup(PMP$rel.coef, rel.starts, PMP$dimensions)
		PMP$rel.coef <- sapply(PMP$rel.coef, reliabilities)
		PMP$rel.coef <- numify(PMP$rel.coef)
		PMP$rel.coef[PMP$rel.coef == "Unavailable"] <- NA
		PMP$rel.coef <- t(PMP$rel.coef)

		model$PMP <- NULL
		model <- c(model, PMP)
		
		}
		


		########GIN#########
		
		do.GIN <- function(GIN) {
			
			GIN <- GIN[4:(length(GIN))]
			GIN <- GIN[-grep("===",GIN)]
			GIN <- GIN[-grep("GIN Number", GIN)]
			GIN <- GIN[-grep("---",GIN)]

			GIN <- gsub("\\s*\t\\s*", "\t", GIN)
			GIN <- gsub(" +","_",GIN)
			GIN <- gsub("^[0-9\\.]+\\.", "", GIN)
			GIN <- read.delim(tempify(GIN),header=FALSE)
			
			GIN <- GIN[colSums(!is.na(GIN))!=0]
			return(make.GIN(GIN))
		}
		
		
		if (!is.null(model$GIN)) {
			model$GIN <- do.GIN(model$GIN)
		}
		
		
		if(!is.null(model$GIN.deltas)) {
			model$GIN.deltas <- do.GIN(model$GIN.deltas)
		}
		
		#return(proc.time()-ptm)
		class(model) <- "CQmodel"

	}



	#######Person Parameters#############
	
	if (!is.null(p.est)) {

		if (is.null(p.type)) {
			eaps <- unlist(strsplit(p.est, "[.]"))
			p.type <- toupper(eaps[length(eaps)])
		}

		p.est <- na.omit(read.table(p.est, stringsAsFactors = FALSE, fill = TRUE))

		if (p.type == "EAP") 
			colperdim = 3
		else colperdim = 4

		model$nDim = floor((length(p.est) - 1)/colperdim)

		if (is.null(model$dimensions) || length(model$dimensions != model$nDim)) 
			model$dimensions <- paste("d", c(1:model$nDim), sep = "")

		if (p.type == "EAP") {
			pp_lab_t <- make.labels(c("est", "error", "pop"), model$dimensions)
		} else {
			pp_lab_t <- c(make.labels(c("sscore", "max"), model$dimensions), make.labels(c("est", "error"), model$dimensions))
		}


		if (length(p.est)%%colperdim == 1) {
			names(p.est) <- c("casenum", pp_lab_t)
		} else if (p.type != "WLE") {
			names(p.est) <- c("casenum", "pid", pp_lab_t)
		} else if (length(p.est)%%colperdim == 3) {
			names(p.est) <- c("casenum","pid", pp_lab_t, "case.fit")
		} else if (any(p.est[4] < p.est[3])) {
			names(p.est) <- c("casenum", pp_lab_t, "case.fit")
		} else {
			names(p.est) <- c("casenum", "pid", pp_lab_t)
		}
		
		

		model$p.est <- p.est
		model$p.est.type <- p.type

	}

	class(model) <- "CQmodel"
	return(model)

}
