setAs("structure3D", "DVH", 
	function(from) {
		return(from$DVH)
	}
)

setAs("zDVH", "DVH", 
	function(from) {
		from$volumes <- apply(from$volumes, 1, sum)
		class(from$volumes) <- "numeric"
		class(from) <- "DVH"
		return(from)
	}
)


setMethod("$", "DVH",
	function (x, name) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			return(NULL)	
		}
		else {
			return(slot(x, name))	
		}
	}
)

setMethod("names", "DVH",
	function (x) {
		return(x$structure.name)
	}
)

setMethod("names<-", "DVH",
 	function (x, value) {
 		x$structure.name <- value
 		return(x)
 	}
)

setMethod("$<-", "DVH",
	function (x, name, value) {
		if (inherits(try(slot(x, name), silent=TRUE), "try-error")) {
			warning("'", name, "' is not a parameter in class 'DVH'")
		}
		else {
			slot(x, name) <- value
		}
		return(x)
	}
)

setMethod("[", "DVH",
	function (x, i, ...) {
		if (!validObject(x)) {
			stop("not a valid object of 'DVH'")
		}
		if (!missing("i")) {
			if (is.null(i)) return(NULL)
			if (length(i) < 1) return(numeric())
			if (is.logical(i)) {
				if (length(i) != length(x@doses)) {
					stop("(subscript 'i') logical subscript does not match length of 'DVH' doses")
				}
				i <- which(i)
			}
			result <- c()
			result.units <- c()
			i <- input <- toupper(as.character(i))
			type <- sub("(V|D).*", "\\1", i)
			volume <- grepl("VOL", i)
			if (any(volume)) {
				type[volume] <- "VOLUME"
			}
			patients <- grepl("PATIENT", i)
			if (any(patients)) {
				type[patients] <- "PATIENT"
			}
			IDs <- grepl("ID", i)
			if (any(IDs)) {
				type[IDs] <- "ID"
			}
			value <- sub("[VD]([-<.0-9]+|MAX|MIN|MEAN|MEDIAN|RX|INTEGRAL).*", "\\1", i)
			type2 <- sub("[VD]([-<.0-9]+|MAX|MIN|MEAN|MEDIAN|RX|INTEGRAL)([%]|GY|CGY|CC)*.*$", "\\2", i)
			type3 <- grepl(".*[(](.*)[)]$", i)
			type4 <- sub(".*[(](.*)[)]$", "\\1", i)
			type5 <- grepl("[-<]", value)
			type6 <- grepl(" (=|==|<=|>=|<|>|!=) [.0-9]+", i)
			type7 <- sub(".* (=|==|<=|>=|<|>|!=) [.0-9]+.*", "\\1", i)
			type8 <- suppressWarnings(as.numeric(sub(".* (=|==|<=|>=|<|>|!=) ([.0-9]+).*([(].*[)])*$", "\\2", i)))
			for (count in 1:length(i)) {
				switch(type[count],
					PATIENT = {
						result <- c(result, x$patient)
						result.units <- c(result.units, "")	
					},
					ID = {
						result <- c(result, x$ID)
						result.units <- c(result.units, "")	
					},
					VOLUME = {
						if (type4[count] == "%") {
							result <- c(result, 100)	
							result.units <- c(result.units, "%")		
						}
						else {
							result <- c(result, x@structure.volume)
							result.units <- c(result.units, "cc")
						}
					},
					V = {
						if (type5[count]) {
							values <- suppressWarnings(as.numeric(unlist(strsplit(value[count], "[-<]"))))
							if (length(values) != 2) {
								warning("Improper format '", input[count], "' (dose must be specified as numeric range, e.g. 'V10-20Gy' or 'V<500cGy')")
								result <- c(result, NA)
								result.units <- c(result.units, NA)
								next
							}
							if (is.na(values[1])) {
								values[1] <- 0
							}
							if (x@type == "differential") {
								values <- convert.DVH(x, type="cumulative")[paste("V", values, type2[count], if (type3[count]) { paste("(", type4[count], ")", sep="")}, sep="")]
							}
							else {
								values <- x[paste("V", values, type2[count], if (type3[count]) { paste("(", type4[count], ")", sep="")}, sep="")]
							}
							result <- c(result, as.numeric(values[1])-as.numeric(values[2]))
							result.units <- c(result.units, names(values)[1])
							next
						}
						else {
							value[count] <- suppressWarnings(as.numeric(value[count]))
						}
						if (is.na(value[count])) {
							warning("Improper format '", input[count], "' (dose must be numeric, e.g. 'V20Gy')")
							result <- c(result, NA)
							result.units <- c(result.units, NA)
							next
						}
						if (length(x@doses) < 1) {
							warning("Empty DVH data")
							result <- c(result, NA)
							result.units <- c(result.units, NA)
							next
						}
						switch(x@type,
							cumulative = {
								switch(x@dose.type,
									absolute = {
										switch(type2[count],
											CGY = if (x@dose.units == "cGy") { TRUE } else { value[count] <- as.numeric(value[count]) / 100 },
											"%" = value[count] <- as.numeric(value[count]) * x@dose.rx / x@rx.isodose,
											GY = if (x@dose.units == "Gy") { TRUE } else { value[count] <- as.numeric(value[count]) * 100 },
											CC = value[count] <- NA,
											value[count] <- NA
										)
									},
									relative = {
										switch(type2[count],
											"%" = TRUE,
											GY = if (x@dose.units == "Gy") { value[count] <- as.numeric(value[count]) * x@rx.isodose / x@dose.rx } else { value[count] <- (as.numeric(value[count]) * x@rx.isodose / x@dose.rx) * 100 },
											CGY = if (x@dose.units == "cGy") { value[count] <- as.numeric(value[count]) * x@rx.isodose / x@dose.rx } else { value[count] <- as.numeric(value[count]) * x@rx.isodose * 0.01 / x@dose.rx },
											CC = value[count] <- NA,
											value[count] <- NA
										)
									},
									value[count] <- NA
								)
								if (is.na(value[count])) {
									warning("Improper format '", input[count], "' (should specify dose as % or cGy or Gy, e.g. 'V20Gy')")
									result <- c(result, NA)
									result.units <- c(result.units, NA)
									next
								}
								switch(type4[count],
									"%" = {
										if (x@volume.type == "absolute") {
											result <- c(result, 100 * approx(x@doses, x@volumes, value[count], yright=0, ties=max)$y / x@structure.volume)
										}
										else {
											result <- c(result, approx(x@doses, x@volumes, value[count], yright=0, ties=max)$y)
										}
										result.units <- c(result.units, "%")										
									},
									CC = {
										if (x@volume.type == "relative") {
											result <- c(result, approx(x@doses, x@volumes, value[count], yright=0, ties=max)$y * x@structure.volume / 100)
										}
										else {
											result <- c(result, approx(x@doses, x@volumes, value[count], yright=0, ties=max)$y)
										}
										result.units <- c(result.units, "cc")
									},
									{
										if (type3[count]) {
											warning("Improper format '", input[count], "' (should specify output volume as % or cc, e.g. 'V__(cc)')")
										}
										if (x@volume.type == "absolute") {
											result <- c(result, approx(x@doses, x@volumes, value[count], yright=0, ties=max)$y)
											result.units <- c(result.units, "cc")
										}
										else {
											result <- c(result, approx(x@doses, x@volumes, value[count], yright=0, ties=max)$y)
											result.units <- c(result.units, "%")										
										}
									}
								)
							},
							differential = {
								warning("No method available to extract volume given differential doses")
								result <- c(result, NA)
								result.units <- c(result.units, NA)
							}
						)
					},
					D = {
						switch(x@type,
							cumulative = {
								if (value[count] %in% c("MAX", "MIN", "MEAN", "MEDIAN", "RX")) {
									switch(value[count],
										MAX = value[count] <- x@dose.max,
										MIN = value[count] <- x@dose.min,
										MEAN = value[count] <- x@dose.mean,
										MEDIAN = value[count] <- x@dose.median,
										RX = if (x@dose.type == "absolute") { value[count] <- x@dose.rx } else { value[count] <- 100 * 100 / x@rx.isodose }
									)
									if ((type2[count] == "%") | (type4[count] == "%")) {
										if (x@dose.type == "relative") {
											result <- c(result, as.numeric(value[count]))
										}
										else {
											result <- c(result, as.numeric(value[count]) * x@rx.isodose / x@dose.rx)
										}
										result.units <- c(result.units, "%")
									}
									else if (type4[count] == "CGY") {
										if (x@dose.type == "relative") {
											value[count] <- as.numeric(value[count]) * x@dose.rx / x@rx.isodose
										}
										if (x@dose.units == "Gy") {
											result <- c(result, as.numeric(value[count]) * 100)
										}
										else {
											result <- c(result, as.numeric(value[count]))
										}
										result.units <- c(result.units, "cGy")
									}
									else if (type4[count] == "GY") {
										if (x@dose.type == "relative") {
											value[count] <- as.numeric(value[count]) * x@dose.rx / x@rx.isodose
										}
										if (x@dose.units == "cGy") {
											result <- c(result, as.numeric(value[count]) / 100)
										}
										else {
											result <- c(result, as.numeric(value[count]))
										}
										result.units <- c(result.units, "Gy")
									}
									else {
										if (type3[count]) {
											warning("Improper format '", input[count], "' (should specify output dose as %, cGy or Gy, e.g. 'Dmax(cGy)')")
										}
										result <- c(result, as.numeric(value[count]))
										result.units <- c(result.units, x@dose.units)
									}
									next
								}
								if (length(x@doses) < 1) {
									warning("Empty DVH data")
									result <- c(result, NA)
									result.units <- c(result.units, NA)
									next
								}																
								if (value[count] == "INTEGRAL") {
									if (type3[count]) {
										if (grepl("(>|<)[.0-9]+([%]|GY|CGY)*$", type4[count])) {
											if (grepl(">", type4[count])) {
												start.i <- as.numeric(sub("(>|<)([.0-9]+)([%]|GY|CGY)*$", "\\2", type4[count]))
												end.i <- Inf
											}
											else {
												start.i <- 0
												end.i <- as.numeric(sub("(>|<)([.0-9]+)([%]|GY|CGY)*$", "\\2", type4[count]))
											}

											units.i <- sub("(>|<)[-.0-9]+([%]|GY|CGY)$", "\\2", type4[count])
										}
										else if (grepl("[.0-9]+-[.0-9]+([%]|GY|CGY)*$", type4[count])) {
											start.i <- as.numeric(sub("([.0-9]+)[-].*", "\\1", type4[count]))
											end.i <- as.numeric(sub(".*[-]([.0-9]+)[^.0-9]*", "\\1", type4[count]))
											if (end.i < start.i) {
												units.i <- end.i
												end.i <- start.i
												start.i <- units.i
											}
											units.i <- sub("[-.0-9]+([%]|GY|CGY)$", "\\1", type4[count])
										}
										else {
											start.i <- 0
											end.i <- Inf
											units.i <- ""
										}
									}
									else {
										start.i <- 0
										end.i <- Inf
										units.i <- ""
									}
									switch(units.i,
										CGY = y <- convert.DVH(x, volume="absolute", dose="absolute", dose.units="cGy"),
										GY = y <- convert.DVH(x, volume="absolute", dose="absolute", dose.units="Gy"),
										"%" = y <- convert.DVH(x, volume="absolute", dose="relative"),
										y <- convert.DVH(x, volume="absolute", dose="absolute")
									)
									dose.bins <- diff(y@doses)
									y <- convert.DVH(y, type="differential")
									start.i <- max(start.i, min(y))				
									end.i <- min(end.i, max(y))				
									if (units.i == "%") {
										result <- c(result,
											y@dose.rx * max(0, integrate(function(dose) {
												return(approx(y@doses, y@volumes*y@doses/dose.bins, dose, yleft=0, yright=0, ties=max)$y)
												}, start.i, end.i, stop.on.error=FALSE, abs.tol=0, rel.tol=100*.Machine$double.eps
											)$value) / y@rx.isodose
										)
									}
									else {
										result <- c(result,
											max(0, integrate(function(dose) {
												return(approx(y@doses, y@volumes*y@doses/dose.bins, dose, yleft=0, yright=0, ties=max)$y)
												}, start.i, end.i, stop.on.error=FALSE, abs.tol=0, rel.tol=100*.Machine$double.eps
											)$value)
										)
									}
									result.units <- c(result.units, paste(y@dose.units, "*cc", sep=""))
									next
								}
								value[count] <- suppressWarnings(as.numeric(value[count]))
								if (is.na(value[count])) {
									warning("Improper format '", input[count], "' (volume must be numeric, e.g. 'D20%')")
									result <- c(result, NA)
									result.units <- c(result.units, NA)
									next
								}
								switch(x@volume.type,
									absolute = {
										switch(type2[count],
											CC = TRUE,
											"%" = value[count] <- as.numeric(value[count]) * x@structure.volume / 100,
											GY = value[count] <- NA,
											CGY = value[count] <- NA,
											value[count] <- NA
										)									
									},
									relative = {
										switch(type2[count],
											"%" = TRUE,
											CC = value[count] <- as.numeric(value[count]) * 100 / x@structure.volume,
											GY = value[count] <- NA,
											CGY = value[count] <- NA,
											value[count] <- NA
										)
									}
								)
								value.count <- as.numeric(value[count])
								if (is.na(value.count)) {
									warning("Improper format '", input[count], "' (should specify volume as % or cc, e.g. 'D__%')")
									result <- c(result, NA)
									result.units <- c(result.units, NA)
									next
								}
								else {
									switch(type4[count],
										"%" = {
											if (value.count > 100) {
												warning("Requested value for volume '", input[count], "' exceeds 100%")
											}
											if (value.count < 0) {
												warning("Requested value for volume '", input[count], "' is less than 0%")
											}
											result <- c(result, x@rx.isodose * min(x@dose.max, max(x@dose.min, approx(x@volumes, x@doses, value.count, ties=max)$y, na.rm=TRUE), na.rm=TRUE) / x@dose.rx)
											result.units <- c(result.units, "%")										
										},
										CGY = {
											if (x@dose.units == "Gy") {
												result <- c(result, min(x@dose.max, max(x@dose.min, approx(x@volumes, x@doses, value.count, ties=max)$y, na.rm=TRUE), na.rm=TRUE) * 100)
											}
											else {
												result <- c(result, min(x@dose.max, max(x@dose.min, approx(x@volumes, x@doses, value.count, ties=max)$y, na.rm=TRUE), na.rm=TRUE))
											}
											result.units <- c(result.units, "cGy")
										},
										GY = {
											if (x@dose.units == "cGy") {
												result <- c(result, min(x@dose.max, max(x@dose.min, approx(x@volumes, x@doses, value.count, ties=max)$y, na.rm=TRUE), na.rm=TRUE) / 100) 
											}
											else {
												result <- c(result, min(x@dose.max, max(x@dose.min, approx(x@volumes, x@doses, value.count, ties=max)$y, na.rm=TRUE), na.rm=TRUE))
											}
											result.units <- c(result.units, "Gy")
										},
										{
											if (type3[count]) {
												warning("Improper format '", input[count], "' (should specify output dose as %, cGy or Gy, e.g. 'D__(cGy)')")
											}
											result <- c(result, min(x@dose.max, max(x@dose.min, approx(x@volumes, x@doses, value.count, ties=max)$y, na.rm=TRUE), na.rm=TRUE))	
											if (x@dose.type == "absolute") {
												result.units <- c(result.units, x@dose.units)
											}
											else {
												result.units <- c(result.units, "%")										
											}
										}		
									)
								}
							},
							differential = {
								if (length(x@doses) < 1) {
									warning("Empty DVH data")
									result <- c(result, NA)
									result.units <- c(result.units, NA)
									next
								}								
								if (value[count] == "INTEGRAL") {
									if (type3[count]) {
										if (grepl("(>|<)[.0-9]+([%]|GY|CGY)*$", type4[count])) {
											if (grepl(">", type4[count])) {
												start.i <- as.numeric(sub("(>|<)([.0-9]+)([%]|GY|CGY)*$", "\\2", type4[count]))
												end.i <- Inf
											}
											else {
												start.i <- 0
												end.i <- as.numeric(sub("(>|<)([.0-9]+)([%]|GY|CGY)*$", "\\2", type4[count]))
											}

											units.i <- sub("(>|<)[-.0-9]+([%]|GY|CGY)$", "\\2", type4[count])
										}
										else if (grepl("[.0-9]+-[.0-9]+([%]|GY|CGY)*$", type4[count])) {
											start.i <- as.numeric(sub("([.0-9]+)[-].*", "\\1", type4[count]))
											end.i <- as.numeric(sub(".*[-]([.0-9]+)[^.0-9]*", "\\1", type4[count]))
											if (end.i < start.i) {
												units.i <- end.i
												end.i <- start.i
												start.i <- units.i
											}
											units.i <- sub("[-.0-9]+([%]|GY|CGY)$", "\\1", type4[count])
										}
										else {
											start.i <- 0
											end.i <- Inf
											units.i <- ""
										}
									}
									else {
										start.i <- 0
										end.i <- Inf
										units.i <- ""
									}
									switch(units.i,
										CGY = y <- convert.DVH(x, volume="absolute", dose="absolute", dose.units="cGy"),
										GY = y <- convert.DVH(x, volume="absolute", dose="absolute", dose.units="Gy"),
										"%" = y <- convert.DVH(x, volume="absolute", dose="relative"),
										y <- convert.DVH(x, volume="absolute", dose="absolute")
									)
									# may be flawed if bin widths variable!!!
									bin.widths <- median(diff(y@doses))
									start.i <- max(start.i, min(y))				
									end.i <- min(end.i, max(y))				
									if (units.i == "%") {
										result <- c(result,
											y@dose.rx * max(0, integrate(function(dose) {
												return(approx(y@doses, y@volumes*y@doses/bin.widths, dose, yleft=0, yright=0, ties=max)$y)
												}, start.i, end.i, stop.on.error=FALSE, abs.tol=0, rel.tol=100*.Machine$double.eps
											)$value) / y@rx.isodose
										)
									}
									else {
										result <- c(result,
											max(0, integrate(function(dose) {
												return(approx(y@doses, y@volumes*y@doses/bin.widths, dose, yleft=0, yright=0, ties=max)$y)
												}, start.i, end.i, stop.on.error=FALSE, abs.tol=0, rel.tol=100*.Machine$double.eps
											)$value)
										)
									}
									result.units <- c(result.units, paste(y@dose.units, "*cc", sep=""))
									next
								}
								warning("No method available to extract dose given differential doses")
								result <- c(result, NA)
								result.units <- c(result.units, NA)
							}
						)
						
					},
					{
						warning("Improper format '", input[count], "' (dose/volume specifier missing, e.g. 'V__' or 'D__')")
						result <- c(result, NA)
						result.units <- c(result.units, NA)
					}
				)
			}
			if (any(type6)) {
				for (count in which(type6)) {
					switch(type7[count],
						"==" = result[count] <- as.logical(result[count] == type8[count]),
						"=" = result[count] <- as.logical(result[count] == type8[count]),
						"<=" = result[count] <- as.logical(result[count] <= type8[count]),
						">=" = result[count] <- as.logical(result[count] >= type8[count]),
						"<" = result[count] <- as.logical(result[count] < type8[count]),
						">" = result[count] <- as.logical(result[count] > type8[count]),
						"!=" = result[count] <- as.logical(result[count] != type8[count])
					)	
					result.units[count] <- "logical"					
				}				
			}
		}
		else {
			return()
		}
		names(result) <- result.units
		return(result)
	}
)


setMethod("c", "DVH",
	function (x, ..., recursive = FALSE) {
		return(c(as(x, "DVH.list"), ..., recursive=recursive))
	}
)

setMethod("sum", "DVH",
	function (x, ..., na.rm = TRUE) {
		return(sum(as(x, "DVH.list"), ..., na.rm = na.rm))
	}
)

setMethod("print", "DVH",
	function (x, ...) {
		if (x@dose.type == "relative") {
			dose.type <- "%"
			dose.min <- x@dose.min * x@rx.isodose / x@dose.rx
			dose.max <- x@dose.max * x@rx.isodose / x@dose.rx
		}
		else {
			dose.type <- x@dose.units
			dose.min <- x@dose.min
			dose.max <- x@dose.max
		}
		print(paste("Structure: ", x@structure.name, " (", sprintf("%.*f", 1, x@structure.volume), "cc), Dose: ", sprintf("%.*f", 2, dose.min), "-", sprintf("%.*f", 2, dose.max), dose.type, " (", x@dose.rx, x@dose.units, " prescribed", if (x@rx.isodose != 100) {paste(" to ", x@rx.isodose, "% isodose line", sep="")}, "), DVH: ", x@type, ", Volume: ", x@volume.type, sep=""))
	}
)

setMethod("show", "DVH",
	function (object) {
		print(object)
	}
)

is.empty <- function (x) {
	if ((length(x@doses) < 1) & (x@structure.volume == 0)) {
		return(TRUE)
	}
	else {
		return(FALSE)
	}
}