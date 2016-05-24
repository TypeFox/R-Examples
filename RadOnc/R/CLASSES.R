############################################
## CLASS: DVH
############################################
setClass("DVH",
	representation(
		patient = "character",
		ID = "character",
		structure.name = "character",
		structure.volume = "numeric",
		type = "character",
		dose.max = "numeric",
		dose.min = "numeric",
		dose.mean = "numeric",
		dose.median = "numeric",
		dose.mode = "numeric",
		dose.STD = "numeric",
		conf.index = "numeric",
		equiv.sphere = "numeric",
		gradient = "numeric",
		plan.sum = "logical",
		dose.rx = "numeric",
		dose.fx = "numeric",
		rx.isodose = "numeric",
		doses = "numeric",
		dose.type = "character",
		dose.units = "character",
		volumes = "numeric",
		volume.type = "character"
	),
	prototype(
		patient = character(),
		ID = character(),
		structure.name = character(),
		structure.volume = numeric(),
		type = character(),
		dose.max = numeric(),
		dose.min = numeric(),
		dose.mean = numeric(),
		dose.median = numeric(),
		dose.mode = numeric(),
		dose.STD = numeric(),
		conf.index = numeric(),
		equiv.sphere = numeric(),
		gradient = numeric(),
		plan.sum = logical(),
		dose.rx = numeric(),
		dose.fx = numeric(),
		rx.isodose = numeric(),
		doses = numeric(),
		dose.type = character(),
		dose.units = character(),
		volumes = numeric(),
		volume.type = character()
	)
)

setMethod("initialize",
	"DVH",
	function (.Object,
		patient = "",
		ID = "",
		structure.name = "",
		structure.volume = numeric(),
		type = c("cumulative", "differential"),
		dose.max = NA,
		dose.min = NA,
		dose.mean = numeric(),
		dose.median = numeric(),
		dose.mode = numeric(),
		dose.STD = numeric(),
		conf.index = numeric(),
		equiv.sphere = numeric(),
		gradient = numeric(),
		plan.sum = FALSE,
		dose.rx = NA,
		dose.fx = numeric(),
		rx.isodose = 100,
		doses = numeric(),
		dose.type = c("absolute", "relative"),
		dose.units = c("cGy", "Gy"),
		volumes = numeric(),
		volume.type = c("relative", "absolute"),
		...
	) {
		.Object@patient <- as.character(patient)
		.Object@ID <- as.character(ID)
		.Object@structure.name <- as.character(structure.name)
		.Object@structure.volume <- max(0, as.numeric(structure.volume), na.rm=TRUE)
		.Object@type <- match.arg(type)
		if (length(doses) > 0) {
			if (is.na(dose.max)) dose.max <- range(doses)[2]
			if (is.na(dose.min)) dose.min <- range(doses)[1]			
			.Object@doses <- doses	
		}
		else {
			.Object@doses <- numeric()
		}
		.Object@dose.max <- max(0, dose.max, na.rm=TRUE)		
		.Object@dose.min <- max(0, dose.min, na.rm=TRUE)
		.Object@dose.mean <- max(0, dose.mean, na.rm=TRUE)
		.Object@dose.median <- max(0, dose.median, na.rm=TRUE)
		.Object@dose.mode <- max(0, dose.mode, na.rm=TRUE)
		.Object@dose.STD <- max(0, dose.STD, na.rm=TRUE)
		.Object@conf.index <- max(0, conf.index, na.rm=TRUE)
		.Object@equiv.sphere <- max(0, equiv.sphere, na.rm=TRUE)
		.Object@gradient <- max(0, gradient, na.rm=TRUE)
		.Object@plan.sum <- plan.sum
		.Object@dose.rx <- max(0, dose.rx, na.rm=FALSE)
		.Object@dose.fx <- max(0, dose.fx, na.rm=TRUE)
		if (is.na(rx.isodose)) rx.isodose <- 100
		.Object@rx.isodose <- max(0, rx.isodose, na.rm=TRUE)
		.Object@dose.type <- match.arg(dose.type)
		.Object@dose.units <- match.arg(dose.units)
		.Object@volume.type <- match.arg(volume.type)
		if (length(volumes) > 0) {
			.Object@volumes <- as.numeric(volumes)
			if ((.Object@structure.volume <= 0) || (is.na(.Object@structure.volume))) {
				if ((.Object@volume.type == "absolute") & (.Object@type == "differential")) {
					.Object@structure.volume <- sum(.Object@volumes)
				}
				else if ((.Object@volume.type == "absolute") & (.Object@type == "cumulative")) {
					.Object@structure.volume <- sum(-diff(.Object@volumes))
				}
				else {
					warning("Cannot infer missing structure volume from relative volumetric data")
				}
				return(.Object)
			}
			if ((.Object@structure.volume < max(.Object@volumes, na.rm=TRUE)) & (.Object@volume.type == "absolute")) {
				.Object@structure.volume <- max(.Object@volumes, na.rm=TRUE)
			}
		}
		else {
			.Object@volumes <- numeric()	
		}
		return(.Object)
	}
)


setValidity("DVH",
	function(object) {
		if (length(object@doses) != length(object@volumes)) return(FALSE)
		if (length(object@doses) == 0) return(TRUE)
		if (any(is.na(object@doses))) return(FALSE)
		if (!is.na(object@dose.rx) & (object@dose.rx <= 0)) return(FALSE)
		if (!object@plan.sum & (is.na(object@rx.isodose) | (object@rx.isodose <= 0))) return(FALSE) 
		if (object@dose.min > object@dose.max) return(FALSE)
		if ((object@dose.mean > object@dose.max) | ((object@dose.mean < object@dose.min) & (object@dose.mean > 0))) return(FALSE)
		if (!identical(order(object@doses, decreasing=FALSE), 1:length(object@doses))) return(FALSE)
		if ((object@dose.type == "relative") & (range(object@doses, na.rm=TRUE)[2] > 250)) return(FALSE)
		if (any(is.na(object@volumes))) return(FALSE)		
		# ENSURE RELATIVE DVH VOLUMES ARE ON SCALE UP TO 100% MAXIMUM (VOLUME SHOULD NEVER BE >100%)		if ((object@volume.type == "relative") & (max(object@volumes, na.rm=TRUE) > 100.00000000001)) return(FALSE)	
		# ENSURE STRUCTURE VOLUME IS SUFFICIENTLY LARGE TO ENCOMPASS ALL LISTED DVH INFORMATION	
		if ((object@volume.type == "absolute") & (object@structure.volume < max(object@volumes, na.rm=TRUE))) return(FALSE)
		# ENSURE CUMULATIVE DOSE HISTOGRAM HAS APPROPRIATE DATA (DOSE RANGE MUST START AT 0)	
		# if ((object@type == "cumulative") & (object@doses[1] != 0)) return(FALSE)
		# ENSURE STRUCTURE VOLUME AND DVH VOLUME DATA ARE EQUIVALENT (TOLERANCE=0.1%)
		if (!grepl("(mean|median)[(].*[)]", object@structure.name)) {
			if ((object@type == "differential") & (object@volume.type == "relative") & (abs(sum(object@volumes, na.rm=TRUE) - 100) > 0.1)) return(FALSE)
			if ((object@type == "differential") & (object@volume.type == "absolute") & (abs(sum(object@volumes, na.rm=TRUE) - object@structure.volume) / object@structure.volume > 0.001)) return(FALSE)
		}
		return(TRUE)
	}
)


############################################
## CLASS: zDVH
############################################
setClass("zDVH",
	contains="DVH"
)


setMethod("initialize",
	"zDVH",
	function (.Object,
		patient = "",
		ID = "",
		structure.name = "",
		structure.volume = numeric(),
		type = c("cumulative", "differential"),
		dose.max = NA,
		dose.min = NA,
		dose.mean = numeric(),
		dose.median = numeric(),
		dose.mode = numeric(),
		dose.STD = numeric(),
		conf.index = numeric(),
		equiv.sphere = numeric(),
		gradient = numeric(),
		plan.sum = FALSE,
		dose.rx = NA,
		dose.fx = numeric(),
		rx.isodose = 100,
		doses = numeric(),
		dose.type = c("absolute", "relative"),
		dose.units = c("cGy", "Gy"),
		volumes = matrix(nrow=0, ncol=0),
		volume.type = c("relative", "absolute"),
		...
	) {
		.Object@patient <- as.character(patient)
		.Object@ID <- as.character(ID)
		.Object@structure.name <- as.character(structure.name)
		.Object@structure.volume <- max(0, as.numeric(structure.volume), na.rm=TRUE)
		.Object@type <- match.arg(type)
		if (length(doses) > 0) {
			if (is.na(dose.max)) dose.max <- range(doses)[2]
			if (is.na(dose.min)) dose.min <- range(doses)[1]			
			.Object@doses <- doses	
		}
		else {
			.Object@doses <- numeric()
		}
		.Object@dose.max <- max(0, dose.max, na.rm=TRUE)		
		.Object@dose.min <- max(0, dose.min, na.rm=TRUE)
		.Object@dose.mean <- max(0, dose.mean, na.rm=TRUE)
		.Object@dose.median <- max(0, dose.median, na.rm=TRUE)
		.Object@dose.mode <- max(0, dose.mode, na.rm=TRUE)
		.Object@dose.STD <- max(0, dose.STD, na.rm=TRUE)
		.Object@conf.index <- max(0, conf.index, na.rm=TRUE)
		.Object@equiv.sphere <- max(0, equiv.sphere, na.rm=TRUE)
		.Object@gradient <- max(0, gradient, na.rm=TRUE)
		.Object@plan.sum <- plan.sum
		.Object@dose.rx <- max(0, dose.rx, na.rm=FALSE)
		.Object@dose.fx <- max(0, dose.fx, na.rm=TRUE)
		if (is.na(rx.isodose)) rx.isodose <- 100
		.Object@rx.isodose <- max(0, rx.isodose, na.rm=TRUE)
		.Object@dose.type <- match.arg(dose.type)
		.Object@dose.units <- match.arg(dose.units)
		.Object@volume.type <- match.arg(volume.type)
		class(volumes) <- c("numeric", "matrix")
		if (is.null(colnames(volumes))) {
			colnames(volumes) <- 1:(dim(volumes)[2])
		}
		.Object@volumes <- volumes
		if (length(volumes) > 0) {
			if ((length(.Object@structure.volume) < 1) || (is.na(.Object@structure.volume)) || (.Object@structure.volume < max(.Object@volumes, na.rm=TRUE))) {
				if ((.Object@volume.type == "absolute") & (.Object@type == "differential")) {
					.Object@structure.volume <- sum(.Object@volumes)
				}
				else if ((.Object@volume.type == "absolute") & (.Object@type == "cumulative")) {
					.Object@structure.volume <- sum(-apply(.Object@volumes, 2, diff))
				}
				else {
					warning("Cannot infer missing structure volume from relative volumetric data")
				}
				return(.Object)
			}
		}
		return(.Object)
	}
)


setValidity("zDVH",
	function(object) {
		if (!is.matrix(object@volumes)) return(FALSE)
		if (length(object@doses) != dim(object@volumes)[1]) return(FALSE)
		if (is.null(colnames(object@volumes))) return(FALSE)
		if (length(object@doses) == 0) return(TRUE)
#		if (length(object@doses) < 2) return(FALSE)
		if (any(is.na(object@doses))) return(FALSE)
		if (!is.na(object@dose.rx) & (object@dose.rx <= 0)) return(FALSE)
		if (!object@plan.sum & (is.na(object@rx.isodose) | (object@rx.isodose <= 0))) return(FALSE) 
		if (object@dose.min > object@dose.max) return(FALSE)
		if ((object@dose.mean > object@dose.max) | (object@dose.mean < object@dose.min)) return(FALSE)
		if (!identical(order(object@doses, decreasing=FALSE), 1:length(object@doses))) return(FALSE)
		if ((object@dose.type == "relative") & (range(object@doses, na.rm=TRUE)[2] > 250)) return(FALSE)		
		if (any(is.na(object@volumes))) return(FALSE)		
		# ENSURE RELATIVE DVH VOLUMES ARE ON SCALE UP TO 100% MAXIMUM (VOLUME SHOULD NEVER BE >100%)		if ((object@volume.type == "relative") & (max(object@volumes, na.rm=TRUE) > 100.00000000001)) return(FALSE)	
		# ENSURE STRUCTURE VOLUME IS SUFFICIENTLY LARGE TO ENCOMPASS ALL LISTED DVH INFORMATION	
		if ((object@volume.type == "absolute") & (object@structure.volume < max(object@volumes, na.rm=TRUE))) return(FALSE)
		# ENSURE CUMULATIVE DOSE HISTOGRAM HAS APPROPRIATE DATA (DOSE RANGE MUST START AT 0)	
		# if ((object@type == "cumulative") & (object@doses[1] != 0)) return(FALSE)
		# ENSURE STRUCTURE VOLUME AND DVH VOLUME DATA ARE EQUIVALENT (TOLERANCE=0.1%)
		if (!grepl("(mean|median)[(].*[)]", object@structure.name)) {
			if ((object@type == "differential") & (object@volume.type == "relative") & (abs(sum(object@volumes, na.rm=TRUE) - 100) > 0.1)) return(FALSE)
			if ((object@type == "differential") & (object@volume.type == "absolute") & (abs(sum(object@volumes, na.rm=TRUE) - object@structure.volume) / object@structure.volume > 0.001)) return(FALSE)
		}
		return(TRUE)
	}
)


############################################
## CLASS: DVH.list
############################################
setClass("DVH.list",
	representation(
		structures = "list"
	),
	prototype(
		structures = list()
	)
)


setMethod("initialize",
	"DVH.list",
	function (.Object,
		structures = list(),
		...
	) {
		if ((length(structures) == 1) & (class(structures) %in% c("DVH", "zDVH"))) {
			structures <- list(structures)
		}
		DVHs <- which(unlist(lapply(structures, class)) %in% c("DVH", "zDVH"))
		if (length(DVHs) >= 1) {
			.Object@structures <- structures[DVHs]
		}
		else {
			.Object@structures <- list()
		}
		return(.Object)
	}
)

setValidity("DVH.list",
	function(object) {
		if (length(object) == 0) return(TRUE)
		if (!all(unlist(lapply(object, class)) %in% c("DVH", "zDVH"))) return(FALSE)
		return(TRUE)
	}
)


############################################
## CLASS: structure3D
############################################
setClass("structure3D",
	representation(
		name = "character",
		volume = "numeric",
		volume.units = "character",
		coordinate.units = "character",
		vertices = "matrix",
		origin = "numeric",
		triangles = "matrix",
		closed.polys = "matrix",
		DVH = "DVH"
	),
	prototype(
		name = character(),
		volume = numeric(),
		volume.units = character(),
		coordinate.units = character(),
		vertices = matrix(),
		origin = numeric(),
		triangles = matrix(),
		closed.polys = matrix(),
		DVH = new("DVH")
	)
)

setMethod("initialize",
	"structure3D",
	function(.Object,
		name = "",
		volume = NULL,
		volume.units = c("cc"),
		coordinate.units = c("cm", "mm"),
		vertices = matrix(nrow=0, ncol=3),
		origin = NULL,
		triangles = matrix(nrow=3, ncol=0),
		closed.polys = matrix(nrow=0, ncol=3),
		DVH = new("DVH")
	) {
		.Object@name <- as.character(name)
		if (is.null(volume)) {
			.Object@volume <- 0
			# calculate volume of structure3D
			# .Object@volume <- as.numeric(volume)
		}
		else {
			.Object@volume <- as.numeric(volume)		
		}		
		if (is.null(vertices)) {
			vertices <- matrix(nrow=0, ncol=3)
		}
		if (is.null(closed.polys)) {
			closed.polys <- matrix(nrow=0, ncol=3)
		}
		if (is.null(origin)) {
			if (dim(vertices)[1] <= 1) {
				origin <- as.numeric(vertices)
			}
			else {
				origin <- apply(vertices, 2, mean)
			}
		}
		if (length(origin) != 3) {
			.Object@origin <- c(0, 0, 0)
		}
		else {
			.Object@origin <- origin
		}
		volume.units <- match.arg(volume.units)
		.Object@volume.units <- as.character(volume.units)
		coordinate.units <- match.arg(coordinate.units)
		.Object@coordinate.units <- as.character(coordinate.units)
		.Object@vertices <- as.matrix(vertices)
		.Object@triangles <- as.matrix(triangles)
		.Object@closed.polys <- as.matrix(closed.polys)
		return(.Object)
	}
)

setValidity("structure3D",
	function(object) {
		if (object@volume < 0) return(FALSE)
		if (!is.matrix(object@vertices)) return(FALSE)
		if (!is.matrix(object@triangles)) return(FALSE)
		if (dim(object@vertices)[2] != 3) return(FALSE)
		if (dim(object@triangles)[1] != 3) return(FALSE)
		if (length(object@origin) != 3) return(FALSE)
		if ((dim(object@triangles)[2] > 0) & (dim(object@vertices)[1] == 0)) return(FALSE)
#		if ((dim(object@vertices)[1] > 0) & (dim(object@triangles)[2] == 0)) return(FALSE)
		if ((dim(object@vertices)[1] > 0) & (dim(object@triangles)[2] > 0)) {
			range.triangles <- suppressWarnings(range(object@triangles))
			if (range.triangles[1] < 1) return(FALSE)
			if (range.triangles[2] > dim(object@vertices)[1]) return(FALSE)			
		}
		return(validObject(object@DVH))
	}
)

############################################
## CLASS: structure.list
############################################
setClass("structure.list",
	representation(
		structures = "list"
	),
	prototype(
		structures = list()
	)
)


setMethod("initialize",
	"structure.list",
	function (.Object,
		structures = list(),
		...
	) {
		if ((length(structures) == 1) & (class(structures) == "structure3D")) {
			structures <- list(structures)
		}
		structs <- which(unlist(lapply(structures, class)) == "structure3D")
		if (length(structs) >= 1) {
			.Object@structures <- structures[structs]
		}
		else {
			.Object@structures <- list()
		}
		return(.Object)
	}
)


setValidity("structure.list",
	function(object) {
		if (!is.list(object)) return(FALSE)
		if (length(object) == 0) return(TRUE)
		if (!all(unlist(lapply(structures, class)) == "structure3D")) return(FALSE)
		return(TRUE)
	}
)


############################################
## CLASS: RTdata
############################################
setClass("RTdata",
	representation(
		name = "character",
		CT = "array",
		dose = "array",
		structures = "structure.list"
	),
	prototype(
		name = character(),
		CT = array(dim=c(0,0,0)),
		dose = array(dim=c(0,0,0)),
		structures = new("structure.list")
	)
)

setMethod("initialize",
	"RTdata",
	function (.Object,
		name = character(),
		CT = array(dim=c(0,0,0)),
		dose = array(dim=c(0,0,0)),
		dose.units = character(),
		structures = new("structure.list"),
		...
	) {
		.Object@name <- name
		if (class(CT) != "array") {
			CT <- array(dim=c(0,0,0))
		}
		.Object@CT <- CT
		if (class(dose) != "array") {
			dose <- array(dim=c(0,0,0))
		}
		.Object@dose <- dose
		attr(.Object@dose, "dose.units") <- character()
		if (class(structures) == "structure.list") {
			.Object@structures <- structures
		}
		else {
			.Object@structures <- new("structure.list")
		}
		return(.Object)
	}
)

setValidity("RTdata",
	function(object) {
		if (class(object@CT) != "array") return(FALSE)
		if (length(dim(object@CT)) != 3) return(FALSE)
		if (class(object@dose) != "array") return(FALSE)
		if (length(dim(object@dose)) != 3) return(FALSE)
		if (is.null(attr(object@dose, "dose.units"))) return(FALSE)
		if (any(dim(object@dose) > 0) & (!attr(object@dose, "dose.units") %in% c("cGy", "Gy"))) return(FALSE)
		if (class(object@structures) != "structure.list") return(FALSE)
		return(TRUE)
	}
)


############################################
## INITIALIZE COMMON GENERIC FUNCTIONS
############################################
setGeneric("print",
	print
)

setGeneric("as.list",
	as.list
)

setGeneric("lapply",
	lapply
)

setGeneric("rev",
	rev
)
