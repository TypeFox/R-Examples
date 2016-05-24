`as.splineMatrix` <- 
function(matrix_argument,
	matrix_argument_name,
	content_area=NULL,
	year=NULL,
	sgp_object=NULL) {

		if (!class(matrix_argument) %in% c("matrix", "splineMatrix")) stop("Supplied object must be of class 'matrix' or 'splineMatrix'.")
		if (class(matrix_argument) == "splineMatrix" && validObject(matrix_argument, test=TRUE)==TRUE) return(matrix_argument)

		### Create relvant variables

			if (class(matrix_argument)=="matrix") tmp.matrix <- matrix_argument else tmp.matrix <- matrix_argument@.Data

			rn <- rownames(tmp.matrix)[-1]
			rn <- gsub("\"", "'", rn) 
	
			rn.knots <- strsplit(rn, "knots = ",)
			rn.knots <- unique(sapply(rn.knots, function(x) strsplit(x[2], ", Boundary.")))
			rn.knots2 <- sapply(rn.knots, function(x) strsplit(x, "knots_"))
			rn.knots2 <- sapply(rn.knots2, function(x) strsplit(x[2], "'"))

			rn.bounds <- strsplit(rn, "Boundary.knots = ",)
			rn.bounds <- sapply(rn.bounds, function(x) strsplit(x[2], ")"))
			rn.bounds <- unique(sapply(rn.bounds, function(x) x[1]))
			rn.bounds2 <- sapply(rn.bounds, function(x) strsplit(x, "boundaries_"))
			rn.bounds2 <- sapply(rn.bounds2, function(x) strsplit(x[2], "'"))

			tmp.last.grade <-  unlist(strsplit(matrix_argument_name, "_"))[2]
			tmp.num.prior <-  unlist(strsplit(matrix_argument_name, "_"))[3]

		### Matrix case ###

		if (class(matrix_argument)=="matrix") {

			if (is.null(sgp_object)) {
				stop("splineMatrix creation with an object of class 'matrix' requires that an sgp_object be supplied.")
			}

			if (is.SGP(sgp_object)) tmp.label <- "sgp_object@SGP$" else tmp.label <- "sgp_object$"


			### Knots

			knots <- list()
			for (i in seq_along(rn.knots)) {
				knots[[i]] <- eval(parse(text=paste(tmp.label, rn.knots[[i]], sep="")))
			}
			names(knots) <- paste("knots", sapply(rn.knots2, function(x) x[1]), sep="_")


			### Boundaries

			boundaries <- list()
			for (i in seq_along(rn.bounds)) {
				boundaries[[i]] <- eval(parse(text=paste(tmp.label, rn.bounds[i], sep="")))
			}
			names(boundaries) <- paste("boundaries", sapply(rn.bounds2, function(x) x[1]), sep="_")


			### Grade Progression

			grade_progression <- as.character(c(rev(sapply(rn.knots2, function(x) x[1])), tmp.last.grade))
			if (!is.numeric(type.convert(grade_progression))) {
				stop("Automatic conversion of older to newer version spline matrices is only available when grade progressions are integers. Please contact package maintainer for help on update of your splineMatrices.")
			}


			### Content Areas

			content_areas <- rep(unlist(strsplit(gsub("'|]]|\"", "", strsplit(rn, "\\[\\[|\\$")[[1]][2]), "[.]"))[1], length(grade_progression))


			### Time Lag

			time_lags <- as.integer(diff(type.convert(grade_progression)))


			### Time

			tmp.time <- unlist(strsplit(gsub("'|]]|\"", "", strsplit(rn, "\\[\\[|\\$")[[1]][2]), "[.]"))[2]
			if (!is.null(year) && tmp.time != year) {
				message("\tNOTE: Year from supplied splineMatrix does not equal year indicated in @SGP[['Coefficient_Matrices']]. Results will proceed based upon @SGP[['Coefficient_Matrices']]")
				tmp.time <- year
			}
			if (tmp.time == "BASELINE") {
				time <- rep("BASELINE", length(grade_progression))
			} else {
				time <- as.character(rev(yearIncrement(tmp.time, -cumsum(c(0, rev(time_lags))))))
			}


			### Version

			version <- list(SGP_Package_Version=as.character(packageVersion("SGP")), Date_Prepared=date())


			### Create new splineMatrix

			new("splineMatrix", 
				.Data=tmp.matrix, 
				Knots=knots, 
				Boundaries=boundaries, 
				Content_Areas=list(content_areas), 
				Grade_Progression=list(grade_progression),
				Time=list(time),
				Time_Lags=list(time_lags), 
				Version=version)

		} ### END if class(matrix_argument)=="matrix"

	
		### splineMatrix case ###

		if (class(matrix_argument)=="splineMatrix") {

			knots <- matrix_argument@Knots
			boundaries <- matrix_argument@Boundaries


			### Grade Progression

			if (.hasSlot(matrix_argument, "Grade_Progression")) {
				grade_progression <- as.character(matrix_argument@Grade_Progression[[1]])
			} else {
				grade_progression <- as.character(c(rev(sapply(rn.knots2, function(x) x[1])), tmp.last.grade))
				if (!is.numeric(type.convert(grade_progression))) {
					stop("Automatic conversion of older to newer version spline matrices is only available when grade progressions are integers. Please contact package maintainer for help on update of your splineMatrices.")
				}
			}


			### Content Areas

			if (.hasSlot(matrix_argument, "Content_Areas")) {
				content_areas <- as.character(matrix_argument@Content_Areas[[1]])
			} else {
				content_areas <- rep(unlist(strsplit(gsub("'|]]|\"", "", strsplit(rn, "\\[\\[|\\$")[[1]][2]), "[.]"))[1], length(grade_progression))
			}


			### Time Lag

			if (.hasSlot(matrix_argument, "Time_Lags")) {
				time_lags <- as.numeric(matrix_argument@Time_Lags[[1]])
			} else {
				time_lags <- as.numeric(diff(type.convert(grade_progression)))
			}


			### Time

			if (.hasSlot(matrix_argument, "Time") && matrix_argument@Version[['SGP_Package_Version']] > "1.0.6.0") {
				time <- as.character(matrix_argument@Time[[1]])
			} else {
				tmp.time <- unlist(strsplit(gsub("'|]]|\"", "", strsplit(rn, "\\[\\[|\\$")[[1]][2]), "[.]"))[2]
				if (!is.null(year) && tmp.time != year) {
					message("\tNOTE: Year from supplied splineMatrix does not equal year indicated in @SGP[['Coefficient_Matrices']]. Results will proceed based upon @SGP[['Coefficient_Matrices']]")
					tmp.time <- year
				}
				if (tmp.time == "BASELINE") {
					time <- rep("BASELINE", length(grade_progression))
				} else {
					time <- as.character(rev(yearIncrement(tmp.time, -cumsum(c(0, rev(time_lags))))))
				}
			}


			### Version

			version <- list(SGP_Package_Version=as.character(packageVersion("SGP")), Date_Prepared=date())


			### Create new splineMatrix

			new("splineMatrix", 
				.Data=tmp.matrix, 
				Knots=knots, 
				Boundaries=boundaries, 
				Content_Areas=list(content_areas), 
				Grade_Progression=list(grade_progression),
				Time=list(time),
				Time_Lags=list(time_lags), 
				Version=version)

		} ### END if class(matrix_argument)=="splineMatrix"
	
} ### END as.splineMatrix
