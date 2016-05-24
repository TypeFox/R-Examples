#
#	Wallraff procedure for comparing angular distances 
#
# Allows to compare the deviation from an angle of interest
# between several data sets. If the angle of interest is
# the mean direction, then it becomes a comparison of
# angular dispersion around the mean.
#
# In essence, it is a rank-based test (Wilcoxon-Mann-Withney
# or Kruskall-Wallis) on the angular distances from the angle
# of interest.
#
# (c) Copyright 2011 Jean-Olivier Irisson
#     GNU General Public License v3
#
#------------------------------------------------------------

# added drop=TRUE 20131106 Claudio

# Generic function
wallraff.test <- function(x, ...) {
	UseMethod("wallraff.test", x)
}

# Default method, for an angle vector and a grouping vector
wallraff.test.default <- function(x, group, ref=NULL, ...) {

	# get data name
	data.name <- paste(deparse(substitute(x)), "by", deparse(substitute(group)))

	# check arguments
	ok <- complete.cases(x, group)
	x <- x[ok]
	group <- group[ok,drop=TRUE]
	if (length(x)==0 | length(table(group)) < 2) {
		stop("No observations or no groups (at least after removing missing values)")
	}
	
	# make sure group is a factor 
	# if not, force it to keep the order in the original vector
	if (!is.factor(group)) {
		group <- factor(group, levels=unique(group))
		if (!is.null(ref)) {
			warning("\"group\" was converted into a factor.\n  The levels were kept in the order of the original vector:\n    ", paste(levels(group), collapse=", "), "\n  Please make sure the elements of \"ref\" match this order")			
		}
	}

	# convert data to the radians/trigonometric case
	if (is.circular(x)) {
		dc <- circularp(x)
		x <- conversion.circular(x, units="radians", zero=0, rotation="counter", modulo="2pi")
		attr(x, "class") <- attr(x, "circularp") <- NULL
	} else {
		dc <- list(type="angles", units="radians", template="none", modulo="asis", zero=0, rotation="counter")
	}
	
	if (is.null(ref)) {
		# if no reference angle is provided, use the mean angle
		xL <- split(x, group)
		ref <- sapply(xL, MeanCircularRad)
	} else {
		# when a reference angle is provided, if it is not of class circular, cast it as circular assuming that it is in the same reference as the original data
		if (!is.circular(ref)) {
			ref <- circular(ref, type=dc$type, units=dc$units, template=dc$template, modulo=dc$modulo, zero=dc$zero, rotation=dc$rotation)
		}
		
		# now that ref necessarily circular, convert it to the radians/trigonometric case
		ref <- conversion.circular(ref, units="radians", zero=0, rotation="counter", modulo="2pi")
		attr(ref, "class") <- attr(ref, "circularp") <- NULL
	}

	# compute concentration parameters and check assumptions
	result <- WallraffTestRad(x, group, ref)
	result$data.name <- data.name

	return(result)
}

# Method for a list
wallraff.test.list <- function(x, ref=NULL, ...) {
	# fecth or fill list names
	k <- length(x)
	if (is.null(names(x))) {
		names(x) <- 1:k
	}
	# get data name
	data.name <- paste(names(x), collapse=" and ")

	# convert into x and group
	ns <- lapply(x, length)
	group <- rep(names(x), times=ns)
	group <- factor(group, levels=unique(group))
	x <- do.call("c", x)
	# NB: unlist() removes the circular attributes here

	# call default method
	result <- wallraff.test.default(x, group, ref)
	result$data.name <- data.name

	return(result)
}

# Method for a formula
wallraff.test.formula <- function(formula, data, ref=NULL, ...) {
	# convert into x and group
	d <- model.frame(as.formula(formula), data)

	# get data name
	data.name <- paste(names(d), collapse=" by ")

	# call default method
	result <- wallraff.test.default(d[,1], d[,2], ref)
	result$data.name <- data.name

	return(result)
}


# Computation in the usual trigonometric space
WallraffTestRad <- function(x, group, ref) {
	
	# consolidate data
	if (length(ref) < nlevels(group)) {
		ref = rep(ref, nlevels(group))
	}
	d <- matrix(c(x, ref=ref[as.numeric(group)]), ncol=2)

	# compute angular distances = ranges
	dists <- apply(d, 1, function(X) RangeCircularRad(X[1:2], test=FALSE) )

	result <- kruskal.test(dists, group)
	# NB: kruskal.test with 2 groups is equivalent to wilcox.test with exact=FALSE and correct=FALSE

	result$method <- "Wallraff rank sum test of angular distance"
	
	return(result)
}
