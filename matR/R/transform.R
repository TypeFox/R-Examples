
#-----------------------------------------------------------------------------------------
#  data transformations of biom objects.
#
#  for matrix transformation, we have a set of "primitives".
#  for biom transformation, we have a transform() method that applies them.
#-----------------------------------------------------------------------------------------

transform.biom <- function (`_data`, ...) {
	x <- `_data`
	li <- list (...)
	if (is.null (names (li))) {
		f.list <- li
		a.list <- replicate (length (li), NULL)
	} else {
		f.list <- a.list <- list()
		length (f.list) <- length (a.list) <- length (li)

		f.list [names (li) == ""] <- li [names (li) == ""]
		a.list [names (li) == ""] <- replicate (sum (names (li) == ""), NULL)

		f.list [names (li) != ""] <- mget (names (li) [names (li) != ""], inherits=TRUE)
		a.list [names (li) != ""] <- li [names (li) != ""]
		}

####  this is:  list(x, list(list of functions, list of argument lists))

	ll <- append(
			list (as.matrix (x, expand=TRUE)),
			mapply (list, f.list, a.list, SIMPLIFY=FALSE))
	pass.to <- function (x, funcWithArgs) {

####  there was an issue here:
####  args being passed in a list capsule should not be.
####  but it is resolved now?

		do.call (funcWithArgs [[1]], append (list (x), funcWithArgs [[2]]))
		}
	data1 <- Reduce (pass.to, ll)

	y <- x [rownames(x) %in% rownames (data1), colnames (x) %in% colnames (data1)]
	y$data <- data1
	y$sparse <- NULL
	y$generated_by <- tagline()
	y$date <- strftime(Sys.time())
	y$id <- paste0 ("derived with ", deparse (match.call(), width.cutoff=500))
	y
	}

t_NA2Zero <- function (x, ...) {
	x [is.na (x)] <- 0;
	x
	}

t_Threshold <- function (x, entry.min=2, row.min=2, col.min=2) {
	x [x < entry.min] <- 0
	x <- x [rowSums (x) >= row.min, ]
	x <- x [, colSums (x) >= col.min]
	x
	}

t_Log <- function (x, ...) {
	log2 (1 + x)
	}

t_ColCenter <- function (x, ...) {
	x - colMeans(x) [col(x)]
	}

t_ColScale <- function (x, ...) {
	sigma <- apply (x, 2, sd)
	sigma [sigma == 0] <- 1
	x / sigma [col(x)]
	}

t_DENorm <- function (x, DEparam, ...) { x }
