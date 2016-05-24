#' Estimates the stem volume 
#' 
#' \code{volume} uses one of the following methods (Smalian, Newton, Huber) to
#' approximate real stem volume. Users should remember they're just 
#' approximations and sample size provide more accurate results them using 
#' different methods.
#' 
#' @param trees a data frame or matrix in format described in dataset inventory
#' (more help \code{\link{inventory}})
#' @param method method used for estimation of the stem volume
#' @return a named vector of volumes, names are defined as same as in first 
#' column
#' @references
#' \url{http://wiki.awf.forst.uni-goettingen.de/wiki/index.php/Stem_volume}
#' @note Newton and Huber methods have small modifications for working just with
#' two mensures (lower and upper diameter). Both of them use mean instead of
#' real middle diameter.
#' @seealso \code{\link{ff}} \code{\link{sf}}
#' @examples
#' example_data <- data.frame(tree_number = 1, 
#'                            dhb = 5, 
#'                            total_height = 20, 
#'                            comercial_height = 15,
#'                            section_height = c(0,5,15),
#'                            section_diameter = 5
#'                            )
#' volume(example_data)
#' #
#' #
#' # A little more complex and common example
#' data(inventory)
#' volume_output <- volume(inventory)
#' summary(volume_output)
#' hist(volume_output)
#' @export
volume <- function(trees, method = 'smalian') {
	if ( !(is.data.frame(trees) || is.matrix(trees) ) ) 
		stop("trees need to be data.frame or matrix")
	methods <- c('smalian', 'newton', 'huber')
	if (is.na(pmatch(method, methods))) 
		stop("invalid method, type ?volume at R console for more information")
	group_trees <- by(trees, trees[, 1], function(x)x)
	r <- numeric(length(group_trees))
	
	for (i in 1:length(group_trees)) {
		v <- numeric(nrow(group_trees[[i]]) - 1)
		for (j in 2:nrow(group_trees[[i]])) {
			# [j] major section (extreme, large) and [j-1] minor section (small)
			l <- group_trees[[i]][j, 5] - group_trees[[i]][j - 1, 5]
			if (pmatch(method, methods) == 3) # huber method
				v[j - 1] <- l * pi * (mean(c(group_trees[[i]][j - 1, 6], 
										group_trees[[i]][j, 6])) / 2) ^ 2
			else if (pmatch(method, methods) == 2) # newton method
				v[j - 1] <- 1/12 * pi * l * (group_trees[[i]][j - 1, 6] ^ 2 + 
										group_trees[[i]][j - 1, 6] * group_trees[[i]][j, 6] + 
										group_trees[[i]][j, 6] ^ 2)
			else # smalian method
				v[j - 1] <- 1/8 * pi * l * (group_trees[[i]][j - 1, 6] ^ 2 + 
										group_trees[[i]][j, 6] ^ 2)
		}
		r[i] <- sum(v)
	}
	names(r) = names(group_trees)
	r
}


#' Factor form for the given volume
#' 
#' This function provide correction for basic volume estimation using cylinder
#' formulation \eqn{v = \frac{DBH^2}{4} \pi H}{ v = DBH^2/4 \pi H}. Factor form 
#' is given by taking ratio between real volume and apparent volume.
#' 
#' @param volume volume of a log, can be the output of \code{\link{volume}}
#' @param dbh diameter at breast height (1.3 meters from floor)
#' @param height commercial height, length of stem or whatever length of log you 
#' used in your estimations of cylinder volume
#' @return form factor ranging form 0-1 (numeric value)
#' @references
#' \url{http://wiki.awf.forst.uni-goettingen.de/wiki/index.php/Stem_shape}
#' @export
ff <- function(volume, dbh, height) {
	# height: comercial height
	volume / (dbh ** 2/4 * pi * height)
}
#' Stacking factor
#' 
#' Ratio between solid cubic meters of wood per stere cubic meter of stacked up 
#' wood.
#' 
#' @param volume real volume of the logs
#' @param height height of the stack (in meters)
#' @param length length of the stack (in meters)
#' @param depth depth of the stack (in meters)
#' @return stacking factor ranging from 0-1
#' @references
#' \url{www.eucalyptus.com.br/capitulos/ENG07.pdf}
#' @export
sf <- function(volume, height, length, depth) {
	volume / (height * length * depth)
}