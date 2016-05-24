#' Sample from a \code{data.frame} according to a stratification variable
#' 
#' The \code{\link{stratified}} function samples from a
#' \code{\link{data.frame}} in which one of the columns can be used as a
#' "stratification" or "grouping" variable. The result is a new
#' \code{data.frame} with the specified number of samples from each group.
#' 
#' 
#' @param df The source \code{data.frame}.
#' @param group Your grouping variables. Generally, if you are using more than
#' one variable to create your "strata", you should list them in the order of
#' \emph{slowest} varying to \emph{quickest} varying. This can be a vector of
#' names or column indexes.
#' @param size The desired sample size. \itemize{ \item If \code{size} is a
#' value between \code{0} and \code{1} expressed as a decimal, size is set to
#' be proportional to the number of observations per group. \item If
#' \code{size} is a single positive integer, it will be assumed that you want
#' the same number of samples from each group. \item If \code{size} is a
#' vector, the function will check to see whether the length of the vector
#' matches the number of groups and use those specified values as the desired
#' sample sizes. The values in the vector should be in the same order as you
#' would get if you tabulated the grouping variable (usually alphabetic order);
#' alternatively, you can name each value to ensure it is properly matched. }
#' @param select A named list containing levels from the "group" variables in
#' which you are interested. The list names must be present as variable names
#' for the input \code{data.frame}.
#' @param seed The seed that you want to use (using \code{\link{set.seed}}
#' \emph{within} the function, if any. Defaults to \code{NULL}.
#' @param \dots Further arguments to be passed to the \code{\link{sample}}
#' function.
#' @note \emph{Slightly different sizes than requested}
#' 
#' Because of how computers deal with floating-point arithmetic, and because R
#' uses a "round to even" approach, the size per strata that results when
#' specifying a proportionate sample may be slightly higher or lower per strata
#' than you might have expected.
#' 
#' \emph{"Seed" argument}
#' 
#' This is different from using \code{\link{set.seed}} before using the
#' function. Setting a seed using this argument is equivalent to using
#' \code{set.seed} each time that you go to take a sample from a different
#' group (in other words, the same seed is used for each group).
#' 
#' The inclusion of a \code{seed} argument is mostly a matter of convenience,
#' to be able to have a single seed with which the samples can be verified
#' later. However, by using the \code{seed} argument, the same seed is used to
#' sample from each group. This may be a problem if there are many groups that
#' have the same number of observations, since it means that the same
#' observation number will be selected from each of those grops. For instance,
#' if group "AA" and "DD" both had the same number of observations (say, 5) and
#' you were sampling 3 cases using a seed of 1, the second, fifth, and fourth
#' observation would be taken from each of those groups. To avoid this, you can
#' set the seed using \code{set.seed} \emph{before} you run the
#' \code{\link{stratified}} function.
#' 
#' As a user, you need to weigh the benefits and drawbacks of setting the seed
#' \emph{before} running the function as opposed to setting the seed
#' \emph{with} the function.  Setting the seed before would be useful if there
#' are several groups with the same number of observations; however, in the
#' slim chance that you need to verify the samples manually, you \emph{may} run
#' into problems.
#' @author Ananda Mahto
#' @references The evolution of this function can be traced at the following
#' links. The version in this package is entirely reworked and does not require
#' an additional package to be loaded.
#' 
#' \itemize{ \item
#' \url{http://news.mrdwab.com/2011/05/15/stratified-random-sampling-in-r-beta/}
#' \item
#' \url{http://news.mrdwab.com/2011/05/20/stratified-random-sampling-in-r-from-a-data-frame/}
#' \item \url{http://stackoverflow.com/a/9714207/1270695} }
#' @examples
#' 
#' # Generate a couple of sample data.frames to play with
#' set.seed(1)
#' dat1 <- data.frame(ID = 1:100,
#'               A = sample(c("AA", "BB", "CC", "DD", "EE"), 100, replace = TRUE),
#'               B = rnorm(100), C = abs(round(rnorm(100), digits=1)),
#'               D = sample(c("CA", "NY", "TX"), 100, replace = TRUE),
#'               E = sample(c("M", "F"), 100, replace = TRUE))
#' dat2 <- data.frame(ID = 1:20,
#'               A = c(rep("AA", 5), rep("BB", 10),
#'                     rep("CC", 3), rep("DD", 2)))
#' # What do the data look like in general?
#' summary(dat1)
#' summary(dat2)
#' 
#' # Let's take a 10% sample from all -A- groups in dat1, seed = 1
#' stratified(dat1, "A", .1, seed = 1)
#' 
#' # Let's take a 10% sample from only "AA" and "BB" groups from -A- in dat1, seed = 1
#' stratified(dat1, "A", .1, select = list(A = c("AA", "BB")), seed = 1)
#' 
#' # Let's take 5 samples from all -D- groups in dat1,
#' #   seed = 1, specified by column number
#' stratified(dat1, group = 5, size = 5, seed = 1)
#' 
#' # Let's take a sample from all -A- groups in dat1, seed = 1,
#' #   where we specify the number wanted from each group
#' stratified(dat1, "A", size = c(3, 5, 4, 5, 2), seed = 1)
#' 
#' # Use a two-column strata: -E- and -D-
#' #   -E- varies more slowly, so it is better to put that first
#' stratified(dat1, c("E", "D"), size = .15, seed = 1)
#' 
#' # Use a two-column strata (-E- and -D-) but only interested in
#' #   cases where -E- == "M"
#' stratified(dat1, c("E", "D"), .15, select = list(E = "M"), seed = 1)
#' 
#' ## As above, but where -E- == "M" and -D- == "CA" or "TX"
#' stratified(dat1, c("E", "D"), .15,
#'      select = list(E = "M", D = c("CA", "TX")), seed = 1)
#' 
#' # Use a three-column strata: -E-, -D-, and -A-
#' s.out <- stratified(dat1, c("E", "D", "A"), size = 2, seed = 1)
#' 
#' list(head(s.out), tail(s.out))
#' 
#' # How many samples were taken from each strata?
#' table(interaction(s.out[c("E", "D", "A")]))
#' 
#' # Can we verify the message about group sizes?
#' names(which(table(interaction(dat1[c("E", "D", "A")])) < 2))
#' 
#' names(which(table(interaction(s.out[c("E", "D", "A")])) < 2))
#' \dontshow{rm(dat1, dat2, s.out)}
#' 
#' @export stratified
stratified <- function(df, group, size, select = NULL, seed = NULL, ...) {
	if (is.null(select)){
		df = df
	} else {
		if (is.null(names(select))) stop("'select' must be a named list")
		if (!all(names(select) %in% names(df))) stop("Please verify your 'select' argument")
		temp <- sapply(names(select), function(x) df[[x]] %in% select[[x]])
		df <- df[rowSums(temp) == length(select), ]
	}

	df.interaction <- interaction(df[group], drop = TRUE)
	df.table <- table(df.interaction)
	df.split <- split(df, df.interaction)
	if (length(size) > 1) {
		if (length(size) != length(df.split))
			stop("Number of groups is ", length(df.split)," but number of sizes supplied is ", length(size))
	if (is.null(names(size))) {
		n <- setNames(size, names(df.split))
		message(sQuote("size"), " vector entered as:\n\nsize = structure(c(", 
			paste(n, collapse = ", "), "),\n.Names = c(",
			paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
	} else {
		ifelse(all(names(size) %in% names(df.split)), 
			n <- size[names(df.split)], 
			stop("Named vector supplied with names ", 
				paste(names(size), collapse = ", "),
					"\n but the names for the group levels are ", 
					paste(names(df.split), collapse = ", ")))
	} 
	} else if (size < 1) {
		n <- round(df.table * size, digits = 0)
	} else if (size >= 1) {
		if (all(df.table >= size)) {
		n <- setNames(rep(size, length.out = length(df.split)),
			names(df.split))
	} else {
		message(
		"Some groups\n---", 
		paste(names(df.table[df.table < size]), collapse = ", "),
		"---\ncontain fewer observations",
		" than desired number of samples.\n",
		"All observations have been returned from those groups.")
		n <- c(sapply(df.table[df.table >= size], function(x) x = size),
				df.table[df.table < size])
		}
	}
	
	seedme <- ifelse(is.null(seed), "No", "Yes")

	temp <- switch(
	seedme,
	No = { temp <- lapply(
		names(df.split), 
		function(x) df.split[[x]][sample(df.table[x], 
							n[x], ...), ]) },
	Yes = { temp <- lapply(
		names(df.split),
		function(x) { set.seed(seed)
			df.split[[x]][sample(df.table[x], 
				n[x], ...), ] })})
	if (!is.null(seed)) {
		rm(.Random.seed, envir=.GlobalEnv) # "resets" the seed
	}

	do.call("rbind", temp)
}