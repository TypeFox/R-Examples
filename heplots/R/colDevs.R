colDevs <- function(x, group, center=mean,  ...) {

	if (!inherits(x, c("data.frame", "matrix")))
		stop("Argument 'x' must be a data.frame or matrix")

	if (missing(group)) {
		group <- factor(rep(1, nrow(x)))
	}
	
	if (!is.factor(group)) {
	warning(deparse(substitute(group)), " coerced to factor.")
	group <- as.factor(group)
	}
	
	nlev <- nlevels(group)
	lev <- levels(group)

	nums <- sapply(x, is.numeric)
	if (!all(nums)) {
		warning("Ignoring ", sum(!nums), " non-numeric column(s)")
		x <- x[, nums, drop=FALSE]
	}
	
	mat <- matrix(0, nrow(x), ncol(x), dimnames=dimnames(x))
	x <- as.matrix(x)
	for (i in 1:nlev) {
		rows <- which(group==lev[i])
		ctr <- apply( x[rows,], 2, center)
		mat[rows, ] <- sweep( x[rows, ], 2, ctr)
	}
	mat
}


TESTME <- FALSE
if(TESTME) {

# multivariate versions of Levene's test
Species <- iris$Species
irisdev <- colDevs(iris[,1:4], Species, mean)
irisdev.mod <- lm( abs(irisdev) ~ Species)
car::Anova(irisdev.mod)
heplot(irisdev.mod)
pairs(irisdev.mod)

irisdev <- colDevs(iris[,1:4], Species, median)
irisdev.mod <- lm( abs(irisdev) ~ Species)
car::Anova(irisdev.mod)

irisdev <- colDevs(iris[,1:4], Species, function(x) mean(x, trim=0.25))
irisdev.mod <- lm( abs(irisdev) ~ Species)
car::Anova(irisdev.mod)

}
