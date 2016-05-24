pickList <- function(data, ranges, nitems=10, nlists=1, replace=FALSE) {
	
	within <- function(x, a, b)
		(!is.na(x)) & (x >= a) & (x <= b)
	
	if (missing(ranges)) {
		vars <- sapply(data, is.numeric)
		ranges <- as.data.frame(apply(data[,vars], 2, function(x) range(na.omit(x))))
	}
	# allow a list of min/max, rather than a data.frame
	if (class(ranges)=='list') ranges <- as.data.frame(ranges)

	names <- colnames(ranges)
	OK <- rep(TRUE, nrow(data))
	for (col in names) {
		OK <- OK & within(data[,col], ranges[1,col], ranges[2,col])
	}
	data <- data[OK,]
	want <- sample(1:nrow(data), size=nitems*nlists, replace=replace)
	data <- data[want,]
	data <- cbind(list=(rep(1:nlists, each=nitems)), data)
	data
}
