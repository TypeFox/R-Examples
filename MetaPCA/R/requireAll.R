requireAll <- function(packages) {
	.packages <- setdiff(packages, installed.packages()[,'Package'])
	if(length(.packages)>0) {
		source("http://bioconductor.org/biocLite.R")
		biocLite(.packages, dependencies=TRUE)
	}
	for(package in packages)
		do.call(require, list(package))
}
