requireAll <- function(packages) {
	.packages <- setdiff(packages, installed.packages()[,'Package'])
	if(length(.packages)>0) {
		suppressWarnings(rm(biocLite))
		source("http://bioconductor.org/biocLite.R")
		biocLite(.packages, dependencies=TRUE)
	}
	for(package in packages)
		suppressPackageStartupMessages(do.call(require, list(package)))
}
