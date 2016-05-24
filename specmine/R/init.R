names = c("genefilter", "xcms", "impute", "MAIT")

for (n in 1:length(names)){
	if (!(names[n] %in% rownames(installed.packages()))){
		source("http://bioconductor.org/biocLite.R")
		biocLite(names[n])
	}
}

#if (!("mzmatch.R" %in% rownames(installed.packages()))){
#	source ("http://puma.ibls.gla.ac.uk/mzmatch.R/install_mzmatch.R")
#}
