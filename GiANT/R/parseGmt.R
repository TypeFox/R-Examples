parseGmt <-function(file){

	input <- readLines(file)
	input <- strsplit(input, split="\t")
	nms <- sapply(input,function(x){
			return(x[1])
		})
	geneSets <- lapply(input,function(x){
			return(x[3:length(x)])
		})
	names(geneSets) <- nms

	return(geneSets)
}
