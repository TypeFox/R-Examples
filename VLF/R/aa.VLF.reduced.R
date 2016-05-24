aa.VLF.reduced <-
function(NA.matrix, sCount, seqlength){
	extra_row <- 1:nrow(NA.matrix)
	convert <- matrix(c(extra_row, NA.matrix), ncol = seqlength + 3)
	list = c()
	for(n in 1:nrow(convert)){
		if(sCount[n] == 0){
			list = c(list,n)
		}
	}
	convert = convert[as.vector((-1)*list),]
	convert = convert[,-1]
	return(convert)
}
