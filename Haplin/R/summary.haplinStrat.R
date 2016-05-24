summary.haplinStrat <- function(object){
##
##
##
cat('List of length ', length(object), ' containing results from a stratified haplin analysis.\n', sep = '')

.summ <- lapply(object, summary)

for(i in seq(along = object)){
	cat("\n\n#### Stratum = ", names(object)[i], " ####---------\n", sep = "")
	print(.summ[[i]])
}
return(invisible(.summ))


}
