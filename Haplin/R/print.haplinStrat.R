print.haplinStrat <- function(object){
##
##
##
cat('List of length ', length(object), ' containing results from a stratified haplin analysis.\nUse "summary" function to get detailed information.\n', sep = '') 

for(i in seq(along = object)){
	cat("\n#### Stratum = ", names(object)[i], " ####\n", sep = "")
	print(object[[i]])
}
return(invisible(object))

}
