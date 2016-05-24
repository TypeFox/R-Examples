summarize.design <- function(){
    command <- paste("summary(",ActiveDataSet(),", brief = TRUE)")
    doItAndPrint(command)
}