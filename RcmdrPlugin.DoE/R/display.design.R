display.design <- function(){
    command <- paste("print(",ActiveDataSet(),")")
    doItAndPrint(command)
}