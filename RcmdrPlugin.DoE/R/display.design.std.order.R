display.design.std.order <- function(){
    command <- paste("print(",ActiveDataSet(),", std.order=TRUE)")
    doItAndPrint(command)
}