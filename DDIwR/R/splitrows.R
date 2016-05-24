splitrows <- function(x, enter, y = 80, spacerep = "") {
    n <- x[1]
    command <- paste(toupper(x), collapse=", ")
    if (nchar(command) > y) {
        command <- precommand <- n
        for (ii in seq(2, length(x))) {
            if (nchar(precommand) > y) {
                precommand <- paste(toupper(x[ii]), ", ", sep="")
                command <- paste(command, ",", enter, spacerep, toupper(x[ii]), sep="")
            }
            else {
                precommand <- paste(precommand, toupper(x[ii]), sep=", ")
                command <- paste(command, toupper(x[ii]), sep=", ")
            }
            
        }
    }
    return(command)
}
