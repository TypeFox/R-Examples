StepUp <- function(adjp.temp) {
    N <- length(adjp.temp) 
    running.min <- adjp.temp[N]
    for(i in N:1) {
        if(adjp.temp[i] > running.min) {
            adjp.temp[i] <- running.min
        } else {
            running.min <- adjp.temp[i]
        }
    }
    adjp.temp[adjp.temp > 1] <- 1
    return(adjp.temp)
}
