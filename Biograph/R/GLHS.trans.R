GLHS.trans <-
function (names)
#Transition matrix of illness-death model (tmat)
# called by GLHS.IllnessDeath() to create J1,J2 and N (3 states)
# based on trans.illdeath function of mstate (Putter et al.)
{
    tmat <- matrix(NA, 3, 3)
    tmat[1, 2:3] <- 1:2
    tmat[2, 3] <- 3
    if (missing(names))
        names <- c("Job1", "Job2", "NoJob")
    else {
        if (length(names) != 3)
            stop("incorrect length of \"names\" argument")
    }
    dimnames(tmat) <- list(from = names, to = names)
    return(tmat)
}
