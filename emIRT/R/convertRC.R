convertRC <- function(.rc, type = "binIRT") {
    if (class(.rc) != "rollcall") {
        stop("Input must be a rollcall object.")
    }
                                        # output
                                        # yea : 1
                                        # nay : -1
                                        # missing : 0
                                        # notInLegis : 9
    ## copy old over
    newrc <- .rc
    ## recode votes
    newrc$votes[, ] <- 0
    for (acode in .rc$codes$yea) {
        idx <- which(.rc$votes == acode)
        newrc$votes[idx] <- 1
    }
    for (acode in .rc$codes$nay) {
        idx <- which(.rc$votes == acode)
        newrc$votes[idx] <- -1
    }
    for (acode in .rc$codes$notInLegis) {
        idx <- which(.rc$votes == acode)
        newrc$votes[idx] <- 9
    }
    ## GENERALIZE THIS
    for (acode in .rc$codes$missing) {
        idx <- which(is.na(.rc$votes))
        newrc$votes[idx] <- 0
    }
    ## document new codes
    newrc$codes$yea <- 1
    newrc$codes$nay <- -1
    newrc$codes$missing <- 0
    newrc$codes$notInLegis <- 9
    ##
    return(newrc)
}
