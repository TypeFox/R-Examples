rakeonvar.factor <-
function(weighton, weightto, weightvec) {
    weighton <- as.factor(weighton)
    if(sum(is.na(weighton))>0){
      levels(weighton) <- c(levels(weighton), "(Missing)")
      weighton[is.na(weighton)] <- "(Missing)"
    }
    lwo <- length(table(weighton))
    lwt <- length(weightto)
    pctchk <- 0
    if (length(names(weightto)) != length(unique(names(weightto)))) {
        stop("only one weighting target can be used per named variable")
    }
    if (sum(weightto) < 1.5) {
        weightto <- weightto * sum(weightvec)
    }
    if (lwo > (lwt + 1)) {
        stop("number of variable levels does not match number of weighting levels")
    }
    if (sum(levels(weighton) %in% c(names(weightto), "(Missing)")) != 
        length(levels(weighton))) {
        stop("variable levels and weighting target levels must match")
    }
    if (sum(names(weightto) %in% levels(weighton)) != length(names(weightto))) {
        stop("variable levels and weighting target levels must match")
    }
    if (lwo == (lwt + 1)) {
        mis <- sum(weightvec[weighton == "(Missing)"])
        nm <- c(names(weightto), "(Missing)")
        weightto <- c(weightto * ((sum(weightvec) - mis)/sum(weightto)), 
            mis)
        names(weightto) <- nm
    }
    for (i in levels(weighton)) {
        weightvec[weighton == i] <- weightvec[weighton == i] * 
            (weightto[names(weightto) == i]/sum(weightvec[weighton == 
                i]))
    }
    weightvec
}

