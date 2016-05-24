# Nash Sutcliffe efficiency calculation
# Wouter Buytaert
# 2007/07/05

NSeff <- function(Qobs,Qsim) {

    # remove NA's:
    
    Qsim <- Qsim[!is.na(Qobs)]
    Qobs <- Qobs[!is.na(Qobs)]
    
    Qobs <- Qobs[!is.na(Qsim)]
    Qsim <- Qsim[!is.na(Qsim)]

    if(length(Qobs) == 0 || length(Qsim) == 0) return(NA)
    
    # calculate efficiency

    NS <- 1 - ( sum((Qobs - Qsim)^2) / sum((Qobs - mean(Qobs))^2) )
    return(NS)
}
