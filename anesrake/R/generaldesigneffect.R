generaldesigneffect <-
function(weightvec) {
    weightvec <- weightvec[!is.na(weightvec)] * sum(!is.na(weightvec))/sum(weightvec[!is.na(weightvec)])
    deff <- sum(weightvec^2)/sum(weightvec)
    deff
}

