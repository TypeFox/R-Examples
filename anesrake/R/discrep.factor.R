discrep.factor <-
function(datavec, targetvec, weightvec) {
    require(Hmisc)
    dat <- sapply(names(targetvec), function(x) {
        sum(weightvec[datavec == x & !is.na(datavec)])/sum(weightvec[!is.na(datavec)])
    })
    out <- targetvec - dat
    out
}

