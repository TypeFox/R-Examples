discrep.default <-
function(datavec, targetvec, weightvec) {
    require(Hmisc)
    dat <- wtd.table(datavec, weightvec)$sum.of.weights/sum(wtd.table(datavec, 
        weightvec)$sum.of.weights)
    out <- targetvec - dat
    out
}

