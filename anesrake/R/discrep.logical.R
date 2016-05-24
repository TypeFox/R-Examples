discrep.logical <-
function(datavec, targetvec, weightvec) {
    require(Hmisc)
    dat <- wtd.table(datavec, weightvec)$sum.of.weights/sum(wtd.table(datavec, 
        weightvec)$sum.of.weights)
    out <- c(targetvec[2] - dat[1], targetvec[1] - dat[2])
    out
}

