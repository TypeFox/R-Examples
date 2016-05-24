whichlevel <-
function (J, filter.number = 10, family = "DaubLeAsymm") 
{
    BigJ <- J + 1
    KeepGoing <- T
    while (KeepGoing) {
        n <- 2^BigJ
        tmpwd <- wd(rep(0, n), filter.number = filter.number, 
            family = family)
        ixvec <- cumsum(2^((BigJ - 1):0))
        somefull <- F
        for (i in 1:J) {
            tmpwd$D <- rep(0, n - 1)
            tmpwd$D[ixvec[i]] <- 1
            vec <- wr(tmpwd)
            somefull <- somefull | all(vec != 0)
            if (somefull) 
                break
        }
        if (!somefull) 
            KeepGoing <- F
        else BigJ <- BigJ + 1
    }
    BigJ
}
