jonckheere <-
function(y, groups)
{
        groups <- as.factor(groups)
        lev <- levels(groups)
        nlev <- length(lev)
 
        jonck <- 0
        ties.p <- F
        for(i in 1:(nlev - 1))
                for(j in (i + 1):nlev) {
                        xx <- y[groups == lev[i]]
                        yy <- y[groups == lev[j]]
                        xy <- c(xx, yy)
                        ties.p <- ties.p | (length(xy) !=length(unique(xy)))
                        sss <- length(xx) * length(yy) + (length(xx) *
                                (length(xx) + 1))/2 - 
sum(rank(xy)[1:length(xx)])
                        jonck <- jonck + sss
                }
#for ties see Lehmann p.235ff
        names(jonck) <- NULL
        ns <- tabulate(as.numeric(groups))
        nn <- length(y)
        expj <- (nn * nn - sum(ns * ns))/4
        if(!ties.p) {
                varj <- (nn * nn * (2 * nn + 3) - sum(ns * ns * (2 * ns + 
3)))/72
        }
        else {
                nun <- as.vector(table(y))
                varj <- (nn * (nn - 1) * (2 * nn + 5) - 
                                sum(ns * (ns - 1) * (2 * ns + 5)) - 
                                sum(nun * (nun - 1) * (2 * nun + 5)))/72 + 
                                (sum(ns * (ns - 1) * (ns -      2)) * sum(nun * 
(nun - 1) * 
                                (nun - 2)))/(36 * nn * (nn - 1) * (nn - 2)) + 
                                (sum(ns * (ns - 1)) * sum(nun * (nun - 1)))/(8 
* nn * (nn - 1))
        }
        c(Jonckheere = jonck, ExpJ = expj, VarJ = varj, 
                P = 1 - pnorm((jonck- expj)/sqrt(varj)))
}
