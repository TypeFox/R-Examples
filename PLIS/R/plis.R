plis <-
function (lis, fdr = 0.001, adjust = TRUE) 
{
    n = length(lis)
    s.lis = sort(lis)
    for (i in 1:n) {
        if (mean(s.lis[1:i]) > fdr) 
            break
    }
    nNonNull = i - 1
    States = rep(0, n)
    if (nNonNull > 0) 
        States[lis <= s.lis[nNonNull]] = 1
    if (adjust) {
        aLIS = sapply(lis, function(cut) mean(lis[which(lis <= 
            cut)]))
        return(list(States = States, aLIS = aLIS))
    }
    else {
        return(list(States = States))
    }
}
