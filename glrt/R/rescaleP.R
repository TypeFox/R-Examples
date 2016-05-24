rescaleP <-
function (pvec, tiny) 
{
    pvec <- ifelse(pvec < tiny, 0, pvec)
    pvec <- pvec/sum(pvec)
    return(pvec)
}