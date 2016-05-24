cohort <-
function(year, age, reset=F){
    n <- length(year); k <- length(age)
    stx <- year[1]-age[k]
    ft <- gl(n, k)
    ret <- array(year[ft]-age, dim=c(k,n), dimnames=list(age, year))
    if (reset) ret-stx+1 else ret
}
