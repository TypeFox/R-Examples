waveletshift.dwt <- function (L, j, N = NULL)
{
    if(L%%2 == 1){
        stop("Filter must be of even length.")
    }
    Lj <- ((2^j)-1)*(L-1) + 1
    if (L == 10 || L == 18) {
        vjH <- (-Lj/2) + 1
    }
    else if (L == 14) {
        vjH <- (-Lj/2) - 1
    }
    else {
        vjH <- -Lj/2
    }
    if(is.null(N) ==  FALSE) {
        shift <- abs(vjH)%%N
    }
    else {
        shift <- abs(vjH)
    }
    shift
}
