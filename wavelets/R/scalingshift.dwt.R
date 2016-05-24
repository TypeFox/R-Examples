scalingshift.dwt <- function (L, j, N = NULL)
{
    if(L%%2 == 1){
        stop("Filter must be of even length.")
    }
    Lj <- ((2^j)-1)*(L-1) + 1
    if (L == 10 || L == 18) {
        vjG <- -((Lj-1)*L)/(2*(L-1))
    }
    else if (L == 14) {
        vjG <- -((Lj-1)*(L-4))/(2*(L-1))
    }
    else {
        vjG <- -((Lj-1)*(L-2))/(2*(L-1))
    }
    if(is.null(N) ==  FALSE) {
        shift <- abs(vjG)%%N
    }
    else {
        shift <- abs(vjG)
    }
    shift
}
