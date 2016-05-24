###
### Projection utilities
###

orth <- function(M, N) {                # O_N M
    ## FIXME: This must be optimised:
    M - N %*% solve(crossprod(N)) %*% crossprod(N, M)
}
## This function is currently not used in the code:
project <- function(M, N) {  #P_N, M)
    ## FIXME: Must be optimised:
    N %*% solve(crossprod(N)) %*% crossprod(N, M)
}
Corth <- function(M, N) {               #C_N M
    ## FIXME: This must be optimised:
    solve(crossprod(N)) %*% crossprod(N, M)
}
