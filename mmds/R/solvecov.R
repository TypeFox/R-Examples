# taken from mrds

solvecov<-function (m, cmax = 1e+10)
# from package fpc
{
    options(show.error.messages = FALSE)
    covinv <- try(solve(m))
    if (class(covinv) != "try-error")
        coll = FALSE
    else {
        p <- nrow(m)
        cove <- eigen(m, symmetric = TRUE)
        coll <- TRUE
        if (min(cove$values) < 1/cmax) {
            covewi <- diag(p)
            for (i in 1:p) if (cove$values[i] < 1/cmax)
                covewi[i, i] <- cmax
            else covewi[i, i] <- 1/cove$values[i]
        }
        else covewi <- diag(1/cove$values, nrow = length(cove$values))
        covinv <- cove$vectors %*% covewi %*% t(cove$vectors)
    }
    options(show.error.messages = TRUE)
    out <- list(inv = covinv, coll = coll)
}
