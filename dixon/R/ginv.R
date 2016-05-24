`ginv` <-
function (m) 
{
    temp <- eigen(m, symmetric = TRUE)
    va <- temp$values
    ve <- temp$vectors
    va <- ifelse((abs(va) < 1e-09), 0, 1/va)
    va2 <- 0 * m
    diag(va2) <- va
    ve %*% va2 %*% t(ve)
}

