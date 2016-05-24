rain.adapt <-
function (x, a, ser) 
{
    a1 <- numeric()
    a2 <- numeric()
    a1bis <- numeric()
    cfreq <- contaleva(x, a, ser)
    a1 <- a - cfreq
    a1bis <- Levaneg(a1)
    cvol <- aggiustam(x, a1bis, ser)
    volx <- sum(x)
    vola <- sum(a1bis)
    if (vola > volx) 
        cvol <- 0 - cvol
    a2 <- aggiusta(a1bis, cvol)
    if (any(a2 < 0)) {
        cat("You can't fix the same volume", "\n")
        a2 = a1bis
    }
    a2
}
