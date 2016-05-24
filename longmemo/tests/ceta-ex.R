library(longmemo)

(H0 <- c(seq(0.5, .999, by = 1/32), 0.98, 0.99))

mform <- function(x, digits=6, wid=9)
    formatC(x, flag="-", width=wid, digits=digits)

.proctime00 <- proc.time()

## CetaARIMA(*) does NOT depend on H !!
for(H in H0) {
    cat("H=", mform(H,wid=7),
        "; CetaFGN:", mform(CetaFGN(eta = c(H = H), m = 256)),
        "; C.ARIMA(*,0,0):",mform(CetaARIMA(eta= c(H = H), m = 256,p=0,q=0)),
        "\n")
}

for(em in 6:18) { ## oops -- becomes slow [specFGN() !] from about em = 13:
    m <- 2^em
    cAA <- CetaARIMA(eta= c(H = 0.7, phi=0.6, psi= -0.3), m = m, p=1,q=1)
    cat("m= 2^",formatC(em, wid=2),": ",
        " C_{FGN}= ", mform(CetaFGN(eta = c(H = 0.7), m = m)),
        "; C_{AR; 1,1}= ", mform(cAA[1,1]),        "\n", sep="")
    print(cAA)
    cat("\n")
}

## Last Line:
cat('Time elapsed: ', proc.time() - .proctime00,'\n')
