VGGFR<-function (L1, L2 = NULL, add = FALSE, lwd = 2, lty = 5, col = "blue", ylim=c(0,1.0), Main="", np = 201) {
	if (is.null(L2)){L2<-L1[2];L1<-L1[1]}  
	if (!is.numeric(L1) | !is.numeric(L2)) {stop("Non-numeric argument to mathematical function")}
	if (!(L1>0) | !(L2>0)) {stop("Shape parameters must be positive")}
	x<-NULL
    Vian <- function(x, L1, L2) {L1 * (1 - abs(x)^L1)^L2/(2 * beta(1/L1, 1 + L2))}
    curve(Vian(x, L1, L2), -1, 1, np, col = col, lwd = lwd, add = add, xlab="Domain of the random variable", ylab="Probability density", ylim=ylim, main=Main)
    Ws1 <- beta(3/L1, L2 + 1)
    Ws2 <- beta(1/L1, L2 + 1)
    Ws3 <- beta(5/L1, L2 + 1)
    s2 <- Ws1/Ws2
    mu4 <- (Ws3/Ws2)
    g2 <- mu4/(s2 * s2) - 3
    h2 <- L1/(2 * Ws2)
    s2 <- round(sqrt(s2), 6)
    g2 <- round(g2, 6)
    h2 <- round(h2, 6)
    cat("Parameters of VGGFR: ", c(L1, L2), "\n")
    cat("Std. Deviation: ", s2, "     Excess Kurtosis: ", g2, 
        "   Ordinate at the mode:", h2, "\n")
    rea <- list(st.dev = s2, kurt = g2, oam = h2)
    return(rea)
}