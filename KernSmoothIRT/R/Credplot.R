Credplot <-function (OBJ, axis, subjects, quants, main, xlab, ylab, ...) 
{
   
    if (missing(ylab)) {
        ylab = "Relative Credibility"
    }
    if (missing(subjects)) {
        subjects <- 1:OBJ$nsubj
    }
    if (missing(main)) {
        main = -1
    }
    par(ask = TRUE)
    for (i in subjects) {
       

        relcred <- OBJ$RCC[[i]]
        if (main == -1) {
            main0 <- paste("Subject: ", i, "\n")
        }
        plot(axis, relcred, type = "l", main = main0, ylab = ylab, 
            xlab = xlab, ...)
        axis(3, at = quants, lab = labels(quants), tck = 0)
        abline(v = quants, col = "blue", lty = 2)
        if (axis[1] == OBJ$evalpoints[1]) {
            abline(v = OBJ$subjtheta[i], col = "red")
        }
        else {
            abline(v = OBJ$subjscore[i], col = "red")
        }
    }
}

