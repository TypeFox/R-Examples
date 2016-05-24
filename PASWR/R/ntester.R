ntester <-
function(actual.data)
{
    Ared <- "#C00000"
    Ablue <- "#0080FF"
    par(mfrow = c(3, 3))
    par(oma = c(1, 0, 2, 0))
    par(mar = c(0, 0, 2, 0))
    par(pty = "s")
    for(i in 1:4) {
        SimData <- rnorm(length(actual.data))
        s <- shapiro.test(SimData)
        qqnorm(SimData, xlab = "", ylab = "", axes = FALSE, col = Ablue,
        main=paste("SimNorm p-val = ", round(s$p.value, 3)),col.main=Ablue)
        box()
        qqline(SimData, col = Ared)
    }
    sx <- shapiro.test(actual.data)
    qqnorm(actual.data, xlab = "", ylab = "", axes =FALSE, col = Ared,
    main = paste("YourData p-val = ", round(sx$p.value, 3)),col.main=Ared)
    box()
    qqline(actual.data, col = Ablue)
    for(i in 6:9) {
        SimData <- rnorm(length(actual.data))
        s <- shapiro.test(SimData)
        qqnorm(SimData, xlab = "", ylab = "", axes = FALSE, col = Ablue,
        main= paste("SimNorm p-val = ", round(s$p.value, 3)),col.main=Ablue)
        box()
        qqline(SimData, col = Ared)
    }
    mtext("Simulated Normal Data on Perimeter - Actual Data in Center",
    side = 3, outer = TRUE, cex = 1.5, col = Ared)
    par(oma = c(0, 0, 0, 0))
    par(mfrow = c(1, 1))
    par(mar = c(5.1, 4.1, 4.1, 2.1))
    par(pty = "m")
}

