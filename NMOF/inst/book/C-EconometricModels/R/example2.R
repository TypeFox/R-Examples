# example2.R -- version 2011-01-06
# ... fit Nelson-Siegel-Svensson to yields
# ... continues example1.R

# set up yield curve (here: artificial data), and plot it
tm <- c(c(1, 3, 6, 9)/12, 1:10)
betaTRUE <- c(5, -2, 5, -5, 1, 6)
yM <- NSS(betaTRUE, tm) 
plot(tm, yM, 
    xlab = "maturities in years",
    ylab = "yields in %")

# collect everything in data
data <- list(
    yM = yM,
    tm = tm, 
    model = NSS,
    min = c( 0,-15,-30,-30,  0,  5),
    max = c(15, 30, 30, 30,  5, 10),
    ww = 1)

algo <- list(
    nP = 100L,
    nG = 500L,
    F = 0.50,
    CR = 0.99,
    min = c( 0,-15,-30,-30, 0, 5),
    max = c(15, 30, 30, 30, 5, 10),
    pen = penalty,
    repair = NULL,
    loopOF = TRUE, 
    loopPen = FALSE,
    loopRepair = TRUE)

system.time(sol <- DEopt(OF = OF, algo = algo, data = data) )
# max. error and objective function value: should be the same
max(abs(data$model(sol$xbest, tm) - data$model(betaTRUE, tm)))
sol$OFvalue
lines(tm,data$model(sol$xbest, tm), col = "blue")

# test: nlminb with random starting value s0
s0 <- algo$min + (algo$max - algo$min) * runif(length(algo$min))
system.time(sol2 <- nlminb(s0, OF, data = data,
        lower = data$min, 
        upper = data$max, 
        control = list(eval.max = 50000L, iter.max = 50000L)))
# max. error and objective function value: should be the same
max(abs(data$model(sol2$par, tm) - data$model(betaTRUE, tm)))
sol2$objective
lines(tm,data$model(sol2$par, tm), col = "green", lty = 2)

legend(x  = "bottom", 
    legend = c("true yields", "DE", "nlminb"),
    col = c("black", "blue", "green"),
    pch = c(1, NA, NA), lty = c(0, 1, 2))