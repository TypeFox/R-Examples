# example1.R -- version 2011-03-18
# ... fit Nelson-Siegel to yields
require(NMOF)

penalty <- function(mP, data) {
    minV <- data$min
    maxV <- data$max
    ww <- data$ww
    # if larger than maxV, element in A is positiv
    A <- mP - as.vector(maxV); A <- A + abs(A)
    # if smaller than minV, element in B is positiv
    B <- as.vector(minV) - mP; B <- B + abs(B)
    # beta 1 + beta2 > 0
    C <- ww * ((mP[1, ] + mP[2, ]) - abs(mP[1, ] + mP[2, ]))
    A <- ww * colSums(A + B) - C
    A
}
OF <- function(param, data) {		
    y   <- data$model(param, data$tm)
    aux <- y - data$yM
    aux <- max(abs(aux))
    if(is.na(aux)) 1e5 else aux 
}

# set up yield curve (here: artificial data), and plot it
tm <- c(c(1, 3, 6, 9)/12, 1:10)
betaTRUE <- c(6, 3, 8, 1)
yM   <- NS(betaTRUE, tm)
plot(tm, yM, 
    xlab = "maturities in years",
    ylab = "yields in %")

# collect everything in data
data <- list(
    yM = yM, 
    tm = tm, 
    model = NS, 
    ww = 0.1,
    min = c( 0,-15,-30, 0), 
    max = c(15, 30, 30,10))

# initialize algo settings
algo <- list(
    nP = 100L, 
    nG = 500L, 
    F = 0.50, 
    CR = 0.99,
    min = c( 0,-15,-30, 0), 
    max = c(15, 30, 30,10),
    pen = penalty, 
    repair = NULL,
    loopOF = TRUE, 
    loopPen = FALSE, 
    loopRepair = TRUE)

system.time(sol <- DEopt(OF = OF, algo = algo, data = data))
# max. error and objective function value: should be the same
max(abs(data$model(sol$xbest, tm) - data$model(betaTRUE, tm)))
sol$OFvalue
lines(tm, data$model(sol$xbest, tm), col = "blue")

# test: nlminb with random starting value s0
s0 <- algo$min + (algo$max - algo$min) * runif(length(algo$min))
system.time(sol2 <- nlminb(s0, OF, data = data,
        lower = data$min, 
        upper = data$max, 
        control = list(eval.max = 50000L, iter.max = 50000L)))
# max. error and objective function value: should be the same
max(abs(data$model(sol2$par, tm) - data$model(betaTRUE, tm)))
sol2$objective	
lines(tm, data$model(sol2$par, tm), col = "green", lty = 2)
legend(x = "bottom", 
    legend = c("true yields", "DE", "nlminb"),
    col = c("black", "blue", "green"),
    pch = c(1, NA, NA), lty = c(0, 1, 2))