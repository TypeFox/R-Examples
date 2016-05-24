##
## Unit testing of Linear Program with non-linear and cone constraints
test.FPLN <- function(){
    if(requireNamespace("numDeriv", quietly = TRUE)){
        ## Creating objective
        q <- c(rep(1, 2), rep(0, 20))
        xnames <- c("W", "H",
                    paste("x", 1:5, sep = ""),
                    paste("y", 1:5, sep = ""),
                    paste("w", 1:5, sep = ""),
                    paste("h", 1:5, sep = "")
                   )
        ## Fixed constants
        gamma <- 5.0
        rho <- 1.0
        Amin <- 100
        ## Inequality constraints
        G <- matrix(0.0, nrow = 26, ncol = 22)
        h <- matrix(0.0, nrow = 26, ncol = 1)
        G[1, 3] <- -1.0                                       ## -x1 <= 0
        G[2, 4] <- -1.0                                       ## -x2 <= 0
        G[3, 6] <- -1.0                                       ## -x4 <= 0
        G[4, c(3, 5, 13)] <- c(1.0, -1.0, 1.0)                ## x1 - x3 + w1 <= -rho
        h[4, 1] <- -rho
        G[5, c(4, 5, 14)] <- c(1.0, -1.0, 1.0)                ## x2 - x3 + w2 <= -rho
        h[5, 1] <- -rho
        G[6, c(5, 7, 15)] <- c(1.0, -1.0, 1.0)                ## x3 - x5 + w3 <= -rho
        h[6, 1] <- -rho
        G[7, c(6, 7, 16)] <- c(1.0, -1.0, 1.0)                ## x4 - x5 + w4 <= -rho
        h[7, 1] <- -rho
        G[8, c(1, 7, 17)] <- c(-1.0, 1.0, 1.0)                ## -W + x5 + w5 <= 0
        G[9, 9] <- -1.0                                       ## -y2 <= 0
        G[10, 10] <- -1.0                                     ## -y3 <= 0
        G[11, 12] <- -1.0                                     ## -y5 <= 0
        G[12, c(8, 9, 19)] <- c(-1.0, 1.0, 1.0)               ## -y1 + y2 + h2 <= -rho
        h[12, 1] <- -rho
        G[13, c(8, 11, 18)] <- c(1.0, -1.0, 1.0)              ##  y1 - y4 + h1 <= -rho
        h[13, 1] <- -rho
        G[14, c(10, 11, 20)] <- c(1.0, -1.0, 1.0)             ##  y3 - y4 + h3 <= -rho
        h[14, 1] <- -rho
        G[15, c(2, 11, 21)] <- c(-1.0, 1.0, 1.0)              ## -H + y4 + h4 <= 0
        G[16, c(2, 12, 22)] <- c(-1.0, 1.0, 1.0)              ## -H + y5 + h5 <= 0
        G[17, c(13, 18)] <- c(-1.0, 1.0 / gamma)              ## -w1 + h1/gamma <= 0
        G[18, c(13, 18)] <- c(1.0, -gamma)                    ##  w1 - gamma * h1 <= 0
        G[19, c(14, 19)] <- c(-1.0, 1.0 / gamma)              ## -w2 + h2/gamma <= 0
        G[20, c(14, 19)] <- c(1.0, -gamma)                    ##  w2 - gamma * h2 <= 0
        G[21, c(15, 19)] <- c(-1.0, 1.0 / gamma)              ## -w3 + h3/gamma <= 0
        G[22, c(15, 20)] <- c(1.0, -gamma)                    ##  w3 - gamma * h3 <= 0
        G[23, c(16, 20)] <- c(-1.0, 1.0 / gamma)              ## -w4  + h4/gamma <= 0
        G[24, c(16, 21)] <- c(1.0, -gamma)                    ##  w4 - gamma * h4 <= 0
        G[25, c(16, 21)] <- c(-1.0, 1.0 / gamma)              ## -w5 + h5/gamma <= 0
        G[26, c(17, 22)] <- c(1.0, -gamma)                    ## w5 - gamma * h5 <= 0
        nno1 <- nnoc(G = G, h = h)
        ## Nonlinear constraints
        f1 <- function(x) -x[13] + Amin / x[18]
        f2 <- function(x) -x[14] + Amin / x[19]
        f3 <- function(x) -x[15] + Amin / x[20]
        f4 <- function(x) -x[16] + Amin / x[21]
        f5 <- function(x) -x[17] + Amin / x[22]
        ## Gradient functions
        g1 <- function(x, func = f1) numDeriv::grad(func = func, x = x)
        g2 <- function(x, func = f2) numDeriv::grad(func = func, x = x)
        g3 <- function(x, func = f3) numDeriv::grad(func = func, x = x)
        g4 <- function(x, func = f4) numDeriv::grad(func = func, x = x)
        g5 <- function(x, func = f5) numDeriv::grad(func = func, x = x)
        ## Hessian functions
        h1 <- function(x, func = f1) numDeriv::hessian(func = func, x = x)
        h2 <- function(x, func = f2) numDeriv::hessian(func = func, x = x)
        h3 <- function(x, func = f3) numDeriv::hessian(func = func, x = x)
        h4 <- function(x, func = f4) numDeriv::hessian(func = func, x = x)
        h5 <- function(x, func = f5) numDeriv::hessian(func = func, x = x)
        ## Initial value
        x0 <- rep(1, 22)
        ## Invoking 'cccp'
        ans <- cccp(q = q, cList = list(nno1), x0 = x0,
                    nlfList = list(f1, f2, f3, f4, f5),
                    nlgList = list(g1, g2, g3, g4, g5),
                    nlhList = list(h1, h2, h3, h4, h5))
        checkTrue(ans$status == "optimal")
    }
    return()
}
