## =============================================================================
## The family of initial value problem solvers. 
## Figure 3.5 from Soetaert, Cash, Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diagram)
par(mar=c(1, 1, 1, 1), mfrow = c(1, 1))
emptyplot()
M <- matrix(nrow = 8, ncol = 8, byrow = TRUE, data = 0)
M[2,1] <- M[3,1] <- M[4,2] <- M[6,2] <- M[5,3] <- M[7,3] <- M[8,3] <-" "
pos <- coordinates(c(1, 2, 5), my = 0.05)

names <- c("Euler", "Runge-Kutta", "Linear Multistep", 
  "Explicit RK",  "Adams", "Implicit RK","BDF","MEBDF")
plotmat(M, pos, box.type = "rect", curve = 0, 
        box.size = 0.09, box.prop = 0.5, name = names,
        cex.txt = 0, arr.length = 0)

treearrow (from = cbind(c(0.005, 0.395), c(0.15, 0.15)), to = c(0.20, 0.1), 
           arr.length = 0)
text(0.20, 0.05, "non-stiff problems", font = 3, cex=1.4)
           
treearrow (from = cbind(c(0.405, 0.995), c(0.15, 0.15)), to = c(0.7, 0.1), 
           arr.length = 0)
text(0.7, 0.05, "stiff problems", font = 3, cex = 1.4)

