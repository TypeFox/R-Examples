## =============================================================================
## Schematic representation of the Jacobian matrix
## Figure 9.5 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

require(diagram)
windows(8, 5)
par (mfrow = c(1, 2))
nn <- 15
A                        <- matrix(nrow = nn, ncol = nn, data = 0)
cc <- nn:1
A[cbind((nn):1, 1:nn)] <- 1
A[cbind(nn:2, 2:nn)] <- 1
A[cbind((nn-1):1, 1:(nn-1))] <- 1
#A[cbind(2:nn, 1:(nn-1))] <- 1

#image(t(A), col = c(0,1))
wA <- which (A != 0, arr.ind = TRUE)
plot(wA, pch = 15, cex = 1.5, xlab = "", axes = FALSE, ylab = "", 
     frame.plot = TRUE, main = "1-D, 15 cells")
writelabel("A")
A[5,6] <- A[6, 5] <- A[11,10] <- A[10, 11] <- 0
A[cbind((10:1),(1:10))] <- 1
A[cbind(15:6, 6:15)] <- 1
wA <- which (A != 0, arr.ind = TRUE)
plot(wA, pch = 15, cex = 1.5, xlab = "", axes = FALSE, ylab = "", 
    frame.plot = TRUE, main = "2-D, 5*3 grid")

abline(h = c(10.5, 5.5))
abline(v = c(10.5, 5.5))

writelabel("B")
