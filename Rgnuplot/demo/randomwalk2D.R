library(rpanel)
# generate the random walk data, 1000 observations of X, Y data for 3 experiments X1, Y1, X2, Y2, X3, Y3 use a seed for reproducibility
set.seed(0)
lenXY <- 1000
matRandXY <- matrix(0, ncol = 6, nrow = lenXY)
matRandXYcs <- matrix(0, ncol = 6, nrow = lenXY)
matRandXY <- apply(matRandXY, 2, function(x) sample(c(-1, 1), size = lenXY, replace = TRUE))
for (n in (1:3 * 2 - 1)) {
    plot(cumsum(matRandXY[, n]), cumsum(matRandXY[, n + 1]), col = n, type = "s", ylim = c(-50, 50), xlim = c(-50, 50))
    matRandXYcs[, n] <- cumsum(matRandXY[, n])
    matRandXYcs[, n + 1] <- cumsum(matRandXY[, n + 1])
    par(new = TRUE)
}
write.table((matRandXY), file = "randwalk3x1000XY.dat", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table((matRandXYcs), file = "randwalk3x1000XYcumsum.dat", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Initialize the gnuplot handle
h1 <- Gpinit()
# set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)

# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
# plot the 2D random walk in 2D
Gpcmd(h1, "#set terminal png;set output \"randomwalk2d1.png\"\nreset\nset xlabel \"X\"\nset ylabel \"Y\"\nset tit \"Random walk 2D\"\nplot \"randwalk3x1000XYcumsum.dat\" using 1:2 w steps notit,\"randwalk3x1000XYcumsum.dat\" using 3:4 w steps notit,\"randwalk3x1000XYcumsum.dat\" using 5:6 w steps notit")
# pause R and gnuplot
Gppause()

# plot the 2D random walk in 3D
Gpcmd(h1, "#set terminal png;set output \"randomwalk2d2.png\"\nreset\nset xlabel \"Time\"\nset ylabel \"X\"\nset zlabel \"Y\"\nset tit \"Random walk 2D\"\nsplot \"randwalk3x1000XYcumsum.dat\" using :1:2 w l notit, \"randwalk3x1000XYcumsum.dat\" using :3:4 w l notit, \"randwalk3x1000XYcumsum.dat\" using :5:6 w l notit")
# pause R and gnuplot
Gppause()

# close gnuplot handle
h1 <- Gpclose(h1) 
