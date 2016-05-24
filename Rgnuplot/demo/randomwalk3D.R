library(rpanel)
# generate the random walk data, 1000 observations of X, Y data for 3 experiments X1, Y1, X2, Y2, X3, Y3 use a seed for reproducibility
set.seed(0)
lenXYZ <- 1000
matRandXYZ <- matrix(0, ncol = 30, nrow = lenXYZ)
matRandXYZcs <- matrix(0, ncol = 30, nrow = lenXYZ)
matRandXYZ <- apply(matRandXYZ, 2, function(x) sample(c(-1, 1), size = lenXYZ, replace = TRUE))
for (n in 1:30) matRandXYZcs[, n] <- cumsum(matRandXYZ[, n])
# 
write.table((matRandXYZ), file = "randwalk10x1000XYZ.dat", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table((matRandXYZcs), file = "randwalk10x1000XYZcumsum.dat", sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

# Initialize the gnuplot handle
h1 <- Gpinit()
# set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)

# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1, "#set terminal png;set output \"randomwalk3d1.png\"\nreset\nset xlabel \"X\"\nset ylabel \"Y\"\nset zlabel \"Z\"\nset tit \"Random walk 3D\";splot \"randwalk10x1000XYZcumsum.dat\" using 1:2:3 w l notit, \"randwalk10x1000XYZcumsum.dat\" using 4:5:6 w l notit, \"randwalk10x1000XYZcumsum.dat\" using 7:8:9 w l notit")
# pause R and gnuplot
Gppause()

# close gnuplot handle
h1 <- Gpclose(h1) 
