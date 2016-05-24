
# Initialize the gnuplot handle
h1 <- Gpinit()
# set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)
# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
# create the fractal data from gnuplot
Gpcmd(h1, "reset\ncomplex(x,y) = x*{1,0}+y*{0,1}\nmandel(x,y,z,n) = (abs(z)>2.0 || n>=10) ? n : mandel(x,y,z*z+complex(x,y),n+1)\nset isosamples 50,50\nset samples 30,30\nset xrange [-1.5:0.5]\nset yrange [-1:1]\nset logscale z\nset hidden3d\nset table \"mandel.dat\"\nsplot mandel(x,y,{0,0},0) notitle\nunset table")
# plot the fractal data on 3D
Gpcmd(h1, "#set terminal png;set output \"fractal1.png\"\nreset\nsplot \"mandel.dat\" w  pm3d notitle")
# plot the fractal data with a heat map
Gpcmd(h1, "#set terminal png;set output \"fractal2.png\"\nreset\nset pm3d map\nsplot \"mandel.dat\" w  pm3d notitle")
# plot the fractal data with a contour plot
Gpcmd(h1, "reset\nunset surface\nset contour\nset pm3d map\nunset key\nsplot \"mandel.dat\" w  pm3d notitle")
# create the fractal data from R calling an R function
mandel <- function(x, y, z, n) if ((abs(z) > 2) | (n >= 10)) n else mandel(x, y, z * z + complex(2, x, y), n + 1)
mandelxy1 <- function(x, y) mandel(x, y, c(0, 0), 0)
GpR2splot(mandelxy1, "mandel2.dat", c(-1.5, 0.5), c(-1, 1), c(50, 50), c(30, 30), TRUE)
Gpcmd(h1, "reset\nsplot \"mandel2.dat\" w  pm3d notitle")
# create the fractal data from R calling a C function
mandelxy2 <- function(x, y) Gpmandel(x, y, c(0, 0), 0, 10)
GpR2splot(mandelxy2, "mandel3.dat", c(-1.5, 0.5), c(-1, 1), c(50, 50), c(30, 30), TRUE)
Gpcmd(h1, "reset\nsplot \"mandel3.dat\" w  pm3d notitle")
# create the fractal data from R calling a C function, with more points and more iterations
mandelxy2 <- function(x, y) Gpmandel(x, y, c(0, 0), 0, 1000)
GpR2splot(mandelxy2, "mandel4.dat", c(-1.5, 0.5), c(-1, 1), c(500, 500), c(300, 300), TRUE)
Gpcmd(h1, "reset\nsplot \"mandel4.dat\" w  pm3d notitle")
# plot the fractal data with a heat map
Gpcmd(h1, "#set terminal png;set output \"fractal3.png\"\nreset\nset pm3d map\nsplot \"mandel4.dat\" w  pm3d notitle")
# pause R and gnuplot
Gppause()
# close gnuplot handles
h1 <- Gpclose(h1)

 
