
# polynomial curve fitting Using Gnuplot to sketch graphs Aaron Titus, 2012 http://linus.highpoint.edu/~atitus/gnuplot/
polnorder <- 7  # order of the polynomial
npoints <- 20  # number of points to plot
xpoints <- (0:npoints) * 0.1  # x values
wpoints <- c(1, 10^-(0:polnorder))  # 'a' to 'h' values
# calculate y values
ypoints <- wpoints[1] + wpoints[2] * xpoints + wpoints[3] * xpoints^2 + wpoints[4] * xpoints^3 + wpoints[5] * xpoints^4 + wpoints[6] * xpoints^5 + wpoints[7] * xpoints^6 + wpoints[8] * 
    xpoints^7
# write the data to a text file
write(t(cbind(xpoints, ypoints)), "data.txt", ncolumns = 2)
# initial values of the parameters
initvalues <- c(1, 1, 1, 1, 1, 1, 1, 1)
# write the initial values of the parameters to a text file
write(t(cbind(letters[1:(polnorder + 1)], "=", initvalues)), "guess.txt", ncolumns = 3, sep = "")
# create gnuplot handle
h1 <- Gpinit()
# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
Gpcmd(h1, "#set terminal postscript eps color;set output \"polynomialfit.eps\"\ny(x)=a+b*x+c*x**2+d*x**3+e*x**4+f*x**5+g*x**6+h*x**7\nfit y(x) \"data.txt\" via \"guess.txt\"\nset xlabel \"x\"\nset ylabel \"y\"\nunset key\nplot y(x), \"data.txt\"")

# pause R and gnuplot
Gppause()
# close gnuplot handle
h1 <- Gpclose(h1) 
