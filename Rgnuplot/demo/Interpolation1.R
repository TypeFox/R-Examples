
# Initialize the gnuplot handle
h1 <- Gpinit()
# use a seed for reproducibility
set.seed(0)
# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)
npoints <- 10  #number of points to plot
randomPoints <- rnorm(npoints, mean = 1, sd = 3)  #generate random points
strgnu <- "\n" %s% paste(1:npoints, randomPoints, sep = "\t", collapse = "\n") %s% "\ne"
Gpcmd(h1, "#set terminal png;set output \"interpolation1.png\"\nreset\nset xrange [0:" %s% (npoints + 1) %s% "]\nset yrange [" %s% (min(randomPoints) - 1) %s% ":" %s% (max(randomPoints) + 
    1) %s% "]\nset key top center\nset multiplot\nset size ratio -1\nplot \"-\" notit with points" %s% strgnu)
Gpcmd(h1, "plot \"-\" smooth csplines notit  w l lc 2" %s% strgnu)
Gpcmd(h1, "plot \"-\" smooth bezier notit  w l lc 3" %s% strgnu)
Gpcmd(h1, "plot NaN notit,  NaN tit \"splines\", NaN tit \"bezier\"")
# pause R and gnuplot
Gppause()
# close gnuplot handle
h1 <- Gpclose(h1) 
