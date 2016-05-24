
# Initialize the gnuplot handle
h1 <- Gpinit()
# use a seed for reproducibility
set.seed(0)
nxpoints <- 10  #number of x points to plot
nypoints <- npoints + sample.int(nxpoints, size = 1)  # add a random number from [1, nxpoints]
xpoints <- sort(sample.int(nxpoints, size = nypoints, replace = TRUE))
ypoints <- trunc(rnorm(nypoints, mean = 1, sd = 3))  #generate random points and truncate them
strgnu <- "\n" %s% paste(xpoints, ypoints, sep = "\t", collapse = "\n") %s% "\ne"
Gpcmd(h1, "#set terminal png;set output \"interpolation2.png\"\nreset\nset xrange [0:" %s% (npoints + 1) %s% "]\nset yrange [" %s% (min(ypoints) - 1) %s% ":" %s% (max(ypoints) + 1) %s% 
    "]\nset key top center # positions the legend\nset multiplot\nset size ratio -1\nplot \"-\" notit with points" %s% strgnu)
Gpcmd(h1, "plot \"-\" smooth csplines notit  w l lc 2" %s% strgnu)
Gpcmd(h1, "plot \"-\" smooth bezier notit  w l lc 3" %s% strgnu)
Gpcmd(h1, "plot \"-\" smooth unique notit  w l lc 4" %s% strgnu)
Gpcmd(h1, "plot \"-\" smooth frequency notit  w l lc 5" %s% strgnu)
Gpcmd(h1, "plot \"-\" smooth sbezier notit  w l lc 6" %s% strgnu)
Gpcmd(h1, "plot \"-\" smooth acsplines notit  w l lc 7" %s% strgnu)
Gpcmd(h1, "plot NaN notit,  NaN tit \"splines\", NaN tit \"bezier\", NaN tit \"unique\", NaN tit \"frequency\", NaN tit \"sbezier\", NaN tit \"acsplines\"")
# pause R and gnuplot
Gppause()
# close gnuplot handle
h1 <- Gpclose(h1) 
