
# based on the example from gplots
data(state)
tmp <- split(state.area, state.region)
means <- sapply(tmp, mean)
stdev <- sqrt(sapply(tmp, var))
n <- sapply(tmp, length)
ciw <- qt(0.975, n) * stdev/sqrt(n)

# Initialize the gnuplot handle
h1 <- Gpinit()
# set gnuplot's additional search directories, to the extdata directory from Rgnuplot (default)
Gpsetloadpath(h1)

# change gnuplot's working directory to be the same as R's working directory (default)
Gpsetwd(h1)

tmpfile <- tempfile()
xtic <- ""
gplabelmean <- ""
gplabeln <- ""
for (h in 1:4) {
    # gpdata <- gpdata %s% h %s% ' ' %s% means[h] %s% ' ' %s% (means[h]-ciw[h]) %s% ' ' %s% (means[h]+ciw[h]) %s% '\n'
    cat(h, " ", means[h], " ", (means[h] - ciw[h]), " ", (means[h] + ciw[h]), "\n", file = tmpfile, append = TRUE)
    xtic <- xtic %s% "\"" %s% names(tmp)[h] %s% "\" " %s% h
    if (h < 4) 
        xtic <- xtic %s% ","
    gplabelmean <- gplabelmean %s% "set label " %s% h %s% " \"" %s% round(means[h], -3) %s% "\" at " %s% h %s% "," %s% means[h] %s% " nopoint  tc rgb \"red\"\n"
    gplabeln <- gplabeln %s% "set label " %s% h %s% " \"n=" %s% n[h] %s% "\" at " %s% h %s% "," %s% 4000 %s% " nopoint  tc rgb \"blue\"\n"
}

# error box, constant width
Gpcmd(h1, "reset\n#set terminal png;set output \"errorbox.png\"\nset xrange [0:5]\nset yrange [ 0:250000]\nset xtics (" %s% xtic %s% ")\nplot \"" %s% tmpfile %s% "\" u 1:2:3:4:(0.2) with boxerror lc rgb \"black\" tit \"Area per state\"")

# pause R and gnuplot
Gppause()

# error bars, mean values and a line
Gpcmd(h1, "reset\n#set terminal png;set output \"errorbar1.png\"\nset xrange [0:5]\nset yrange [ 0:250000]\nset xtics (" %s% xtic %s% ")\nset xtics out\nset ytics out\n" %s% gplabelmean %s% 
    "\nplot \"" %s% tmpfile %s% "\" with errorbars tit \"Area per state\", \"" %s% tmpfile %s% "\" w l notit")

# pause R and gnuplot
Gppause()

# error bars, mean values, a line and more labels
Gpcmd(h1, "reset\n#set terminal png;set output \"errorbar2.png\"\nset xrange [0:5]\nset yrange [ 0:250000]\nset xtics (" %s% xtic %s% ")\nset xtics out\nset ytics out\nset xlabel \"Region\"\nset ylabel \"Area\"\n" %s% 
    gplabeln %s% "\nplot \"" %s% tmpfile %s% "\" with errorbars tit \"Area per state\", \"" %s% tmpfile %s% "\" w l notit")

# pause R and gnuplot
Gppause()

# close gnuplot handles
h1 <- Gpclose(h1) 
