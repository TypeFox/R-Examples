siardemo <-
function(siarversion=0) {

cat("==================== Demo ==================== \n \n")
cat("This is a simple demo using the Geese data provided \n")
cat("with the package. The data are loaded into the R workspace \n")
cat("but are also available as text files in the package directory. \n")
cat("The data are called geese1demo, sourcesdemo, correctionsdemo \n")
cat("and concdepdemo. \n \n")
cat("This example deals with data with 1 group and 2 isotopes. \n")

cat("Press <Enter> to see the data \n")
readline()
invisible()

geese1demo <- matrix(c(10.22,10.37,10.44,10.52,10.19,10.45,9.91,11.27,9.34,-11.36,-11.88,-10.6,-11.25,-11.66,-10.41,-10.88,-14.73,-11.52),ncol=2,nrow=9)
colnames(geese1demo) <- c("d15NPl","d13CPl")
sourcesdemo <- matrix(c(6.48898447009942,4.4321601032975,11.1926127956759,9.8162797508688,1.45946324319852,2.2680708954525,1.112438464109,0.827103932159427,-11.1702276840033,-30.8798439532492,-11.1709000370306,-14.0570109233264,1.21495615708948,0.641318210364447,1.95933055832669,1.17246767177711),ncol=4,nrow=4)
sourcesdemo <- data.frame(c("Zostera","Grass","U.lactuca","Enteromorpha"),sourcesdemo)
colnames(sourcesdemo) <- c("Sources","Meand15N","SDd15N","Meand13C","SDd13C")
correctionsdemo <- matrix(c(rep(3.54,4),rep(0.74,4),rep(1.63,4),rep(0.63,4)),ncol=4,nrow=4)
correctionsdemo <- data.frame(c("Zostera","Grass","U.lactuca","Enteromorpha"),correctionsdemo)
colnames(correctionsdemo) <- c("Source","Mean15N","Sd15N","Mean13C","Sd13C")
concdepdemo <- matrix(c(0.0297,0.0355,0.0192,0.0139,0.0097,0.0063,0.0053,0.0057,0.3593,0.4026,0.2098,0.1844,0.0561,0.038,0.0327,0.1131),ncol=4,nrow=4)
concdepdemo <- data.frame(c("Zostera","Grass","U.lactuca","Enteromorpha"),concdepdemo)
colnames(concdepdemo) <- c("Sources","Meand15N","SDd15N","Meand13C","SDd13C")


cat(" The target isotope data is called geese1demo and \n")
cat(" has the following format: \n")
print(geese1demo)
cat("Press <Enter> to continue... \n")
readline()
invisible()
cat("\n The source isotope data is called sourcesdemo and looks like this: \n")
print(sourcesdemo)
cat("Press <Enter> to continue... \n")
readline()
invisible()
cat("\n The fraction correction data is called correctionsdemo: \n")
print(correctionsdemo)
cat("Press <Enter> to continue... \n")
readline()
invisible()
cat("\n The concentration depdendence data is called concdepdemo: \n")
print(concdepdemo)
cat("Press <Enter> to continue... \n")
readline()
invisible()
cat("\n")

cat("The data can be loaded by typing siarmenu() at the command prompt \n")
cat("and following the options to load in data (option 1), then load \n")
cat("in R objects and follow the instructions. \n \n")
cat("Option 2 runs SIAR for a single group, \n")
cat("This will run for ~10 seconds ... \n")

cat("Press <Enter> to run the model \n")
readline()
invisible()

out <- siarmcmcdirichletv4(geese1demo,sourcesdemo,correctionsdemo,concdepdemo)

cat("Press <Enter> to continue \n")
readline()
invisible()

cat("The data can be plotted and looks like this... \n")
cat("Press <Enter> to continue \n")
readline()
invisible()

siarplotdata(out)

cat("Press <Enter> to continue \n")
readline()
invisible()

cat("From the options menu you can now choose a plot, such as this \n")
cat("density plot... \n")

cat("Press <Enter> to see the plot \n")
readline()
invisible()

siarhistograms(out)

cat("Press <Enter> to continue")
readline()
invisible()

cat("With more complicated data sets (see geese2demo), you can fit the model \n")
cat("to multiple groups and produce different types of plots \n")
cat("For advanced users, the function siarmcmcdirichletv4() will allow runs \n")
cat("with different run parameters (such as the number of iterations). \n")
cat("Type help(siarmcmcdirichletv4) for more details. \n ")
cat("See the webpage http://mathsci.ucd.ie/~parnell_a/siar.html for further help. \n \n")

cat("Good luck using the software. \n")
cat("Please report bugs to Andrew.Parnell@ucd.ie \n")

}
