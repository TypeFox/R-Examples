`print.ecespa.getis` <-
function(x,...)
{

cat("Getis local density function of the dataset", x$dataname,  "\n computed for a radius R=",x$R, 
", with a grid of", x$nx, "x",x$ny,".\n")
     
cat("Plot it to see the result.\n")
}

