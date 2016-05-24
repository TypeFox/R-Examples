`print.ecespa.marksum` <-
function(x,...)
{
cat("Mark-sum measure of the dataset", x$dataname,  "\n computed for a radius R=",x$R, 
", with a grid of", x$nx, "x",x$ny,".\n")
cat("Plot it to see the result.\n")
}

