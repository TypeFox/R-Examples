`print.ecespa.kmm` <-
function(x, ...)
{
    if(!is.null(x$nsim)){
cat("mark-weighted K function computed on",x$dataname,"\nand tested with",
               x$nsim,"random permutations of marks.\n")
    }
    if(is.null(x$nsim)){
cat("mark-weighted K function computed on",x$dataname,"\n")
   }
cat("plot it to see the result.\n")
}

