qqplot.t <-
function(x, dof, print = FALSE){

length.x <- length(x)
range.x  <- range(x)
seq.line <- seq(range.x[1], range.x[2], by = 0.01)

theoretical.quantiles <- qt(((1:length.x) - 0.5)/length.x, df = dof)
sample.quantiles      <- sort(x)

plot(theoretical.quantiles, sample.quantiles, xlab = "Theoretical quantiles", ylab = "Sample quantiles")
lines(seq.line, seq.line, col = "black", lwd = 2)

if(print == TRUE){

output <- list()

output$theoretical <- theoretical.quantiles
output$sample      <- sample.quantiles 

output

}

}
