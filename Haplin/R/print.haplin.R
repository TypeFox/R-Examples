print.haplin <- function(x,...){
cat('This is the result of a haplin run.\n')
cat("Number of data lines used: ", round(x$result$ntri), " | Number of haplotypes used: ", round(x$result$nall), "\n", sep = "")
#invisible(print.tri.glm(x$result,...))
cat('Please use the "summary", "plot", "haptable" or "output" functions to obtain more details.\n')
}
