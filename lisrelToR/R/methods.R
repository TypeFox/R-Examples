
print.lisrel <- function(x, ...)
{
  cat("\n",paste(x$output,collapse="\n"),"\n")
}
# 
# # Plot function from semPlot:
# plot.lisrel <- function(x, ...) 
# {
#   if (!require("semPlot")) stop("Plot method requires 'semPlot' package to be installed")
#   semPaths(x, ...)
# }
# 
# # Summary function returns RAM:
# summary.lisrel <- function(object, ...) 
# {
#   if (!require("semPlot")) stop("Summary method requires 'semPlot' package to be installed")
#   semPlotModel(object)@RAM
# }