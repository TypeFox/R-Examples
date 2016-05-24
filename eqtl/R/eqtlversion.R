eqtlversion <- function()
{
	u <- strsplit(library(help = eqtl)[[3]][[1]][4], " ")[[1]]
	u[length(u)]
}
