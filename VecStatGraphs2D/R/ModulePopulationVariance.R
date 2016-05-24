ModulePopulationVariance <- function (modules) 
{
    n = length(modules)
	mean = sum(modules) / n
    return((sum((modules - mean)^2)) / n)
}
