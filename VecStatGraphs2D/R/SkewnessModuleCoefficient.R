SkewnessModuleCoefficient <- function (modules) 
{
    n = length(modules)
	mean = sum(modules) / n
    s = ModuleStandardDeviation(modules)
    return((n / ((n - 1) * (n - 2))) * (sum(((modules - mean)/s)^3)))
}
