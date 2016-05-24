SkewnessModuleCoefficient3D <- function (modules) 
{
    n = length(modules)
	mean = sum(modules) / n
    s = ModuleStandardDeviation3D(modules)
    return((n / ((n - 1) * (n - 2))) * (sum(((modules - mean)/s)^3)))
}
