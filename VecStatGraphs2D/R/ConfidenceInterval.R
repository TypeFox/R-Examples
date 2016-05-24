ConfidenceInterval <- function (n, mean_azimuth, mean_module, von_mises) 
{
    z = 1.96
    vm = 1/(sqrt(n * von_mises * mean_module))
    ci1 = mean_azimuth + ToSexagesimal(asin(z * vm))
    ci2 = mean_azimuth - ToSexagesimal(asin(z * vm))
    return(c(ci1, ci2))
}
