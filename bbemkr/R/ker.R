ker <- function(u, kerntype = c("Gaussian", "Epanechnikov", "Quartic", "Triweight", "Triangular", "Uniform"))
{
    kerntype = match.arg(kerntype)
    if(kerntype == "Gaussian")
    {
        result = 1/(sqrt(2 * pi)) * exp(-0.5 * (u^2))
    }
    else
    {
        lenu = length(u)
        result = vector(, lenu)
        for(j in 1:lenu)
        {
            if(abs(u[j])<=1)
            {
                if(kerntype == "Epanechnikov")
                {
                    result[j] = 3/4 * (1 - u[j]^2)
                }
                if(kerntype == "Quartic")
                {
                    result[j] = 15/16 * ((1 - u[j]^2)^2)
                }
                if(kerntype == "Triweight")
                {
                    result[j] = 35/32 * ((1 - u[j]^2)^3)
                }
                if(kerntype == "Triangular")
                {
                    result[j] = (1 - abs(u[j]))
                }
                if(kerntype == "Uniform")
                {
                    result[j] = 1/2
                }
            }
            else
            {
                result[j] = 0
            }
        }
    }
    return(result)
}

