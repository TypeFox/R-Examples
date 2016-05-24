RayleighTest <- function (azimuths) 
{
    n = length(azimuths)
    m_module = MeanModule(azimuths)
    z = n * (m_module^2)
    p = exp(-z)
    if (round(p, 3) == 0) {
        print(paste("Rayleigh Test: P-value for the hypothesis of uniformity < 0.001"))
    }
    else {
        print(paste("Rayleigh Test: P-value for the hypothesis of uniformity =", 
            round(p, 3)))
    }
}
