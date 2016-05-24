cchart.S <- function(x, type = "n", m = NULL)
{
    if(type == "n")
        qcc(x, type = "S")
    else
    {
        if(type == "c" && !is.null(m))
            qcc(x, type = "S", limits = c((sqrt(qchisq(0.00135, m - 1) / (m - 1))) * sd.S(x),(sqrt(qchisq(0.99865, m - 1) / (m - 1))) * sd.S(x)))
        else
            if(type == "c" && is.null(m))
            {
                sprintf("WARNING: The sample size m wasn't specified, so a S ??? control chart was plotted instead.")
                qcc(x, type = "S")
            }
    }
}