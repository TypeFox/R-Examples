cchart.Xbar <- function(x1 = NULL, n1 = NULL, x2 = NULL, n2 = NULL, x2bars = NULL, sigma = NULL)
{
    if(!is.null(x1) && !is.null(n1))
        OK1 = TRUE
    else
        OK1 = FALSE
    if(!is.null(x2) && !is.null(n2) && (OK1 || (!is.null(x2bars) && !is.null(sigma))))
        OK2 = TRUE
    else
        OK2 = FALSE
#-- Error messages
    if(!OK1 && !OK2)
    {
        if(is.null(x1) && is.null(n1))
            return("Phase I data and samples sizes are missing")
        else
        {
            if(!is.null(n1))
                return("Phase I data is missing")
            else
                return("Phase I samples sizes not specified")
        }
    }
    if(!OK2)
    {
        if(is.null(x2) && !is.null(n2))
            return("Phase II data is missing")
        if(!is.null(x2) && is.null(n2))
            return("Phase II samples sizes not specified")
    }

#-- Phase I
    if(OK1 && !OK2)
    {
        a <- rowMeans(x1)
        x2bars <- mean(a)
        sigma <- sd.xbar(x1)
        qcc(x1, type = "xbar", n1)
    }
#-- Phase II
    if(OK2)
    {
        if(is.null(x2bars))
        {
            a <- rowMeans(x1)
            x2bars <- mean(a)
        }
        if(is.null(sigma))
            sigma <- sd.xbar(x1)
        qcc(x2, type = "xbar", center = x2bars, std.dev = sigma)
    }
}