cchart.p <- function(x1 = NULL, n1 = NULL, type = "norm", p1 = NULL, x2 = NULL, n2 = NULL, phat = NULL, p2 = NULL)
{
    if((!is.null(n1)) && (!is.null(x1) || !is.null(p1)))
        OK1 = TRUE
    else
        OK1 = FALSE
    if(!is.null(n2) && (!is.null(x2) || !is.null(p2)) && (OK1 || !is.null(phat)))
        OK2 = TRUE
    else
        OK2 = FALSE
    
#-- Error messages
    if(!OK1 && !OK2)
    {
        if(is.null(x1) && is.null(n1) && is.null(p1))
            return("Phase I data and samples sizes are missing")
        else
        {
            if(is.null(n1))
                return("Phase I samples sizes not specified")
            else
                return("Phase I data is missing")
        }
    }
    if(!OK2)
    {
        if(is.null(n2) && (!is.null(x2) || !is.null(p2)))
            return("Phase II samples sizes not specified")
        if(!is.null(n2) && (is.null(x2) && is.null(p2)))
            return("Phase II data is missing")
        if(!is.null(x2) && !is.null(n2) && !is.null(p2))
            return("Information about phase I is missing")
    }

#-- Phase I
    if(OK1 && !OK2)
    {
        if(!is.null(x1))
        {
            m1 <- length(x1)
            if(length(n1) != length(x1))
                return("The arguments x1 and n1 must have the same length")
        }
        if(!is.null(p1))
        {
            m1 <- length(p1)
            if(length(n1) != length(p1))
                return("The arguments p1 and n1 must have the same length")
        }
        if(is.null(p1))
            p1 <- x1 / n1
        if(is.null(x1))
            x1 <- p1 * n1
        phat <- mean(p1)
        l <- matrix(nrow = m1, ncol = 1)
#------ Shewhart
        if(type == "norm")
        {
            u <- matrix(nrow = m1, ncol = 1)
            for(i in 1:m1)
            {
                UCL <- phat + (3 * sqrt((phat * (1 - phat)) / n1[i]))
                u[i, ] <- UCL
                LCL <- phat - (3 * sqrt((phat * (1 - phat)) / n1[i]))
                l[i, ] <- LCL
            }
            qcc(x1, type = "p", n1, limits = c(l, u), center = phat, title = "Shewhart p-chart (phase I)")
        }
#------ Cornish-Fisher
        if(type == "CF")
        {
            u <- matrix(nrow = m1, ncol = 1)
            for(i in 1:m1)
            {
                UCL <- phat + (3 * sqrt((phat * (1 - phat)) / n1[i])) + (4 * (1 - 2 * phat) / (3 * n1[i]))
                u[i, ] <- UCL
                LCL <- phat - (3 * sqrt((phat * (1 - phat)) / n1[i])) + (4 * (1 - 2 * phat) / (3 * n1[i]))
                l[i, ] <- LCL
            }
            qcc(x1, type = "p", n1, limits = c(l, u), center = phat, title = "Cornish-Fisher p-chart (phase I)")
        }
#------ Standardized
        if(type == "std")
        {
            for(i in 1:m1)
            {
                z <- (p1[i] - phat) / sqrt((phat * (1 - phat)) / n1[i])
                l[i, ] <- z
            }
            std <- l * n1
            qcc(std, type = "p", n1, center = 0, limits = c(-3, 3), title = "Standardized p-chart (phase I)")
        }
    }
    
#-- Phase II
    if(OK2)
    {
        if(!is.null(x2))
        {
            m2 <- length(x2)
            if(length(n2) != length(x2))
                return("The arguments x2 and n2 must have the same length")
        }
        if(!is.null(p2))
        {
            m2 <- length(p2)
            if(length(n2) != length(p2))
                return("The arguments p2 and n2 must have the same length")
        }
        if(is.null(p2))
            p2 <- x2 / n2
        if(is.null(x2))
            x2 <- p2 * n2
        if(is.null(phat))
        {
            if(is.null(p1))
                p1 <- x1 / n1
            phat <- mean(p1)
        }
        l <- matrix(nrow = m2, ncol = 1)
#------ Shewhart
        if(type == "norm")
        {
            u <- matrix(nrow = m2, ncol = 1)
            for(i in 1:m2)
            {
                UCL <- phat + (3 * sqrt((phat * (1 - phat)) / n2[i]))
                u[i, ] <- UCL
                LCL <- phat - (3 * sqrt((phat * (1 - phat)) / n2[i]))
                l[i, ] <- LCL
            }
            qcc(x2, type = "p", n2, limits = c(l, u), center = phat, title = "Shewhart p-chart (phase II)")
        }
#------ Cornish-Fisher
        if(type == "CF")
        {
            u <- matrix(nrow = m2, ncol = 1)
            for(i in 1:m2)
            {
                UCL <- phat + (3 * sqrt((phat * (1 - phat)) / n2[i])) + (4 * (1 - 2 * phat) / (3 * n2[i]))
                u[i, ] <- UCL
                LCL <- phat - (3 * sqrt((phat * (1 - phat)) / n2[i])) + (4 * (1 - 2 * phat) / (3 * n2[i]))
                l[i, ] <- LCL
            }
            qcc(x2, type = "p", n2, limits = c(l, u), center = phat, title = "Cornish-Fisher p-chart (phase II)")
        }
#------ Standardized
        if(type == "std")
        {
            for(i in 1:m2)
            {
                z <- (p2[i] - phat) / sqrt((phat * (1 - phat)) / n2[i])
                l[i, ] <- z
            }
            std <- l * n2
            qcc(std, type = "p", n2, center = 0, limits = c(-3, 3), title = "Standardized p-chart (phase II)")
        }
    }
}
