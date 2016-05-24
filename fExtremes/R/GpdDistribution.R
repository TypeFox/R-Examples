
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the 
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General 
# Public License along with this library; if not, write to the 
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, 
# MA  02111-1307  USA

# Copyrights (C)
# for this R-port: 
#   1999 - 2007, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# FUNCTION:               GPD DISTRIBUTION FAMILY:
#  dgpd                    Density for the Generalized Pareto DF [USE FROM EVIS]
#   pgpd                    Probability for the Generalized Pareto DF
#   qgpd                    Quantiles for the Generalized Pareto DF
#   rgpd                    Random variates for the Generalized Pareto DF
# FUNCTION:               GPD MOMENTS:
#  gpdMoments              Computes true statistics for GPD distribution
# FUNCTION:               GPD SLIDER:
#  gpdSlider               Displays distribution and rvs for GPD distribution
# FUNCTION:               INTERNAL GPD DISTRIBUTION FAMILY:
#  .depd                   Density for the Generalized Pareto DF 
#  .pepd                   Probability for the Generalized Pareto DF
#  .qepd                   Quantiles for the Generalized Pareto DF
#  .repd                   Random variates for the Generalized Pareto DF
################################################################################


dgpd <- 
    function(x, xi = 1, mu = 0, beta = 1, log = FALSE)
{   
    # A function written by Diethelm Wuertz

    # Description:
    #   Density for the Generalized Pareto DF
    
    # Arguments:
    
    # FUNCTION:
    
    # Transform:
    shape = xi
    location = mu
    scale = beta
    
    # Density:
    d = .depd(x, location, scale, shape, log)
    
    # Add Control:
    attr(d, "control") = data.frame(xi = xi, mu = mu, 
        beta = beta[1], log = log, row.names = "")
        
    # Return Value:
    d
}


# ------------------------------------------------------------------------------


pgpd <- 
    function(q, xi = 1, mu = 0, beta = 1, lower.tail = TRUE)
{   
    # A function written by Diethelm Wuertz

    # Description:
    #   Probability for the Generalized Pareto DF
    
    # Arguments:
    
    # FUNCTION:
    
    # Transform:
    shape = xi
    location = mu
    scale = beta
    
    # Probability:
    p = .pepd(q, location, scale, shape, lower.tail)
    
    # Add Control:
    attr(p, "control") = data.frame(xi = xi, mu = mu, 
        beta = beta[1], lower.tail = lower.tail, row.names = "")
        
    # Return Value:
    p
}


# ------------------------------------------------------------------------------


qgpd <- 
    function(p, xi = 1, mu = 0, beta = 1, lower.tail = TRUE)
{   
    # A function written by Diethelm Wuertz

    # Description:
    #   Quantiles for the Generalized Pareto DF
    
    # Arguments:
    
    # FUNCTION:
    
    # Transform:
    shape = xi
    location = mu
    scale = beta
    
    # Quantiles:
    q = .qepd(p, location, scale, shape, lower.tail)
    
    # Add Control:
    attr(q, "control") = data.frame(xi = xi, mu = mu, 
        beta = beta[1], lower.tail = lower.tail, row.names = "")
        
    # Return Value:
    q
}


# ------------------------------------------------------------------------------


rgpd <- 
    function(n, xi = 1, mu = 0, beta = 1)
{   
    # A function written by Diethelm Wuertz

    # Description:
    #   Random variates for the Generalized Pareto DF
    
    # Arguments:
    
    # FUNCTION:
    
    # Transform:
    shape = xi
    location = mu
    scale = beta
    
    # Random Variates:
    r = .repd(n, location, scale, shape)
    
    # Add Control:
    attr(r, "control") = data.frame(xi = xi, mu = mu, 
        beta = beta[1], row.names = "")
        
    # Return Value:
    r
}


# ------------------------------------------------------------------------------


gpdMoments <- 
    function(xi = 1, mu = 0, beta = 1)
{   
    # A function implemented by Diethelm Wuertz
 
    # Description:
    #   Compute true statistics for Generalized Pareto distribution
    
    # Arguments:
    
    # Value:
    #   Returns true mean of Generalized Pareto distribution 
    #   for xi < 1 else NaN
    #   Returns true variance of Generalized Pareto distribution 
    #   for xi < 1 else NaN

    # FUNCTION: 
    
    # MEAN: Returns 1 for x <= 0 and -Inf's's else
    a = c(1, NaN, NaN)
    gpdMean = mu + beta/(1-xi)*a[sign(xi-1)+2]
    
    # VAR: Rreturns 1 for x <= 0 and -Inf's's else
    a = c(1, NaN, NaN)
    gpdVar = beta*beta/(1-xi)^2/(1-2*xi) * a[sign(2*xi-1)+2]
    
    # Result:
    param = c(xi = xi, mu = mu, beta = beta)
    ans = list(param = param, mean = gpdMean, var = gpdVar)      
    
    # Return Value:
    ans         
}


# ------------------------------------------------------------------------------


gpdSlider <- 
    function(method = c("dist", "rvs"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays distribution and rvs for GPD distribution
    
    # Arguments:
    
    # FUNCTION:
    
    # Settings:
    method = match.arg(method)
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        N = .sliderMenu(no = 1)
        xi = .sliderMenu(no = 2)
        mu = .sliderMenu(no = 3)
        beta = .sliderMenu(no = 4)
        
        # Compute Data:  
        pmin = 0.00  
        pmax = 0.99 
        xmin = round(qgpd(pmin, xi, mu, beta), digits = 2)
        xmax = round(qgpd(pmax, xi, mu, beta), digits = 2)
        s = seq(xmin, xmax, length = N)
        y1 = dgpd(s, xi, mu, beta)
        y2 = pgpd(s, xi, mu, beta)
        Moments = gpdMoments(xi, mu, beta)
        Mean = round(Moments$mean, 2)
        Var = round(Moments$var, 2)
        mText = paste("Mean =", Mean, " | Variance = ", Var)
        main1 = paste("GPD Density\n", 
            "xi = ", as.character(xi), " | ",
            "mu = ", as.character(mu), " | ",
            "beta = ", as.character(beta) )
        main2 = paste("GPD Probability\n",
            "xmin [0.00] = ", as.character(xmin), " | ",
            "xmax [0.99] = ", as.character(xmax) )  
        Median = qgpd(0.5, xi, mu, beta)   
            
        # Frame:
        par(mfrow = c(2, 1), cex = 0.7)
        
        # Density:
        if (method == "rvs") {
            x = rgpd(N, xi, mu, beta)
            hist(x, probability = TRUE, col = "steelblue", border = "white",
                xlim = c(xmin, xmax), ylim = c(0, 1.1*max(y1)), main = main1,
                breaks = "FD" )
            lines(s, y1, col = "orange")
            mtext(mText, side = 4, col = "grey", cex = 0.7)  
        } else {
            plot(s, y1, type = "l", xlim = c(xmin, xmax), col = "steelblue")
            abline(h = 0, lty = 3)
            abline(v = Median, lty = 3, col = "red")
            abline(v = Mean, lty = 3, col = "darkgreen")
            title(main = main1)  
            mtext(mText, side = 4, col = "grey", cex = 0.7)  
        }    
        
        # Probability:
        plot(s, y2, type = "l", xlim = c(xmin, xmax), ylim = c(0, 1),
            col = "steelblue" )
        abline(h = 0, lty = 3)
        abline(h = 0.5, lty = 3, col = "red")
        abline(v = Median, lty = 3, col = "red")
        abline(v = Mean, lty = 3, col = "darkgreen")
        title(main = main2) 
        mtext(mText, side = 4, col = "grey", cex = 0.7)  
        
        # Reset Frame:
        par(mfrow = c(1, 1), cex = 0.7)
    }
  
    # Open Slider Menu:
    .sliderMenu(refresh.code,
       names =       c(   "N", "xi",  "mu",  "beta"),
       minima =      c(   50,  0.00, -5.00,   0.10 ),
       maxima =      c( 1000,  1.50, +5.00,   5.00 ),
       resolutions = c(   50,  0.01,  0.10,   0.10 ),
       starts =      c(  500,  1.00,  0.00,   1.00 )
    )
}


################################################################################


.depd <- 
    function(x, location = 0, scale = 1, shape = 0, log = FALSE) 
{
    # Description:
    
    # Arguments:
    
    # FUNCTION:
    
    # Check:
    stopifnot(min(scale) > 0) 
    stopifnot(length(shape) == 1) 
    
    # Density:
    d <- (x - location)/scale
    nn <- length(d)
    scale <- rep(scale, length.out = nn)
    index <- (d > 0 & ((1 + shape * d) > 0)) | is.na(d)
    if (shape == 0) {
        d[index] <- log(1/scale[index]) - d[index]
        d[!index] <- -Inf
    } else {
        d[index] <- log(1/scale[index]) - (1/shape+1)*log(1+shape*d[index])
        d[!index] <- -Inf
    }
    
    # Log:
    if (!log) 
        d <- exp(d)
    
    # Add Control:
    attr(d, "control") = data.frame(location = location[1], scale = scale[1], 
        shape = shape[1], log = log, row.names = "")
           
    # Return Value: 
    d
}


# ------------------------------------------------------------------------------


.pepd <- 
    function(q, location = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{
    # Description:
    
    # Arguments:
    
    # FUNCTION:
    
    # Check:
    stopifnot(min(scale) > 0) 
    stopifnot(length(shape) == 1) 
        
    # Probability:
    q <- pmax(q - location, 0)/scale
    if (shape == 0) 
        p <- 1 - exp(-q)
    else {
        p <- pmax(1 + shape * q, 0)
        p <- 1 - p^(-1/shape)
    }
    
    # Lower Tail:
    if (!lower.tail) 
        p <- 1 - p
       
    # Add Control:
    attr(p, "control") = data.frame(location = location[1], scale = scale[1], 
        shape = shape[1], lower.tail = lower.tail, row.names = "")
        
    # Return Value: 
    p
}


# ------------------------------------------------------------------------------


.qepd <- 
    function(p, location = 0, scale = 1, shape = 0, lower.tail = TRUE) 
{ 
    # Description:
    
    # Arguments:
    
    # FUNCTION:
    
    # Check:
    stopifnot(min(scale) > 0) 
    stopifnot(length(shape) == 1)
    stopifnot(min(p, na.rm = TRUE) >= 0)
    stopifnot(max(p, na.rm = TRUE) <= 1) 
        
    # Lower Tail:
    if (lower.tail) 
        p <- 1 - p
        
    # Quantiles:
    if (shape == 0) {
        q = location - scale * log(p)
    } else {
        q = location + scale * (p^(-shape) - 1)/shape
    }
       
    # Add Control:
    attr(q, "control") = data.frame(location = location[1], scale = scale[1], 
        shape = shape[1], lower.tail = lower.tail, row.names = "")
        
    # Return Value: 
    q
}


# ------------------------------------------------------------------------------


.repd <- 
    function(n, location = 0, scale = 1, shape = 0) 
{
    # Description:
    
    # Arguments:
    
    # FUNCTION:
    
    # Check:
    stopifnot(min(scale) > 0) 
    stopifnot(length(shape) == 1)
        
    # Random Variates:
    if (shape == 0) {
        r = location + scale * rexp(n)
    } else {
        r = location + scale * (runif(n)^(-shape) - 1)/shape
    }
    
    # Add Control:
    attr(r, "control") = data.frame(location = location[1], scale = scale[1], 
        shape = shape[1], row.names = "")
    
    # Return Value:
    r
}


################################################################################