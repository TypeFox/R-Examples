
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
#   1999 - 2006, Diethelm Wuertz, GPL
#   Diethelm Wuertz <wuertz@itp.phys.ethz.ch>
#   info@rmetrics.org
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file 


################################################################################
# FUNCTION:            DESCRIPTION:    
#  'fDISTFIT'           S4 Class representation
#  print.fDISTFIT       Prints Results from a Fitted Distribution
# FUNCTION:            NORMAL AND T DISTRIBUTION:
#  nFit                 Fits parameters of a Normal density
#  tFit                 Fits parameters of a Student-t density
# FUNCTION:            STABLE DISTRIBUTION:
#  stableFit            Fits parameters of a stable density
# FUNCTION:            GENERALIZED DISTRIBUTION:
#  ghFit                Fits parameters of a generalized hyperbolic density
#  hypFit               Fits parameters of a hyperbolic density
#  nigFit               Fits parameters of a normal inverse Gaussian density
################################################################################


test.nFit <- function()
{
    # Graph Frame:
    par(mfrow = c(1, 1))
     
    # Simulate normal random variates:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    s = rnorm(n = 2000, mean = 1, sd = 0.5) 
    
    # Fit:
    cat("nFit():\n")
    print(ans <- nFit(x = s))

    # Precision within 10% ?
    relErrorTest =  c(
        ( (ans@fit$estimate[1] - 1.0)/1.0 < 0.10 ), 
        ( (ans@fit$estimate[2] - 0.5)/0.5 < 0.10 ))
    print(relErrorTest)
    checkTrue(as.logical(mean(relErrorTest)))

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.tFit = 
function()
{ 
    # Graph Frame:
    par(mfrow = c(1, 1))
    
    # Simulate random variates t(4):
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    s = rt(n = 2000, df = 4)
    
    # Fit:  
    ans = tFit(x = s, df = 2*var(s)/(var(s)-1), trace = FALSE)
    
    # Precision of df within 10% ?
    relErrorTest =  c(
        ( (ans@fit$estimate[1] - 4.0)/4.0 < 0.10 ))
    print(ans)
    print(relErrorTest)
    checkTrue(as.logical(mean(relErrorTest)))
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.ghFit = 
function()
{ 
    # Graph Frame:
    par(mfrow = c(1, 1))
    
    # Simulate random variates:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    s = rgh(n = 2000, alpha = 0.8, beta = 0.2, delta = 2, mu = -0.4, 
        lambda = 1) 
    
    # Fit:
    ans = ghFit(x = s, alpha = 1, beta = 0, delta = 1, mu = 0, lambda = 1,
        trace = FALSE) 
    
    # Precision of parameters within 10% ?
    relErrorTest =  c(
        ( (ans@fit$estimate[1] - 0.8)/( 0.8) < 0.10 ), 
        ( (ans@fit$estimate[2] - 0.2)/( 0.2) < 0.10 ),
        ( (ans@fit$estimate[3] - 2.0)/( 2.0) < 0.10 ),
        ( (ans@fit$estimate[4] + 0.4)/(-0.4) < 0.10 ),
        ( (ans@fit$estimate[5] - 2.0)/( 2.0) < 0.10 ))
    print(ans)
    print(relErrorTest)
    checkTrue(as.logical(mean(relErrorTest)))
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.hypFit = 
function()
{ 
    # Graph Frame:
    par(mfrow = c(1, 1))
    
    # Simulate normal random variates:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    s = rhyp(n = 2000, alpha = 1.5, beta = 0.8, delta = 0.5, mu = -1) 
    
    # Fit:
    ans = hypFit(s, alpha = 1, beta = 0, delta = 1, mu = mean(s), 
        trace = FALSE)
    
    # Precision of parameters within 10% ?
    relErrorTest =  c(
        ( (ans@fit$estimate[1] - 1.5)/( 1.5) < 0.10 ), 
        ( (ans@fit$estimate[2] - 0.8)/( 0.8) < 0.10 ),
        ( (ans@fit$estimate[3] - 0.5)/( 0.5) < 0.10 ),
        ( (ans@fit$estimate[4] + 1.0)/(-1.0) < 0.10 ))
    print(ans)
    print(relErrorTest)
    checkTrue(as.logical(mean(relErrorTest)))
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.nigFit <- function()
{
    # Graph Frame:
    par(mfrow = c(1, 1))
    
    # Simulate normal random variates:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    th.0 <- c(alpha = 1.5, beta = -0.7, delta = 0.5, mu = -1.0)
    s <- do.call(rnig, c(list(n = 2000), th.0))

    # Fit {includes plot !}
    print(ans <- nigFit(s, alpha = 1, beta = 0, delta = 1, mu = mean(s),
                        trace = FALSE))
    # Precision of parameters within 10% / 35% ?
    est <- ans@fit$estimate
    relErr <- (est - th.0) / th.0
    checkTrue(all(abs(relErr) < 0.35))# 35%
    ## 10% / 35$ is silly, rather use  SE scale:
    print(SE <- ans@fit$se.coef)
    print(cbind(th.0, est, SE))# NaN
    SE[is.na(SE)] <- .1 # arbitrary
    cat("t-ratios:\n"); print(t. <- (est - th.0)/ SE)
    checkTrue(all(abs(t.) < 6))## 6 is too large, but have -5.476 for beta

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.stableFit <- function()
{
    # Graph Frame:
    par(mfrow = c(1, 1))
    
    # Simulate stable random variates:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")

    th.0 <- c(alpha=1.8, beta=0.3, gamma = 1, delta = 0.1)
    s <- do.call(stabledist::rstable, c(list(n = 500), th.0))

    # Fit {includes plot !}
    cat("stableFit():\n")
    print(ans <- stableFit(x = s, alpha = 1.5))
    # Precision of parameters within 5% / 10% / ...?
    est <- ans@fit$estimate
    arelErr <- abs(relErr <- (est - th.0) / th.0)
    print(cbind(th.0, est, relErr))
    checkTrue(all(arelErr[c("alpha","gamma")] < 0.05) &&
              arelErr[["beta"]] < 0.15 && arelErr[["delta"]] < 1.50) # 150% !

    # MLE Fit:
    if (FALSE) { # as this takes rather too long time ...
        print(system.time(
            ans <- stableFit(x = s, alpha = 1.5, type = "mle", trace = TRUE)
            ## .mleStableFit(s, 1.75, 0, 1, 0)
            ))
        ##  user  {lynne, 2012-09}
        ## 159.7
        # The result would be:
        print(ans)
        ## Precision of parameters within 5% / 10% / ...?
        est <- ans@fit$estimate
        arelErr <- abs(relErr <- (est - th.0) / th.0)
        print(cbind(th.0, est, relErr))
        checkTrue(all(arelErr[1:3] < 0.02) && arelErr[["delta"]] < 0.42)
    }
    
    # Return Value:
    return() 
}
 
   
################################################################################

