
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
# FUNCTION:            AUGMENTED DICKEY FULLER DATA TABLES:
# adfTable              Finite sample p values for the Dickey-Fuller test
# .adfPlot              Plots sample p values for the Dickey-Fuller test
# padf                  Returns probabilities  for the ADF Test given quantiles
# qadf                  Returns quantiles for the ADF Test given probabilities
################################################################################


test.adfTable = 
function()
{ 
    # Dickey-Fuller Test Tables:
    #   adfTable(trend = c("nc", "c", "ct"), statistic = c("t", "n"),
    #       includeInf = TRUE)
    
    Table = adfTable("nc", "t")
    print(Table)
    target = sum(Table$z)
    print(target)
    current = -15.11
    checkEqualsNumeric(target, current)
    
    Table = adfTable("c", "t")
    print(Table)
    target = sum(Table$z)
    print(target) 
    current = -70.55
    checkEqualsNumeric(target, current)
    
    Table = adfTable("ct", "t")
    print(Table)
    target = sum(Table$z)
    print(target)
    current = -104.68
    checkEqualsNumeric(target, current)
    
    Table = adfTable("nc", "n")
    print(Table)
    target = sum(Table$z)
    print(target)
    current = -184.07
    checkEqualsNumeric(target, current)
    
    Table = adfTable("c", "n")
    print(Table)
    target = sum(Table$z)
    print(target)
    current = -357.11 
    checkEqualsNumeric(target, current)
    
    Table = adfTable("ct", "n")
    print(Table)
    target = sum(Table$z)
    print(target)
    current = -582.60
    checkEqualsNumeric(target, current)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.adfPlot =
function()
{
    if (FALSE) {
        
        require(akima)
        
        # .adfPlot(trend = c("nc", "c", "ct"), statistic = c("t", "n"))
    
        par(mfrow = c(1, 1))
        for (trend in c("nc", "c", "ct")) {
            for (statistic in c("t", "n")) {
                .adfPlot(trend, statistic)
            }
        }
        
    }
    
    # Return Value:
    return()  
    
}


# ------------------------------------------------------------------------------


test.adfQuantiles = 
function()
{    
    if (FALSE) {
        
        require(akima)
        
        # padf(q, n.sample, trend = c("nc", "c", "ct"), statistic = c("t", "n")) 
        # qadf(p, n.sample, trend = c("nc", "c", "ct"), statistic = c("t", "n"))
        
        p = 0.984
        n.sample = 78
        for (trend in c("nc", "c", "ct")) {
            for (statistic in c("t", "n")) {
                cat(trend, statistic, ": ")
                Q = qadf(p, n.sample, trend, statistic)
                target = P = padf(Q, n.sample, trend, statistic)
                cat(c(Q, P), 
                    checkEqualsNumeric(target, current = p, tolerance = 1e-3), "\n")
            }
        }
        
    }
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.interpolationQuantiles = 
function()
{      
    if (FALSE) {
        
        require(akima)
        
        # Extrapolation: Quantiles
        adfTable()
        check = c(NA, NA, NA, -2.66, NA, 2.16, NA, NA)
        ans = NULL
        for (p in c(0.005, 0.010, 0.990, 0.995)) {
            for (n.sample in c(10, 25)) {
                Q = qadf(p, n.sample)
                ans = rbind(ans, c(p, n.sample, Q))
            }
        }
        ans = cbind(ans, check)
        ans
        checkEqualsNumeric(target = ans[, 3], current = ans[, 4])
    
        # (Extra)Interpolation: Probabilities - uses linearInterpp()
        A = padf(q = -2.70, N = 100)  # NA
        print(A)
        B = padf(q =  2.03, N = 100)  # 0.99 
        print(B)
        C = padf(q = -1.50, N =  10)  # NA
        print(C)
        D = padf(q = -1.60, N = Inf)  # 0.1064
        print(D)
        target = c(A[[1]], B[[1]], C[[1]], D[[1]])
        current = c(NA, 0.99, NA, 0.1063745)
        checkEqualsNumeric(target, current, tolerance = 1e-4)
        
    }
    
    # Return Value:
    return()    
}


################################################################################

