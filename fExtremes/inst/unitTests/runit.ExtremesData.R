
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
# FUNCTION             EXPLORATIVE DATA ANALYSIS:
#  emdPlot              Creates an empirical distribution plot
#  qqparetoPlot         Creates exploratory QQ plot for EV analysis
#  mePlot               Creates a sample mean excess function plot
#   mxfPlot             Creates another view of a sample mean excess plot
#   mrlPlot             Returns a mean residual life plot with confidence levels
#  recordsPlot          Plots records development
#   ssrecordsPlot       Plots records development of data subsamples
#  msratioPlot          Plots ratio of maximums and sums
#  sllnPlot             Verifies Kolmogorov's Strong Law of large numbers
#  lilPlot              Verifies Hartman-Wintner's Law of the iterated logarithm
#  xacfPlot             Plots autocorrelations of exceedences
################################################################################


test.emd = 
function()
{
    # emdPlot - Creates an empirical distribution plot
    
    # Artificial Data Set:
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    x = rgpd(1000)
    # Empirical distribution plot:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    emdPlot(x, plottype = "xy")
    emdPlot(x, plottype = "x")
    emdPlot(x, plottype = "y")
    # emdPlot(x, plottype = " ")                # CHECK !!!
    
    # Artificial Data Set:
    x = rt(1000, df = 3) 
    # Empirical distribution plot:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    emdPlot(x, plottype = "xy")
    emdPlot(x, plottype = "x")
    emdPlot(x, plottype = "y")
    # emdPlot(x, plottype = " ")                # CHECK !!!
    
    # Artificial Data Set:
    x = rnorm(1000) 
    # Empirical distribution plot:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    emdPlot(x, plottype = "xy")
    emdPlot(x, plottype = "x")
    emdPlot(x, plottype = "y")
    # emdPlot(x, plottype = " ")                # CHECK !!!
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.qqpareto = 
function()
{
    # qqparetoPlot - Creates exploratory QQ plot for EV analysis

    # Artificial Data Set - 
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    r0 = rgpd(n = 1000, xi = 0)
    r1 = rgpd(n = 1000, xi = 1)
    
    # Graph Frame:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
       
    # Empirical Pareto Distribution Plot:
    qqparetoPlot(x = r0, xi = 0)
    qqparetoPlot(x = r1, xi = 1)
    
    # Empirical Normal Distribution Plot:
    qqnormPlot(x = r0)
    qqnormPlot(x = r1)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.me = 
function()
{
    # mePlot - Creates a sample mean excess function plot
    # mxfPlot - Creates another view of a sample mean excess plot
    # mrlPlot - Returns a mean residual life plot with confidence levels
    
    # Artificial Data Set - 
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    r = rgpd(n = 1000)
     
    # Mean Excess Function Plot:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    mePlot(x = r)           # Check, the largest point is missing ...
    mxfPlot(x = r)
    mrlPlot(x = r)
    
    # No Labels:
    par(mfrow = c(2, 2))
    par(ask = FALSE)
    mePlot(x = r, labels = FALSE)
    mxfPlot(x = r, labels = FALSE)
    mrlPlot(x = r, labels = FALSE)
      
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.records = 
function()
{
    #  recordsPlot - Plots records development
    #  ssrecordsPlot - Plots records development of data subsamples

    # Artificial Data Set - 
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    r = rgpd(n = 1000)
    
    # Records Plot:
    par(mfrow = c(2, 2))
    par(ask = FALSE)
    recordsPlot(x = r)
    recordsPlot(x = r, ci = 0.99)
    ans = recordsPlot(x = r, labels = FALSE)
    print(ans)
    
    # Subrecords Plot:
    set.seed(1985)
    r = rgpd(n = 10000)
    par(mfrow = c(2, 2))
    par(ask = FALSE)
    recordsPlot(r)
    ssrecordsPlot(r, subsamples = 1)
    ssrecordsPlot(r, subsamples = 1, plottype = "log")
    ans = ssrecordsPlot(r, subsamples = 1, plottype = "lin")
    print(ans)
    
    # Subrecords Plot:
    set.seed(1985)
    r = rgpd(n = 10000)
    par(mfrow = c(2, 2))
    par(ask = FALSE)
    ssrecordsPlot(r, subsamples = 10)
    ssrecordsPlot(r, subsamples = 50)
    ssrecordsPlot(r, subsamples = 10, plottype = "log")
    ans = ssrecordsPlot(r, subsamples = 50, plottype = "log", labels = FALSE)
    print(ans)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.msratio = 
function()
{
    # msratioPlot - Plots ratio of maximums and sums
    
    # Artificial Data Set - 
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry") 
    r = rgpd(n = 1000)
    
    # Mean Excess Function Plot:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    msratioPlot(x = r, p = 1:4)
    ans = msratioPlot(x = r, p = 1:4, labels = FALSE)
    print(head(ans))
     
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.laws = 
function()
{
    # sllnPlot - Verifies Kolmogorov's Strong Law of large numbers
    # lilPlot - Verifies Hartman-Wintner's Law of the iterated logarithm

    # Artificial Data Set - 
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    r = rgpd(n = 1000)
    
    # Strong Law of Large Numbers:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    sllnPlot(x = r)
    ans = sllnPlot(x = r, labels = FALSE)
    print(ans)
    
    # Law of the Iterated Logarithm:
    lilPlot(x = r)
    ans = lilPlot(x = r, labels = FALSE)
    print(ans)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.xacf = 
function()
{
    # xacfPlot - Plots autocorrelations of exceedences
    
    # Create an Artificial Data Set: 
    RNGkind(kind = "Marsaglia-Multicarry", normal.kind = "Inversion")
    set.seed(4711, kind = "Marsaglia-Multicarry")
    r = rgpd(n = 1000)
    
    # ACF of Exceedances Plot:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    ans = xacfPlot(x = r)
    print(ans)
    
    # ACF of Exceedances Plot:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    xacfPlot(x = r, labels = FALSE)
    
    # ACF of Exceedances Plot:
    par(ask = FALSE)
    par(mfrow = c(2, 2))
    xacfPlot(x = r, labels = FALSE, which = 1); title(main = "1")
    xacfPlot(x = r, labels = FALSE, which = 2); title(main = "2")
    xacfPlot(x = r, labels = FALSE, which = "3"); title(main = "3")
    xacfPlot(x = r, labels = FALSE, which = "4"); title(main = "4")
      
    # Return Value:
    return()    
}


################################################################################

