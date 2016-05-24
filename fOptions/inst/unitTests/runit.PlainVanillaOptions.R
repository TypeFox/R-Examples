
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
# FUNCTION:                  DESCRIPTION:
#  'fOPTION'                  S4 Class Representation
# FUNCTION:                  DESCRIPTION:
#  NDF                        Normal distribution function
#  CND                        Cumulative normal distribution function
#  CBND                       Cumulative bivariate normal distribution  
# FUNCTION:                  DESCRIPTION:
#  GBSOption                  Computes Option Price from the GBS Formula
#  GBSCharacteristics         Computes Option Price and all Greeks of GBS Model
#  BlackScholesOption         Synonyme Function Call to GBSOption
#  GBSGreeks                  Computes one of the Greeks of the GBS formula
# FUNCTION:                  DESCRIPTION:
#  Black76Option              Computes Prices of Options on Futures
#  MiltersenSchwartzOption    Pricing a Miltersen Schwartz Option
# S3 METHODS:                DESCRIPTION:
#  print.option               Print Method
#  summary.otion              Summary Method
################################################################################


test.NDF =
function()
{
    # NDF:                        
    #   Normal distribution function
    
    # Arguments:
    #   NDF(x)
    
    # NDF:
    x = (-3):3
    NDF(x)
    dnorm(x)
    NDF(x)-dnorm(x)

    # Return Value:
    return()    
}

# ------------------------------------------------------------------------------


test.CND =
function()
{
    # CND:                        
    #   Cumulative normal distribution function 
    
    # Arguments:
    #   CND(x)
    
    # CND:
    # NDF:
    x = (-3):3
    CND(x)
    pnorm(x)
    CND(x)-pnorm(x)

    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.CBND =
function()
{
    # CBND:                   
    #   Cumulative bivariate normal distribution  
    
    # Arguments:
    #   CBND(x1, x2, rho)
    
    # CBND:
    CBND(0, 0, 1/2)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------

test.GBSOption =
function()
{
    # GBSOption:                   
    #   Computes Option Price from the GBS Formula

    # Arguments:
    #   GBSOption(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, 
    #   title = NULL, description = NULL) 
    
    # GBSOption:
    GBSOption("c", 100, 100, 1, 0.10, 0.10, 0.30)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.GBSCharacteristics =
function()
{
    # GBSCharacteristics:          
    #   Computes Option Price and all Greeks of GBS Model

    # Arguments:
    #   GBSCharacteristics(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma) 
    
    # GBSCharacteristics:
    GBSCharacteristics("c", 100, 100, 1, 0.10, 0.10, 0.30)
    
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.BlackScholesOption =
function()
{
    # BlackScholesOption:        
    #   Synonyme Function Call to GBSOption
   
    # Arguments:
    #   BlackScholesOption(...)
  
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.Black76Option =
function()
{ 
    # Black76Option              
    #   Computes Prices of Options on Futures
    
    # Arguments:
    #   Black76Option = (TypeFlag = c("c", "p"), FT, X, Time, r, sigma, 
    #   title = NULL, description = NULL)
    
    # Black76Option:
    Black76Option(FT = 95, X = 80, Time = 1/2, r = 0.05, sigma = 0.266)
        
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.MiltersenSchwartzOption =
function()
{ 
    # MiltersenSchwartzOption    
    #   Pricing a Miltersen Schwartz Option

    # Arguments:
    #   MiltersenSchwartzOption(TypeFlag = c("c", "p"), Pt, FT, X, time, Time, 
    #   sigmaS, sigmaE, sigmaF, rhoSE, rhoSF, rhoEF, KappaE, KappaF, 
    #   title = NULL, description = NULL)
    
    # MiltersenSchwartzOption:
    MiltersenSchwartzOption(TypeFlag = "c", Pt = exp(-0.05/4), FT = 95, 
        X = 80, time = 1/4, Time = 1/2, sigmaS = 0.2660, sigmaE = 0.2490, 
        sigmaF = 0.0096, rhoSE = 0.805, rhoSF = 0.0805, rhoEF = 0.1243, 
        KappaE = 1.045, KappaF = 0.200)
        
    # Return Value:
    return()    
}


# ------------------------------------------------------------------------------


test.print =
function()
{
    # GBSOption(TypeFlag = c("c", "p"), S, X, Time, r, b, sigma, 
    #   title = NULL, description = NULL) 
    GBS = GBSOption("c", 100, 100, 1, 0.10, 0.10, 0.30)
    
    # Print Method:
    show(GBS)
    print(GBS)
    
    # Summary Method:
    summary(GBS)
    
    # Return Value:
    return()    
}

        
# ------------------------------------------------------------------------------
    
    
test.GBSOptionSlider = 
function()
{  
    .GBSOptionSlider = 
    function(TypeFlag = "c", S = 100, X = 100, Time = 1, r = 0.10, b = 0.10, 
    sigma = 0.25, span = 0.25, N = 40) 
    {
        # Internal Function:
        refresh.code = function(...)
        {
            # Sliders:
            S     = .sliderMenu(no = 1)
            X     = .sliderMenu(no = 2)
            Time  = .sliderMenu(no = 3)
            sigma = .sliderMenu(no = 4) 
            r     = .sliderMenu(no = 5) 
            b     = .sliderMenu(no = 6) 
            theta = .sliderMenu(no = 7)
            phi   = .sliderMenu(no = 8)
            
            TypeFlagText = c(c = "Call:", p = "Put:")
            
            if (r != rNow | b != bNow) {
                for (j in 1:nY)
                    z[j, ] <<- GBSOption(TypeFlag, sOption, xOption, 
                        timeOption[j], r = rNow, b = bNow, sigmaOption)@price
                rNow <<- r
                bNow <<- b
            }
    
            persp(x, y, z,
                theta = theta, phi = phi,
                ticktype = "detailed",
                col = "steelblue", 
                shade = 0.5, 
                border = TRUE) -> Option  
                
            ZZ = GBSOption(TypeFlag, S, X, Time, r=rNow, b=bNow, sigma)@price     
            XX <<- sigma^2*Time
            YY <<- S/X 
               
            points(trans3d(XX, YY, ZZ, pm = Option), pch = 19, col = "orange")
            title(main = paste(
                TypeFlagText[TypeFlag], as.character(signif(ZZ, 5))))
            mS = signif(S, 3)
            mX = signif(X, 3)
            mSigma = round(sigma, digits = 2)
            mTime = round(Time, digits = 2)
            mText = paste(
                "S =", mS,
                "| X =", mX,
                "| Time =", mTime,
                "| sigma =", mSigma)
            mtext(mText)
        }
        
        # Initialization:
        TypeFlag <<- TypeFlag
        rNow <<- r
        bNow <<- b
        N <<- N
        
        Smin = S*(1-span)
        Smax = S*(1+span)
        Sres = (Smax-Smin)/N
        Son  = (Smin+Smax)/2
        
        Xmin = X*(1-span)
        Xmax = X*(1+span)
        Xres = (Xmax-Xmin)/N
        Xon  = (Xmin+Xmax)/2
        
        sOption <<- seq(Smin, Smax, by = Sres)
        xOption <<- Xon
        nX <<- length(sOption)
        timeOption <<- seq(1e-6, 3, length = N)
        sigmaOption <<- 0.25
        nY <<- length(timeOption)
        z <<- matrix(rep(0, nX*nY), ncol = nX) 
        for (j in 1:nY)
            z[j, ] <<- GBSOption(TypeFlag, sOption, xOption, timeOption[j], 
            r = rNow, b = bNow, sigmaOption)@price
        x <<- sigmaOption^2*timeOption
        y <<- sOption/xOption 
        
        # Open Slider Menu:
        plot.names = c("Plot - theta", "... phi")
        .sliderMenu(refresh.code,
           names =       c( "S",  "X", "Time", "sigma",  "r",  "b", plot.names),
           minima =      c(Smin, Xmin,   1e-6,   0.005, 0.01, 0.01,  -180,   0),
           maxima =      c(Smax, Xmax,   3.00,   0.500, 0.20, 0.20,   180, 360),
           resolutions = c(Sres, Xres,   0.10,   0.005, 0.01, 0.01,     2,   2),
           starts =      c( Son,  Xon,   1.00,   0.250, 0.10, 0.10,   -40,  30))
    }
    
    
    # Try:
    # .GBSOptionSlider("p")
 
    # Return Value:
    return()    
}

 
# ------------------------------------------------------------------------------
    

test.GBSGreeksSlider = 
function()
{         
    .GBSGreeksSlider = 
    function(TypeFlag = "c", S = 100, X = 100, Time = 1, 
    r = 0.10, b = 0.10, sigma = 0.25, span = 0.25, N = 40) 
    {
        # Internal Function:
        refresh.code = function(...)
        {
            # Sliders:
            S     = .sliderMenu(no = 1)
            X     = .sliderMenu(no = 2)
            Time  = .sliderMenu(no = 3)
            sigma = .sliderMenu(no = 4) 
            r     = .sliderMenu(no = 5) 
            b     = .sliderMenu(no = 6) 
            theta = .sliderMenu(no = 7)
            phi   = .sliderMenu(no = 8)
            
            Selection = "Gamma"
            TypeFlagText = c(c = "Call:", p = "Put:")
        
            if (r != rNow | b != bNow) {
                for (j in 1:nY)
                    z[j, ] <<- GBSGreeks(Selection, TypeFlag, sOption, xOption, 
                        timeOption[j], r = rNow, b = bNow, sigmaOption) 
                rNow <<- r
                bNow <<- b
            }
    
            persp(x, y, z,
                theta = theta, phi = phi,
                ticktype = "detailed",
                col = "steelblue", 
                shade = 0.5, 
                border = TRUE) -> Option  
                
            ZZ = GBSGreeks(Selection, TypeFlag, S, X, Time, 
                r = rNow, b = bNow, sigma) 
            XX <<- sigma^2*Time
            YY <<- S/X 
               
            points(trans3d(XX, YY, ZZ, pm = Option), pch = 19, col = "orange")
            title(main = paste(
                TypeFlagText[TypeFlag], as.character(signif(ZZ, 5))))
            mS = signif(S, 3)
            mX = signif(X, 3)
            mSigma = round(sigma, digits = 2)
            mTime = round(Time, digits = 2)
            mText = paste(
                "S =", mS,
                "| X =", mX,
                "| Time =", mTime,
                "| sigma =", mSigma)
            mtext(mText)
        }
        
        # Initialization:
        TypeFlag <<- TypeFlag
        rNow <<- r
        bNow <<- b
        N <<- N
        
        Smin = S*(1-span)
        Smax = S*(1+span)
        Sres = (Smax-Smin)/N
        Son  = (Smin+Smax)/2
        
        Xmin = X*(1-span)
        Xmax = X*(1+span)
        Xres = (Xmax-Xmin)/N
        Xon  = (Xmin+Xmax)/2
        
        sOption <<- seq(Smin, Smax, by = Sres)
        xOption <<- Xon
        nX <<- length(sOption)
        timeOption <<- seq(0, 3, length = N+1)[-1]
        sigmaOption <<- 0.25
        nY <<- length(timeOption)
        z <<- matrix(rep(0, nX*nY), ncol = nX) 
        for (j in 1:nY)
            z[j, ] <<- GBSGreeks("Gamma", TypeFlag, sOption, xOption, 
                timeOption[j], r = rNow, b = bNow, sigmaOption) 
        x <<- sigmaOption^2*timeOption
        y <<- sOption/xOption 
        
        # Open Slider Menu:
        plot.names = c("Plot - theta", "... phi")
        .sliderMenu(refresh.code,
           names =       c( "S",  "X", "Time", "sigma",  "r",  "b", plot.names),
           minima =      c(Smin, Xmin,   1e-6,   0.005, 0.01, 0.01,  -180,   0),
           maxima =      c(Smax, Xmax,   3.00,   0.500, 0.20, 0.20,   180, 360),
           resolutions = c(Sres, Xres,   0.10,   0.005, 0.01, 0.01,     2,   2),
           starts =      c( Son,  Xon,   1.00,   0.250, 0.10, 0.10,   -40,  30))
    }
        
    # Try
    # .GBSGreeksSlider("c")

    
    # Return Value:
    return()    
}


################################################################################

