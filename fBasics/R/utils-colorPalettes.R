
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
# You should have received A copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 COLOR PALETTES:
#  rainbowPalette            Creates a rainbow color palette          
#  heatPalette               Creates a heat color palette
#  terrainPalette            Creates a terrain color palette
#  topoPalette               Creates a topo color palette 
#  cmPalette                 Creates a cm color palette
#  greyPalette               Creates a grey palette
#  timPalette                Creates a cyan, yellow, to orange palette
# FUNCTION:                 COLOR RAMPS:
#  rampPalette               Creates a color ramp palette
#  seqPalette                Creates a sequential color palette
#  divPalette                Creates a diverging color palette
#  qualiPalette              Creates a qualitative color palette 
#  focusPalette              Creates a focus color palette
#  monoPalette               Creates a mono color palette
################################################################################


################################################################################
# FUNCTION:                 DESCRIPTION:
#  rainbowPalette            Creates a rainbow color palette          
#  heatPalette               Creates a heat color palette
#  terrainPalette            Creates a terrain color palette
#  topoPalette               Creates a topo color palette 
#  cmPalette                 Creates a cm color palette


rainbowPalette <- 
function(n = 64, ...) 
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Return Value:
    rainbow(n = n, ...)
}


# ------------------------------------------------------------------------------


heatPalette <- 
function(n = 64, ...) 
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Return Value:
    heat.colors(n, ...)
}


# ------------------------------------------------------------------------------


terrainPalette <- 
function(n = 64, ...) 
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Return Value:
    terrain.colors(n, ...)
}


# ------------------------------------------------------------------------------


topoPalette <- 
function(n = 64, ...) 
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Return Value:
    topo.colors(n, ...)
}


# ------------------------------------------------------------------------------


cmPalette <- 
function(n = 64, ...) 
{
    # A function implemented by Diethelm Wuertz
    
    # FUNCTION:
    
    # Return Value:
    cm.colors(n, ...)
}


################################################################################
# FUNCTION:                 COLOR/GREY PALETTES:
#  greyPalette               Creates a grey palette
#  timPalette                Creates a cyan, yellow, to orange palette


greyPalette <- 
function(n = 64, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Create a vector of n gamma-corrected gray colors. 

    # Arguments:
    #   n - the number of greys to be constructed
    #   start, end - the range of the color palette
    #   gamma a gamma-correction 
    
    # Value:
    #   returns a grey palette like rainbow does
    #   for color palettes

    # FUNCTION:

    # Compose:
    ans = gray.colors(n, ...)


    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


timPalette <- 
function(n = 64)
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a cyan, yellow, to orange palette
    
    # Notes:
    #   'Tim.colors' in 'fields' package goes from blue to red, and passes
    #   through the colors cyan, yellow, and orange. Also known as Jet
    #   color-map in Matlab. You can also easily design your own color map
    #   using 'rgb' function from 'gdDevices'.
    #   From:  <Jaroslaw.W.Tuszynski@saic.com>

    # FUNCTION:
    
    orig = c(
        "#00008F", "#00009F", "#0000AF", "#0000BF", "#0000CF",
        "#0000DF", "#0000EF", "#0000FF", "#0010FF", "#0020FF",
        "#0030FF", "#0040FF", "#0050FF", "#0060FF", "#0070FF",
        "#0080FF", "#008FFF", "#009FFF", "#00AFFF", "#00BFFF",
        "#00CFFF", "#00DFFF", "#00EFFF", "#00FFFF", "#10FFEF",
        "#20FFDF", "#30FFCF", "#40FFBF", "#50FFAF", "#60FF9F",
        "#70FF8F", "#80FF80", "#8FFF70", "#9FFF60", "#AFFF50",
        "#BFFF40", "#CFFF30", "#DFFF20", "#EFFF10", "#FFFF00",
        "#FFEF00", "#FFDF00", "#FFCF00", "#FFBF00", "#FFAF00",
        "#FF9F00", "#FF8F00", "#FF8000", "#FF7000", "#FF6000",
        "#FF5000", "#FF4000", "#FF3000", "#FF2000", "#FF1000",
        "#FF0000", "#EF0000", "#DF0000", "#CF0000", "#BF0000",
        "#AF0000", "#9F0000", "#8F0000", "#800000")
    if (n == 64) return(orig)
    rgb.tim = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , 64)
    xg = seq(0, 1, , n)
    for (k in 1:3) { 
        hold = spline(x, rgb.tim[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[, k] = round(hold)
    }
    ans = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    
    # Return Value:
    ans
}


################################################################################
# Package: colorRamps
# Type: Package
# Title: Builds pleasing color tables
# Version: 1.0
# Date: 2007-04-05
# Author: Tim Keitt
# Maintainer: Tim Keitt <tkeitt@gmail.com>
# Description: Builds single and double gradient color maps
# License: GPL
# Packaged: Thu Apr  5 16:34:42 2007; tkeitt


.blue2redPalette <-  
function(n)
{
    # A copy from contributed R-package colorRamps
    
    # FUNCTION:
    
    # Color Ramp:
    n2 = ceiling(n / 2)
    red = rep(c(0, 1), each = n2)[1:n]
    green = 1 - abs(seq(-1, 1, length.out = n))
    blue = rev(red)
    
    # Return Value:
    rgb(red, green, blue)
}


# ------------------------------------------------------------------------------


.green2redPalette <-  
function(n)
{
    # A copy from contributed R-package colorRamps
    
    # FUNCTION:
    
    # Color Ramp:
    n2 = ceiling(n / 2)
    red = rep(c(0, 1), each = n2)[1:n]
    blue = 1 - abs(seq(-1, 1, length.out = n))
    green = rev(red)
    
    # Return Value:
    rgb(red, green, blue)
  }


# ------------------------------------------------------------------------------


.blue2greenPalette <-  
function(n)
{
    # A copy from contributed R-package colorRamps
    
    # FUNCTION:
    
    # Color Ramp:
    n2 = ceiling(n / 2)
    green = rep(c(0, 1), each = n2)[1:n]
    red = 1 - abs(seq(-1, 1, length.out = n))
    blue = rev(green)
    
    # Return Value:
    rgb(red, green, blue)
}


# ------------------------------------------------------------------------------


.purple2greenPalette <-  
function(n)
{
    # A copy from contributed R-package colorRamps
    
    # FUNCTION:
    
    # Color Ramp:
    red = rep(0.5, length.out = n)
    green = seq(0, 1, length.out = n)
    blue = rev(green)
    
    # Return Value:
    rgb(red, green, blue)
}


# ------------------------------------------------------------------------------


.blue2yellowPalette <-  
function(n)
{
    # A copy from contributed R-package colorRamps
    
    # FUNCTION:
    
    # Color Ramp:
    red = seq(0, 1, length.out = n)
    green = red
    blue = rev(red)
    
    # Return Value:
    rgb(red, green, blue)
}


# ------------------------------------------------------------------------------


.cyan2magentaPalette <-  
function(n)
{
    # A copy from contributed R-package colorRamps
    
    # FUNCTION:
    
    # Color Ramp:
    red = seq(0, 1, length.out = n)
    green = rev(red)
    blue = rep(1, n)
    
    # Return Value:
    rgb(red, green, blue)
}


# ------------------------------------------------------------------------------


rampPalette <-
function(n, name = c("blue2red", "green2red", "blue2green",     
    "purple2green", "blue2yellow", "cyan2magenta"))
{
    # Description:
    #   Creates a color ramp palette
    
    # FUNCTION:
    
    # Color Ramp:
    name = match.arg(name)
    funPalette = match.fun(paste(".", name, "Palette", sep = ""))
    ans = funPalette(n)
    
    # Return Value:
    ans
}  


################################################################################
# Package: RColorBrewer
# Version: 1.0-2
# Date: 2007-10-21
# Title: ColorBrewer palettes
# Author: Erich Neuwirth <erich.neuwirth@univie.ac.at>
# Maintainer: Erich Neuwirth <erich.neuwirth@univie.ac.at>
# Depends: R (>= 2.0.0)
# Description: The packages provides palettes for drawing nice maps
#   shaded according to a variable.
# License: Apache License 2.0

            
seqPalette <-
function(n, name = c("Blues", "BuGn", "BuPu", "GnBu", "Greens", 
    "Greys", "Oranges", "OrRd", "PuBu", "PuBuGn", "PuRd", "Purples", 
    "RdPu", "Reds", "YlGn", "YlGnBu", "YlOrBr", "YlOrRd"))
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a sequential color palette
    
    # FUNCTION:
    
    # Color Sets:
    #   Blues,  BuGn, BuPu,    GnBu, Greens, Greys, Oranges, OrRd,   PuBu,
    #   PuBuGn, PuRd, Purples, RdPu, Reds,   YlGn,  YlGnBu,  YlOrBr, YlOrRd.
    Blues = 
        rgb(c(247,222,198,158,107,66,33,8,8),
            c(251,235,219,202,174,146,113,81,48),
            c(255,247,239,225,214,198,181,156,107),maxColorValue=255)    
    BuGn =  
        rgb(c(247,229,204,153,102,65,35,0,0),
            c(252,245,236,216,194,174,139,109,68),
            c(253,249,230,201,164,118,69,44,27),maxColorValue=255)        
    BuPu =  
        rgb(c(247,224,191,158,140,140,136,129,77),
            c(252,236,211,188,150,107,65,15,0),
            c(253,244,230,218,198,177,157,124,75),maxColorValue=255)     
    GnBu = 
        rgb(c(247,224,204,168,123,78,43,8,8),
            c(252,243,235,221,204,179,140,104,64),
            c(240,219,197,181,196,211,190,172,129),maxColorValue=255)        
    Greens = 
        rgb(c(247,229,199,161,116,65,35,0,0),
            c(252,245,233,217,196,171,139,109,68),
            c(245,224,192,155,118,93,69,44,27),maxColorValue=255)   
    Greys =  
        rgb(c(255,240,217,189,150,115,82,37,0),
            c(255,240,217,189,150,115,82,37,0),
            c(255,240,217,189,150,115,82,37,0),maxColorValue=255)         
    Oranges =  
        rgb(c(255,254,253,253,253,241,217,166,127),
            c(245,230,208,174,141,105,72,54,39),
            c(235,206,162,107,60,19,1,3,4),maxColorValue=255)         
    OrRd =  
        rgb(c(255,254,253,253,252,239,215,179,127),
            c(247,232,212,187,141,101,48,0,0),
            c(236,200,158,132,89,72,31,0,0),maxColorValue=255)        
    PuBu = 
        rgb(c(255,236,208,166,116,54,5,4,2),
            c(247,231,209,189,169,144,112,90,56),
            c(251,242,230,219,207,192,176,141,88),maxColorValue=255)          
    PuBuGn = 
        rgb(c(255,236,208,166,103,54,2,1,1),
            c(247,226,209,189,169,144,129,108,70),
            c(251,240,230,219,207,192,138,89,54),maxColorValue=255)       
    PuOr =  
        rgb(c(127,179,224,253,254,247,216,178,128,84,45),
            c(59,88,130,184,224,247,218,171,115,39,0),
            c(8,6,20,99,182,247,235,210,172,136,75),maxColorValue=255)       
    PuRd = 
        rgb(c(247,231,212,201,223,231,206,152,103),
            c(244,225,185,148,101,41,18,0,0),
            c(249,239,218,199,176,138,86,67,31),maxColorValue=255)         
    Purples =  
        rgb(c(252,239,218,188,158,128,106,84,63),
            c(251,237,218,189,154,125,81,39,0),
            c(253,245,235,220,200,186,163,143,125),maxColorValue=255)    
    RdPu = 
        rgb(c(255,253,252,250,247,221,174,122,73),
            c(247,224,197,159,104,52,1,1,0),
            c(243,221,192,181,161,151,126,119,106),maxColorValue=255)          
    Reds =  
        rgb(c(255,254,252,252,251,239,203,165,103),
            c(245,224,187,146,106,59,24,15,0),
            c(240,210,161,114,74,44,29,21,13),maxColorValue=255)      
    YlGn = 
        rgb(c(255,247,217,173,120,65,35,0,0),
            c(255,252,240,221,198,171,132,104,69),
            c(229,185,163,142,121,93,67,55,41),maxColorValue=255)         
    YlGnBu = 
        rgb(c(255,237,199,127,65,29,34,37,8),
            c(255,248,233,205,182,145,94,52,29),
            c(217,177,180,187,196,192,168,148,88),maxColorValue=255)         
    YlOrBr = 
        rgb(c(255,255,254,254,254,236,204,153,102),
            c(255,247,227,196,153,112,76,52,37),
            c(229,188,145,79,41,20,2,4,6),maxColorValue=255)         
    YlOrRd =
        rgb(c(255,255,254,254,253,252,227,189,128),
            c(255,237,217,178,141,78,26,0,0),
            c(204,160,118,76,60,42,28,38,38),maxColorValue=255)        
                                  
    # Compose:
    name = match.arg(name)
    orig =  eval(parse(text = name))
    rgb = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , length(orig))
    xg = seq(0, 1, , n)
    for (k in 1:3) { 
        hold = spline(x, rgb[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[, k] = round(hold)
    }
    palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    
    # Return Value:
    palette     
}


# ------------------------------------------------------------------------------
            

divPalette <-
function(n, name = c("BrBG", "PiYG", "PRGn", "PuOr", "RdBu", "RdGy", 
    "RdYlBu", "RdYlGn", "Spectral"))
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a diverging color palette
    
    # FUNCTION:
    
    # Color Sets:
    # BrBG, PiYG, PRGn, PuOr, RdBu, RdGy, RdYlBu, RdYlGn, Spectral.
    BrBG = 
        rgb(c(84,140,191,223,246,245,199,128,53,1,0),
            c(48,81,129,194,232,245,234,205,151,102,60),
            c(5,10,45,125,195,245,229,193,143,94,48),maxColorValue=255) 
    PiYG = 
        rgb(c(142,197,222,241,253,247,230,184,127,77,39),
            c(1,27,119,182,224,247,245,225,188,146,100),
            c(82,125,174,218,239,247,208,134,65,33,25),maxColorValue=255)         
    PRGn =  
        rgb(c(64,118,153,194,231,247,217,166,90,27,0),
            c(0,42,112,165,212,247,240,219,174,120,68),
            c(75,131,171,207,232,247,211,160,97,55,27),maxColorValue=255)  
    PuOr =  
        rgb(c(127,179,224,253,254,247,216,178,128,84,45),
            c(59,88,130,184,224,247,218,171,115,39,0),
            c(8,6,20,99,182,247,235,210,172,136,75),maxColorValue=255)       
    RdBu = 
        rgb(c(103,178,214,244,253,247,209,146,67,33,5),
            c(0,24,96,165,219,247,229,197,147,102,48),
            c(31,43,77,130,199,247,240,222,195,172,97),maxColorValue=255) 
    RdGy = 
        rgb(c(103,178,214,244,253,255,224,186,135,77,26),
            c(0,24,96,165,219,255,224,186,135,77,26),
            c(31,43,77,130,199,255,224,186,135,77,26),maxColorValue=255)           
    RdYlBu = 
        rgb(c(165,215,244,253,254,255,224,171,116,69,49),
            c(0,48,109,174,224,255,243,217,173,117,54),
            c(38,39,67,97,144,191,248,233,209,180,149),maxColorValue=255)         
    RdYlGn =
        rgb(c(165,215,244,253,254,255,217,166,102,26,0),
            c(0,48,109,174,224,255,239,217,189,152,104),
            c(38,39,67,97,139,191,139,106,99,80,55),maxColorValue=255)           
    Spectral = 
        rgb(c(158,213,244,253,254,255,230,171,102,50,94),
            c(1,62,109,174,224,255,245,221,194,136,79),
            c(66,79,67,97,139,191,152,164,165,189,162),maxColorValue=255) 
                   
    # Compose:
    name = match.arg(name)
    orig =  eval(parse(text = name))
    rgb = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , length(orig))
    xg = seq(0, 1, , n)
    for (k in 1:3) { 
        hold = spline(x, rgb[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[, k] = round(hold)
    }
    palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    
    # Return Value:
    palette   
}


# ------------------------------------------------------------------------------


qualiPalette <-
function(n, name = c("Accent", "Dark2", "Paired", "Pastel1", 
    "Pastel2", "Set1", "Set2", "Set3"))
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a qualitative color palette
    
    # FUNCTION:
    
    # Color Sets:      
    Accent = 
        rgb(c(127,190,253,255,56,240,191,102),
            c(201,174,192,255,108,2,91,102),
            c(127,212,134,153,176,127,23,102),maxColorValue=255)                       
    Dark2 = 
        rgb(c(27,217,117,231,102,230,166,102),
            c(158,95,112,41,166,171,118,102),
            c(119,2,179,138,30,2,29,102),maxColorValue=255)
    Paired = 
        rgb(c(166,31,178,51,251,227,253,255,202,106,255,177),
            c(206,120,223,160,154,26,191,127,178,61,255,89),
            c(227,180,138,44,153,28,111,0,214,154,153,40),maxColorValue=255)        
    Pastel1 =
        rgb(c(251,179,204,222,254,255,229,253,242),
            c(180,205,235,203,217,255,216,218,242),
            c(174,227,197,228,166,204,189,236,242),maxColorValue=255)         
    Pastel2 =  
        rgb(c(179,253,203,244,230,255,241,204),
            c(226,205,213,202,245,242,226,204),
            c(205,172,232,228,201,174,204,204),maxColorValue=255)                         
    Set1 = 
        rgb(c(228,55,77,152,255,255,166,247,153),
            c(26,126,175,78,127,255,86,129,153),
            c(28,184,74,163,0,51,40,191,153),maxColorValue=255)        
    Set2 = 
        rgb(c(102,252,141,231,166,255,229,179),
            c(194,141,160,138,216,217,196,179),
            c(165,98,203,195,84,47,148,179),maxColorValue=255)          
    Set3 = 
        rgb(c(141,255,190,251,128,253,179,252,217,188,204,255),
            c(211,255,186,128,177,180,222,205,217,128,235,237),
            c(199,179,218,114,211,98,105,229,217,189,197,111),maxColorValue=255)          
             
    # Compose:
    name = match.arg(name)
    orig =  eval(parse(text = name))
    rgb = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , length(orig))
    xg = seq(0, 1, , n)
    for (k in 1:3) { 
        hold = spline(x, rgb[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[, k] = round(hold)
    }
    palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    
    # Return Value:
    palette   
}


################################################################################
# Package: PerformanceAnalytics
# Type: Package
# Title: Econometric tools for performance and risk analysis.
# Version: 0.9.5
# Date: 2007-06-29
# Author: Peter Carl, Brian G. Peterson
# Maintainer: Brian G. Peterson <brian@braverock.com>
# Description: Library of econometric functions for performance and risk 
#   analysis. This library aims to aid practitioners and researchers in 
#   utilizing the latest research in analysis of non-normal return streams.  
#   In general, this library is most tested on return (rather than price) 
#   data on a monthly scale, but most functions will work with daily or 
#   irregular return data as well.
# Depends: R (>= 2.4.0), fExtremes, fPortfolio, quadprog, tseries, Hmisc
# License: GPL
# URL: http://braverock.com/R/
# Packaged: Tue Jul 10 04:30:47 2007

        
focusPalette <-
function(n, name = c("redfocus", "greenfocus", "bluefocus"))
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a color palette for graphs
    
    # Source:
    #   Contributed R package PerformanceAnalytics
    
    # Details:
    #   This is not a function, per se, but a way to set up specific color 
    #   pallets for use in the charts we use. These pallets have been 
    #   designed to create readable, comparable line and bar graphs with 
    #   specific objectives outlined before each category below.
    #   We use this approach rather than generating them on the fly for two 
    #   reasons: 1) fewer dependencies on libraries that don't need to be 
    #   called dynamically; and 2) to guarantee the color used for the n-th 
    #   column of data.
    #
    #   FOCUS PALETTE
    #   Colorsets designed to provide focus to the data graphed as the first 
    #   element. This palette is best used when there is clearly an important 
    #   data set for the viewer to focus on, with the remaining data being 
    #   secondary, tertiary, etc. Later elements graphed in diminishing 
    #   values of gray. These were generated with RColorBrewer, using the 8 
    #   level "grays" palette and replacing the darkest with the focus color.
    #   For best results, replace the highlight color with the first color 
    #   of the equal weighted palette from below. This will coordinate 
    #   charts with different purposes.
    # 
    #   MONOCHROME PALETTES: 
    #   Colorsets for monochrome color displays.
  
    # FUNCTION:
    
    # Match Arguments:
    name = match.arg(name)
    
    # Focus Palettes:
    redfocus = c(
        "#CB181D", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", 
        "#D9D9D9", "#F0F0F0")
    greenfocus = c(
        "#41AB5D", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", 
        "#D9D9D9", "#F0F0F0")
    bluefocus = c(
        "#0033FF", "#252525", "#525252", "#737373", "#969696", "#BDBDBD", 
        "#D9D9D9", "#F0F0F0")

    # Compose:
    orig =  eval(parse(text = name))
    rgb = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , length(orig))
    xg = seq(0, 1, , n)
    for (k in 1:3) { 
        hold = spline(x, rgb[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[, k] = round(hold)
    }
    palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    
    # Return Value:
    palette   
}


# ------------------------------------------------------------------------------


monoPalette <-
function(n, name = c("redmono", "greenmono", "bluemono"))
{
    # A function implemented by Diethelm Wuertz
    
    # Description:
    #   Creates a mono color palette
    
    # Source:
    #   Contributed R package PerformanceAnalytics
  
    # FUNCTION:
    
    # Match Arguments:
    name = match.arg(name)
          
    # Monochrome Palettes:
    redmono = c(
        "#99000D", "#CB181D", "#EF3B2C", "#FB6A4A", "#FC9272", "#FCBBA1", 
        "#FEE0D2", "#FFF5F0")
    greenmono = c(
        "#005A32", "#238B45", "#41AB5D", "#74C476", "#A1D99B", "#C7E9C0", 
        "#E5F5E0", "#F7FCF5")
    bluemono = c(
        "#084594", "#2171B5", "#4292C6", "#6BAED6", "#9ECAE1", "#C6DBEF", 
        "#DEEBF7", "#F7FBFF")

    # Compose:
    orig =  eval(parse(text = name))
    rgb = t(col2rgb(orig))
    temp = matrix(NA, ncol = 3, nrow = n)
    x = seq(0, 1, , length(orig))
    xg = seq(0, 1, , n)
    for (k in 1:3) { 
        hold = spline(x, rgb[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[, k] = round(hold)
    }
    palette = rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    
    # Return Value:
    palette   
}


################################################################################

    