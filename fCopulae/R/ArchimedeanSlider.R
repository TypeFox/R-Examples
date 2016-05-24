
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


################################################################################
# FUNCTION:                  ARCHIMEDEAN COPULAE SLIDERS:
#  rarchmSlider               Displays interactively Archimedean probability
#  parchmSlider               Displays interactively Archimedean probability 
#  .parchmPerspSlider          Perspective Archimedean probability slider
#  .parchmContourSlider        Contour Archimedean probability slider
#  darchmSlider                Displays interactively archimedean density 
#  .darchmPerspSlider          Perspective Archimedean density slider
#  .darchmContourSlider        Contour Archimedean density slider
################################################################################


################################################################################
rarchmSlider <- 
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Displays interactively perspective plots of probability
    
    # FUNCTION:
    
    # Graphic Frame:
    par(mfrow = c(1, 1))
    
    # Internal Function:
    refresh.code <- function(...)
    {
        # Sliders:
        #           1       5        10        15        20
        Counter = c(1,2,3,2,4,2,5,2,5,5,6,2,7,2,2,7,4,8,7,7,2,5)
        Copula = as.integer(.sliderMenu(no = 1))
        No = Counter[Copula]
        N = .sliderMenu(no = 2)
        alpha = .sliderMenu(no = No+2)
                     
        # There is no known Copula for the following bounds:
        eps = 1.0e-6
        if (Copula == 11) if (alpha == 0.5) alpha = 0.5 - eps
        if (Copula == 13) if (alpha == 0.0) alpha = eps
        
        # Title:
        Names = c(
            "- Clayton", "", 
            "- Ali-Mikhail-Hag", 
            "- Gumbel-Hougard", 
            "- Frank",
            "- Joe-Frank", "", "", 
            "- Gumbel-Barnett", "", "", "", "", "", 
            "- Genest-Ghoudi", "", "", "", "", "", "", "")      
        Title = paste("Archimedean Copula No:", as.character(Copula), 
            Names[Copula], "\nalpha = ", as.character(alpha)) 
        
        # Plot: 
        R = rarchmCopula(n = N, alpha = alpha, type = as.character(Copula))
        plot(R, xlab = "U", ylab = "V", pch = 19, col = "steelblue")
        grid()
        title(main = Title)
                           
        # Reset Frame:
        par(mfrow = c(1, 1))
    }
  
    # Open Slider Menu:
    C2 = "2-4-6-8-12-14-15-21"
    C = c("1", C2, "3", "5-17", "7-9-10-22", "11", "13-16-19-20","18")                                           
    L = c( -1,  1,  -1,     -B,           0,    0,            0,   2 )
    U = c(  B,  B,   1,      B,           1,  0.5,            B,   B )
    A = c(0.5,  2, 0.5,      1,         0.5,  0.2,            1,   3 ) 
    V = rep(0.1, 8)
    .sliderMenu(refresh.code,
        names       = c("Copula",  "N", C),
        minima      = c(       1,  100, L),
        maxima      = c(      22, 1000, U),
        resolutions = c(       1,  100, V),
        starts      = c(       1,  100, A)) 
}


################################################################################
parchmSlider <- 
    function(type = c("persp", "contour"), B = 10)
{   
    # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Displays interactively plots of probability
    
    # Description:
    #   Displays interactively plots of probability
    
    # Arguments:
    #   type - a character string specifying the plot type.
    #       Either a perspective plot which is the default or
    #       a contour plot with an underlying image plot will
    #       be created.
    #   B - the maximum slider menu value when the boundary
    #       value is infinite. By default this is set to 10.
    
    # FUNCTION:
    
    # Match Arguments:
    type = match.arg(type)
    
    # Plot:
    if (type[1] == "persp")
        .parchmPerspSlider(B = B)
    if (type[1] == "contour")
        .parchmContourSlider(B = B)
        
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
.parchmPerspSlider =
function(B = 5)
{   # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Displays interactively perspective plots of probability
    
    # FUNCTION:
    
    # Graphic Frame:
    par(mfrow = c(1, 1))
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        Counter = c(1,2,3,2,4,2,5,2,5,5,6,2,7,2,2,7,4,8,7,7,2,5)
        Copula = as.integer(.sliderMenu(no = 1))
        No = Counter[Copula]
        N = .sliderMenu(no = 2)
        alpha = .sliderMenu(no = No+2)
        theta = .sliderMenu(no = 11)
        phi = .sliderMenu(no = 12)
          
        # Skip:
        if (Copula == 11) if (alpha == 0.5) return(invisible())
        if (Copula == 13) if (alpha == 0)  return(invisible())
        
        # Do we have a strict Copula?
        strict = c(
            "Yes","No","Yes","Yes","Yes","Yes","No","No","Yes","Yes",
            "No","Yes","Yes","Yes","No","Yes","Yes","No","Yes","Yes", 
            "No","Yes")[Copula]
        if (alpha < 0 & Copula == 1) strict[1] = "No"
        if (alpha == 0 & Copula == 16) strict[16] = "No"
        
        # What is the Range?
        RANGE = c(
            "[-1|Inf)", "[1|Inf)", "[-1|1)", "(-Inf|Inf)", "(0|1]", 
            "(0|0.5]", "(0|Inf)", "[2|Inf)")[No]
                 
        # Which one is the Limit Copula?
        limitTitle = rep("NA", times = 22)
        if (alpha == -1) 
            limitTitle = c(
                "W ", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                "NA", "NA", "NA", "NA", "NA", "NA", "Pi", "NA", "NA", "NA",
                "NA", "NA")   
        if (alpha == 0) 
            limitTitle = c(
                "Pi", "NA", "Pi", "NA", "Pi", "NA", "W ", "NA", "Pi", "Pi",
                "Pi", "NA", "NA", "NA", "NA", "W ", "NA", "NA", "L ", "Pi",
                "NA", "Pi")    
        if (alpha == 1) 
            limitTitle = c(
                "L ", "W ", "L ", "Pi", "NA", "Pi", "Pi", "W ", "NA", "NA",
                "NA", "L ", "Pi", "L ", "W ", "NA", "NA", "NA", "NA", "NA",
                "W ", "NA")
        limitTitle = limitTitle[Copula]
        if (limitTitle == "NA") {
            limitTitle = " "
        } else {
            limitTitle = paste("  Copula = ", limitTitle[1])
        }
        
        # Tau/Rho:
        Tau = round(approx(.ALPHA[, Copula], .TAU[, Copula], xout = alpha)$y,
            digits = 3)
        Rho = round(approx(.ALPHA[, Copula], .RHO[, Copula], xout = alpha)$y,
            digits = 3)
  
        # Title:
        Names = c(
            "- Clayton", "", "- Ali-Mikhail-Hag", "- Gumbel-Hougard", "- Frank",
            "- Joe-Frank", "", "", "- Gumbel-Barnett", "",
            "", "", "", "", "- Genest-Ghoudi", "", "", "", "", "", "", "")      
        Title = paste("Archimedean Copula No:", as.character(Copula), 
            Names[Copula], "\n", RANGE, " alpha =", as.character(alpha), 
            " tau =", as.character(Tau), " rho =", as.character(Rho)) 
        
        # Plot: 
        uv = grid2d(x = (0:N)/N)
        P = .parchm1Copula(u = uv, alpha = alpha, type = Copula, 
            output = "list")
        persp(P, theta = theta, phi = phi, col = "steelblue", shade = 0.5,
            ticktype = "detailed", cex = 0.5, xlab = "u", ylab = "v",
            zlab = "C(u,v)" )
        title(main = Title)
                           
        # Reset Frame:
        par(mfrow = c(1, 1))
    }
  
   
    # Open Slider Menu:
    B = 5
    C1 = "1: [-1,Inf]"
    C2 = "2-4-6-8-12-14-15-21: [1,Inf)"
    C3 = "3: [-1,1)"
    C4 = "5-17: (-Inf,Inf)|{0}"
    C5 = "7-9-10-22: (0,1]"
    C6 = "11: (0, 1/2]"
    C7 = "13-16-19-20: (0,Inf)"
    C8 = "18: [2, Inf)" 
    C = c(   C1, C2,   C3,   C4,   C5,   C6,   C7,  C8 )  
    L = c(   -1,  1,   -1,   -B,    0,    0,    0,   2 )
    U = c(    B,  B,    1,    B,    1,  0.5,    B,   B )
    A = c(  0.5,  2,  0.5,    1,  0.5,  0.2,    1,   3 ) 
    V = rep(0.01, 8)
    plot.names = c("Plot - theta", "... phi")
    .sliderMenu(refresh.code,
        names       = c("Copula", "N", C, plot.names),
        minima      = c(       1,  10, L, -180,    0),
        maxima      = c(      22, 100, U,  180,  360),
        resolutions = c(       1,  10, V,    1,    1),
        starts      = c(       1,  10, A,  -40,   30)) 
}


# ------------------------------------------------------------------------------
.parchmContourSlider <- 
    function(B = 5)
{   
    # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Displays interactively perspective plots of probability
    
    #FUNCTION:
    
    # Graphic Frame:
    par(mfrow = c(1, 1))
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        #           1       5        10        15        20
        Counter = c(1,2,3,2,4,2,5,2,5,5,6,2,7,2,2,7,4,8,7,7,2,5)
        Copula = as.integer(.sliderMenu(no = 1))
        No = Counter[Copula]
        N = .sliderMenu(no = 2)
        alpha = .sliderMenu(no = No+2)
        n.lev = .sliderMenu(no = 11)
        n.col = .sliderMenu(no = 12)
          
        # Skip:
        if (Copula == 11) if (alpha == 0.5) return(invisible())
        if (Copula == 13) if (alpha == 0)  return(invisible())
        
        # Do we have a strict Copula?
        strict = c(
            "Yes","No","Yes","Yes","Yes","Yes","No","No","Yes","Yes",
            "No","Yes","Yes","Yes","No","Yes","Yes","No","Yes","Yes", 
            "No","Yes")[Copula]
        if (alpha < 0 & Copula == 1) strict[1] = "No"
        if (alpha == 0 & Copula == 16) strict[16] = "No"
        
        # What is the Range?
        RANGE = c(
            "[-1|Inf)", "[1|Inf)", "[-1|1)", "(-Inf|Inf)", "(0|1]", 
            "(0|0.5]", "(0|Inf)", "[2|Inf)")[No]
                 
        # Which one is the Limit Copula?
        limitTitle = rep("NA", times = 22)
        if (alpha == -1) 
            limitTitle = c(
                "W ", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                "NA", "NA", "NA", "NA", "NA", "NA", "Pi", "NA", "NA", "NA",
                "NA", "NA")   
        if (alpha == 0) 
            limitTitle = c(
                "Pi", "NA", "Pi", "NA", "Pi", "NA", "W ", "NA", "Pi", "Pi",
                "Pi", "NA", "NA", "NA", "NA", "W ", "NA", "NA", "L ", "Pi",
                "NA", "Pi")    
        if (alpha == 1) 
            limitTitle = c(
                "L ", "W ", "L ", "Pi", "NA", "Pi", "Pi", "W ", "NA", "NA",
                "NA", "L ", "Pi", "L ", "W ", "NA", "NA", "NA", "NA", "NA",
                "W ", "NA")
        limitTitle = limitTitle[Copula]
        if (limitTitle == "NA") {
            limitTitle = " "
        } else {
            limitTitle = paste("  Copula = ", limitTitle[1])
        }
        
        # Tau/Rho:
        Tau = round(approx(.ALPHA[, Copula], .TAU[, Copula], xout = alpha)$y,
            digits = 3)
        Rho = round(approx(.ALPHA[, Copula], .RHO[, Copula], xout = alpha)$y,
            digits = 3)
  
        # Title:
        Names = c(
            "- Clayton", "", "- Ali-Mikhail-Hag", "- Gumbel-Hougard", "- Frank",
            "- Joe-Frank", "", "", "- Gumbel-Barnett", "",
            "", "", "", "", "- Genest-Ghoudi", "", "", "", "", "", "", "")      
        Title = paste("Archimedean Copula No:", as.character(Copula), 
            Names[Copula], "\n", RANGE, " alpha =", as.character(alpha), 
            " tau =", as.character(Tau), " rho =", as.character(Rho)) 
        
        # Plot:   
        uv = grid2d(x = (0:N)/N)
        P = .parchm1Copula(u = uv, alpha = alpha, type = Copula, 
            output = "list")
        image(P, col = heat.colors(n.col) )
        contour(P, xlab = "u", ylab = "v", nlevels = n.lev, add = TRUE)
        title(main = Title)
                           
        # Reset Frame:
        par(mfrow = c(1, 1))
    }
  
    # Open Slider Menu:
    C1 = "1: [-1,Inf]"
    C2 = "2-4-6-8-12-14-15-21: [1,Inf)"
    C3 = "3: [-1,1)"
    C4 = "5-17: (-Inf,Inf)|{0}"
    C5 = "7-9-10-22: (0,1]"
    C6 = "11: (0, 1/2]"
    C7 = "13-16-19-20: (0,Inf)"
    C8 = "18: [2, Inf)" 
    C = c(   C1, C2,   C3,   C4,   C5,   C6,   C7,  C8 )  
    L = c(   -1,  1,   -1,   -B,    0,    0,    0,   2 )
    U = c(    B,  B,    1,    B,    1,  0.5,    B,   B )
    A = c(  0.5,  2,  0.5,    1,  0.5,  0.2,    1,   3 ) 
    V = rep(0.01, 8)
    plot.names = c("Plot - levels", "... colors")
    .sliderMenu(refresh.code,
        names       = c("Copula", "N", C, plot.names),
        minima      = c(       1,  10, L,    5,   12),
        maxima      = c(      20, 100, U,  100,  256),
        resolutions = c(       1,  10, V,    5,    1),
        starts      = c(       1,  10, A,   10,   12)) 
}


################################################################################
darchmSlider <- 
    function(type = c("persp", "contour"), B = 10)
{   
    # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Displays interactively plots of density
    
    # Arguments:
    #   type - a character string specifying the plot type.
    #       Either a perspective plot which is the default or
    #       a contour plot with an underlying image plot will
    #       be created.
    #   B - the maximum slider menu value when the boundary
    #       value is infinite. By default this is set to 10.
    
    # FUNCTION:
    
    # Match Arguments:
    type = match.arg(type)
    
    # Plot:
    if (type == "persp")
        .darchmPerspSlider(B = B)
    if (type == "contour")
        .darchmContourSlider(B = B)
        
    # Return Value:
    invisible()
}


# ------------------------------------------------------------------------------
.darchmPerspSlider <- 
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Displays interactively perspective plots of density
    
    # FUNCTION:
    
    # Graphic Frame:
    par(mfrow = c(1, 1))
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        Counter = c(1,2,3,2,4,2,5,2,5,5,6,2,7,2,2,7,4,8,7,7,2,5)
        Copula = as.integer(.sliderMenu(no = 1))
        No = Counter[Copula]
        N = .sliderMenu(no = 2)
        alpha = .sliderMenu(no = No+2)
        theta = .sliderMenu(no = 11)
        phi = .sliderMenu(no = 12)
          
        # Skip:
        if (Copula == 11) if (alpha == 0.5) return(invisible())
        if (Copula == 13) if (alpha == 0)  return(invisible())
        
        # Do we have a strict Copula?
        strict = c(
            "Yes","No","Yes","Yes","Yes","Yes","No","No","Yes","Yes",
            "No","Yes","Yes","Yes","No","Yes","Yes","No","Yes","Yes", 
            "No","Yes")[Copula]
        if (alpha < 0 & Copula == 1) strict[1] = "No"
        if (alpha == 0 & Copula == 16) strict[16] = "No"
        
        # What is the Range?
        RANGE = c(
            "[-1|Inf)", "[1|Inf)", "[-1|1)", "(-Inf|Inf)", "(0|1]", 
            "(0|0.5]", "(0|Inf)", "[2|Inf)")[No]
                 
        # Which one is the Limit Copula?
        limitTitle = rep("NA", times = 22)
        if (alpha == -1) 
            limitTitle = c(
                "W ", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                "NA", "NA", "NA", "NA", "NA", "NA", "Pi", "NA", "NA", "NA",
                "NA", "NA")   
        if (alpha == 0) 
            limitTitle = c(
                "Pi", "NA", "Pi", "NA", "Pi", "NA", "W ", "NA", "Pi", "Pi",
                "Pi", "NA", "NA", "NA", "NA", "W ", "NA", "NA", "L ", "Pi",
                "NA", "Pi")    
        if (alpha == 1) 
            limitTitle = c(
                "L ", "W ", "L ", "Pi", "NA", "Pi", "Pi", "W ", "NA", "NA",
                "NA", "L ", "Pi", "L ", "W ", "NA", "NA", "NA", "NA", "NA",
                "W ", "NA")
        limitTitle = limitTitle[Copula]
        if (limitTitle == "NA") {
            limitTitle = " "
        } else {
            limitTitle = paste("  Copula = ", limitTitle[1])
        }
        
        # Tau/Rho:
        Tau = round(approx(.ALPHA[, Copula], .TAU[, Copula], xout = alpha)$y,
            digits = 3)
        Rho = round(approx(.ALPHA[, Copula], .RHO[, Copula], xout = alpha)$y,
            digits = 3)
  
        # Title:
        Names = c(
            "- Clayton", "", "- Ali-Mikhail-Hag", "- Gumbel-Hougard", "- Frank",
            "- Joe-Frank", "", "", "- Gumbel-Barnett", "",
            "", "", "", "", "- Genest-Ghoudi", "", "", "", "", "", "", "")      
        Title = paste("Archimedean Copula No:", as.character(Copula), 
            Names[Copula], "\n", RANGE, " alpha =", as.character(alpha), 
            " tau =", as.character(Tau), " rho =", as.character(Rho)) 
        
        # Plot: 
        uv = grid2d(x = (1:(N-1))/N)
        D = .darchm1Copula(u = uv, alpha = alpha, type = as.character(Copula), 
            output = "list")
        persp(D, theta = theta, phi = phi, col = "steelblue", shade = 0.5,
            ticktype = "detailed", cex = 0.5, xlab = "u", ylab = "v",
            zlab = "C(u,v)" )
        title(main = Title)
                           
        # Reset Frame:
        par(mfrow = c(1, 1))
    }
  
    # Open Slider Menu:
    B = 5
    C1 = "1: [-1,Inf]"
    C2 = "2-4-6-8-12-14-15-21: [1,Inf)"
    C3 = "3: [-1,1)"
    C4 = "5-17: (-Inf,Inf)|{0}"
    C5 = "7-9-10-22: (0,1]"
    C6 = "11: (0, 1/2]"
    C7 = "13-16-19-20: (0,Inf)"
    C8 = "18: [2, Inf)" 
    C = c(   C1, C2,   C3,   C4,   C5,   C6,   C7,  C8 )  
    L = c(   -1,  1,   -1,   -B,    0,    0,    0,   2 )
    U = c(    B,  B,    1,    B,    1,  0.5,    B,   B )
    A = c(  0.5,  2,  0.5,    1,  0.5,  0.2,    1,   3 ) 
    V = rep(0.1, 8)
    plot.names = c("Plot - theta", "... phi")
    .sliderMenu(refresh.code,
        names       = c("Copula", "N", C, plot.names),
        minima      = c(       1,  10, L, -180,    0),
        maxima      = c(      22, 100, U,  180,  360),
        resolutions = c(       1,  10, V,    1,    1),
        starts      = c(       1,  20, A,  -40,   30)) 
}


# ------------------------------------------------------------------------------
.darchmContourSlider <- 
    function(B = 10)
{   
    # A function implemented by Diethelm Wuertz
        
    # Description:
    #   Displays interactively perspective plots of density
    
    #FUNCTION:
    
    # Graphic Frame:
    par(mfrow = c(1, 1))
    
    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        Counter = c(1,2,3,2,4,2,5,2,5,5,6,2,7,2,2,7,4,8,7,7,2,5)
        Copula = as.integer(.sliderMenu(no = 1))
        No = Counter[Copula]
        N = .sliderMenu(no = 2)
        alpha = .sliderMenu(no = No+2)
        n.lev = .sliderMenu(no = 11)
        n.col = .sliderMenu(no = 12)
          
        # Skip:
        if (Copula == 11) if (alpha == 0.5) return(invisible())
        if (Copula == 13) if (alpha == 0)  return(invisible())
        
        # Do we have a strict Copula?
        strict = c(
            "Yes","No","Yes","Yes","Yes","Yes","No","No","Yes","Yes",
            "No","Yes","Yes","Yes","No","Yes","Yes","No","Yes","Yes", 
            "No","Yes")[Copula]
        if (alpha < 0 & Copula == 1) strict[1] = "No"
        if (alpha == 0 & Copula == 16) strict[16] = "No"
        
        # What is the Range?
        RANGE = c(
            "[-1|Inf)", "[1|Inf)", "[-1|1)", "(-Inf|Inf)", "(0|1]", 
            "(0|0.5]", "(0|Inf)", "[2|Inf)")[No]
                 
        # Which one is the Limit Copula?
        limitTitle = rep("NA", times = 22)
        if (alpha == -1) 
            limitTitle = c(
                "W ", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA", "NA",
                "NA", "NA", "NA", "NA", "NA", "NA", "Pi", "NA", "NA", "NA",
                "NA", "NA")   
        if (alpha == 0) 
            limitTitle = c(
                "Pi", "NA", "Pi", "NA", "Pi", "NA", "W ", "NA", "Pi", "Pi",
                "Pi", "NA", "NA", "NA", "NA", "W ", "NA", "NA", "L ", "Pi",
                "NA", "Pi")    
        if (alpha == 1) 
            limitTitle = c(
                "L ", "W ", "L ", "Pi", "NA", "Pi", "Pi", "W ", "NA", "NA",
                "NA", "L ", "Pi", "L ", "W ", "NA", "NA", "NA", "NA", "NA",
                "W ", "NA")
        limitTitle = limitTitle[Copula]
        if (limitTitle == "NA") {
            limitTitle = " "
        } else {
            limitTitle = paste("  Copula = ", limitTitle[1])
        }
        
        # Tau/Rho:
        Tau = round(approx(.ALPHA[, Copula], .TAU[, Copula], xout = alpha)$y,
            digits = 3)
        Rho = round(approx(.ALPHA[, Copula], .RHO[, Copula], xout = alpha)$y,
            digits = 3)
  
        # Title:
        Names = c(
            "- Clayton", "", "- Ali-Mikhail-Hag", "- Gumbel-Hougard", "- Frank",
            "- Joe-Frank", "", "", "- Gumbel-Barnett", "",
            "", "", "", "", "- Genest-Ghoudi", "", "", "", "", "", "", "")      
        Title = paste("Archimedean Copula No:", as.character(Copula), 
            Names[Copula], "\n", RANGE, " alpha =", as.character(alpha), 
            " tau =", as.character(Tau), " rho =", as.character(Rho)) 
        
        # Plot:   
        uv = grid2d(x = (1:(N-1)/N))
        D = .darchm1Copula(u = uv, alpha = alpha, type = as.character(Copula), 
            output = "list")
        image(D, xlim = c(0, 1), ylim = c(0,1), col = heat.colors(n.col) )
        contour(D, xlab = "u", ylab = "v", nlevels = n.lev, add = TRUE)
        title(main = Title)
                           
        # Reset Frame:
        par(mfrow = c(1, 1))
    }
  
    # Open Slider Menu:
    B = 5
    C1 = "1: [-1,Inf]"
    C2 = "2-4-6-8-12-14-15-21: [1,Inf)"
    C3 = "3: [-1,1)"
    C4 = "5-17: (-Inf,Inf)|{0}"
    C5 = "7-9-10-22: (0,1]"
    C6 = "11: (0, 1/2]"
    C7 = "13-16-19-20: (0,Inf)"
    C8 = "18: [2, Inf)" 
    C = c(   C1, C2,   C3,   C4,   C5,   C6,   C7,  C8 )  
    L = c(   -1,  1,   -1,   -B,    0,    0,    0,   2 )
    U = c(    B,  B,    1,    B,    1,  0.5,    B,   B )
    A = c(  0.5,  2,  0.5,    1,  0.5,  0.2,    1,   3 ) 
    V = rep(0.1, 8)
    plot.names = c("Plot - levels", "... colors")
    .sliderMenu(refresh.code,
        names       = c("Copula", "N", C, plot.names),
        minima      = c(       1,  10, L,   10,   12),
        maxima      = c(      22, 100, U,  100,  256),
        resolutions = c(       1,  10, V,   10,    1),
        starts      = c(       1,  30, A,   30,   64)) 
}

