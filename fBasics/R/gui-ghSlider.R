
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
# FUNCTION:             DESCRIPTION:
#  ghSlider              Displays generalized hyperbolic distribution function
################################################################################


ghSlider <-
function()
{
    # A function implemented by Diethelm Wuertz

    # Generalized Hyperbolic Distribution:
    #   dhyp(x, alpha = 1, beta = 0, delta = 1, lamba = 1)

    # FUNCTION:

    # Internal Function:
    refresh.code = function(...)
    {
        # Sliders:
        N      = .sliderMenu(no = 1)
        alpha  = .sliderMenu(no = 2)
        beta   = .sliderMenu(no = 3)
        delta  = .sliderMenu(no = 4)
        mu     = .sliderMenu(no = 5)
        lambda = .sliderMenu(no = 6)

        # Plot Data:
        xmin = round(qgh(0.01, alpha, beta, delta, mu, lambda), digits = 2)
        xmax = round(qgh(0.99, alpha, beta, delta, mu, lambda), digits = 2)
        s = seq(xmin, xmax, length = N)
        y1 = dgh(s, alpha, beta, delta, mu, lambda)
        y2 = pgh(s, alpha, beta, delta, mu, lambda)
        main1 = paste("GH Density\n",
            "alpha = ", as.character(alpha), " | ",
            "beta = ", as.character(beta), " | ",
            "delta = ", as.character(delta), " | ",
            "mu = ", as.character(mu), "|",
            "lambda = ", as.character(lambda) )
        main2 = paste("GH Probability\n",
            "xmin 0.01% = ", as.character(xmin), " | ",
            "xmax 0.99% = ", as.character(xmax))

        # Frame
        par(mfrow = c(2, 1), cex = 0.7)

        # Density:
        plot(s, y1, type = "l", xlim = c(xmin, xmax), col = "steelblue")
        abline (h = 0, lty = 3)
        title(main = main1)

        # Probability:
        plot(s, y2, type = "l", xlim = c(xmin, xmax), ylim = c(0, 1),
            col = "steelblue" )
        abline(h = 0.0, lty = 3)
        abline(h = 1.0, lty = 3)
        abline(h = 0.5, lty = 3)
        abline(v = mu, lty = 3, col = "red")
        title(main = main2)

        # Reset Frame:
        par(mfrow = c(1, 1), cex = 0.7)
    }

    # Open Slider Menu:
    .sliderMenu(refresh.code,
       names =       c( "N","alpha","beta","delta", "mu","lambda"),
       minima =      c(  50,  0.00, -2.00,   0.00, -5.0,  -4),
       maxima =      c(1000,  2.00, +2.00,   5.00, +5.0,   4),
       resolutions = c(  50,  0.20,  0.20,   1.00,  1.0,  0.5),
       starts =      c(  50,  1.00,  0.00,   1.00,  0.0,   1))
}


################################################################################

