
# This R package is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             STABLE SLIDERS:
#  stableSlider          Displays stable distribution function
################################################################################


stableSlider <- 
    function(col= "steelblue", col.med = "gray30")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Displays the stable distribution

    # FUNCTION:

    # Internal Function:
    refresh.code <- function(...)
    {
        # Sliders:
        N     = .sliderMenu(no = 1)
        alpha = .sliderMenu(no = 2)
        beta  = .sliderMenu(no = 3)
        gamma = .sliderMenu(no = 4)
        delta = .sliderMenu(no = 5)
        pm    = .sliderMenu(no = 6)

        # Compute Data:
        x.rng <- round(
          stabledist::qstable(c(1,99)/100, alpha, beta, gamma, delta, pm),
          digits = 2)
        xmin <- x.rng[1]; xmax <- x.rng[2]
        s = seq(xmin, xmax, length = N)
        f = stabledist::dstable(s, alpha, beta, gamma, delta, pm)
        F = stabledist::pstable(s, alpha, beta, gamma, delta, pm)
        med <- stabledist::qstable(0.5, alpha, beta, gamma, delta, pm)
        main1 <- paste("Stable Density\n",
            "alpha = ", as.character(alpha), " | ",
            "beta = ", as.character(beta), " | ",
            "gamma = ", as.character(gamma), " | ",
            "delta = ", as.character(delta))
        main2 <- paste("Stable Probability\n",
            "xmin [ 1%] = ", as.character(xmin), " | ",
            "xmax [99%] = ", as.character(xmax), " | ",
            "pm = ", as.character(pm))

        # Frame:
        op <- par(mfrow = c(2, 1), cex = 0.7) ; on.exit(par(op))

        # Density:
        plot(s, f, type = "l", main = main1, xlim = x.rng, col = col)
        abline (h = 0, lty = 3)

        # Probability:
        plot(s, F, type = "l", xlim = x.rng, ylim = c(0, 1),
            col = col, main = main2)
        abline(h = 0:1, lty = 3)
        lines(c(par("usr")[1],med,med),
            c(0.5 ,0.5,  0), lty = 2, col=col.med)
        text(med, 0.1, "median", adj=0, col=col.med)
        axis(1, labels=expression(delta), at = delta,
            col = "red", col.axis="red", lwd=1.5,
            line = .5, tck = 1/8, hadj = -1, padj = -4)
    }

    # Open Slider Menu:
    .sliderMenu(refresh.code,
       names = c("no. points", "alpha","beta","gamma","delta", "pm"),
       minima =      c(   50,    0. ,  -1. ,    0. ,    -5. ,    0),
       maxima =      c( 1000,    2. ,  +1. ,    5. ,    +5. ,    2),
       resolutions = c(   50,    0.1,   0.1,    0.5,     0.5,    1),
       starts =      c(   50,    1.8,   0. ,    1. ,     0. ,    0))
}


################################################################################


