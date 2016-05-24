# Copyright (C) 2014 Mohammad H. Ferdosi
#
# HSPhase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# HSPhase program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http:#www.gnu.org/licenses/>.

rplot <- function(x, distance, start = 1, end = ncol(x), maximum = 100, overwrite = FALSE, method = "constant")
{
    if (!is.matrix(x)) 
        stop("The phase matrix must be matrix")
    y <- pm(bmh(x), method=method)
    y <- apply(y, 2, sum)
    if (missing(maximum))
    {
        maximum <- max(y) + 1
    }
    i <- which(y[start:end] < maximum)
    y <- y[i]
    if (overwrite == FALSE)
    {
        plot(distance[i], y, ylab = "Sum of the probability of recombinations", xlab = "Markers", col = "1", 
            type = "l")
    }
    else
    {
        color <- c("#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#A65628", "#F781BF", "#999999")
        lines(i, y, col = sample(color, 1), type = "l", lty = 4, lwd = 2)
    }
} 
