#    Function: Plot routine for object of class mrm.
#    Copyright (C) 2011  David Preinerstorfer
#    david.preinerstorfer@univie.ac.at
#
#    This program is free software; you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation; either version 2 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details. A copy may be obtained at
#    http://www.r-project.org/Licenses/

plot.mrm <- function(x, ...){
c <- length(x$class.size)
par(mfrow = c(1, 2))
matplot(x$beta, type = "o", main = "Item Easiness Parameters", xlab = "Items", ylab = "Parameters", pch = 20:(20+c), bg = 1:c)
abline(h = 0, lty = 2)
matplot(x$pi.r.c, type = "l", main = "Latent Score Probabilities", xlab = "Score", ylab = "")
legend("topright", legend = paste("class", 1:c, sep = " "), fill = 1:c)
}
