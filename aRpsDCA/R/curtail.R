# aRpsDCA
# Copyright (C) 2016 dwt | terminus data science, LLC
# <dwt [at] terminusdatascience.com>

# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.

# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.

# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301
# USA

curtailed.q <- function (decl, t.curtail, t)
{
    q <- rep(decl$qi, length(t))
    q[t > t.curtail] <- arps.q(decl, t[t > t.curtail] - t.curtail)
    q
}

curtailed.Np <- function (decl, t.curtail, t)
{
    Np <- decl$qi * t
    Np[t > t.curtail] <- decl$qi * t.curtail +
      arps.Np(decl, t[t > t.curtail] - t.curtail)
    Np
}

curtailed.D <- function (decl, t.curtail, t)
{
    D <- rep(0, length(t))
    D[t > t.curtail] <- arps.D(decl, t[t > t.curtail] - t.curtail)
    D
}

curtail <- function (decl, t.curtail)
{
    res <- list(arps=decl, t.curtail=t.curtail)
    class(res) <- c("curtailed", "arps")
    res
}

arps.q.curtailed <- function(decl, t)
    curtailed.q(decl$arps, decl$t.curtail, t)

arps.Np.curtailed <- function(decl, t)
    curtailed.Np(decl$arps, decl$t.curtail, t)

arps.D.curtailed <- function(decl, t)
    curtailed.D(decl$arps, decl$t.curtail, t)

format.curtailed <- function(x, ...)
{
    paste("Curtailed ",
          format(x$arps, ...),
          " with t.curtail = ",
          format(x$t.curtail, ...),
          sep="")
}
