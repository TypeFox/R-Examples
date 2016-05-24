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

arps.decline <- function(qi, Di, b=NA, Df=NA)
{
    if (is.na(b) && !is.na(Df))
        stop("Must specify b for hyp2exp decline.")

    if (is.na(b)) {
        res <- list(qi=qi, D=Di)
        class(res) <- c("exponential", "arps")
    } else if (is.na(Df)) {
        res <- list(qi=qi, Di=Di, b=b)
        class(res) <- c("hyperbolic", "arps")
    } else {
        res <- list(qi=qi, Di=Di, b=b, Df=Df)
        class(res) <- c("hyp2exp", "arps")
    }

    res
}

arps.q <- function(decl, t)
{
    UseMethod("arps.q")
}

arps.Np <- function(decl, t)
{
    UseMethod("arps.Np")
}

arps.D <- function(decl, t)
{
    UseMethod("arps.D")
}

arps.q.arps <- function(decl, t) stop("Need specific decline class.")
arps.Np.arps <- function(decl, t) stop("Need specific decline class.")
arps.D.arps <- function(decl, t) stop("Need specific decline class.")

arps.q.exponential <- function(decl, t) exponential.q(decl$qi, decl$D, t)
arps.q.hyperbolic <- function(decl, t) hyperbolic.q(decl$qi, decl$Di, decl$b, t)
arps.q.hyp2exp <- function(decl, t) hyp2exp.q(decl$qi, decl$Di, decl$b, decl$Df, t)

arps.Np.exponential <- function(decl, t) exponential.Np(decl$qi, decl$D, t)
arps.Np.hyperbolic <- function(decl, t) hyperbolic.Np(decl$qi, decl$Di, decl$b, t)
arps.Np.hyp2exp <- function(decl, t) hyp2exp.Np(decl$qi, decl$Di, decl$b, decl$Df, t)

arps.D.exponential <- function(decl, t) decl$D
arps.D.hyperbolic <- function(decl, t) hyperbolic.D(decl$Di, decl$b, t)
arps.D.hyp2exp <- function(decl, t) hyp2exp.D(decl$Di, decl$b, decl$Df, t)

format.arps <- function(x, ...)
{
    paste("Arps decline:", format(unclass(x), ...), sep="\n")
}

format.exponential <- function(x, ...)
{
    paste("Arps exponential decline: <qi = ",
          format(x$qi, ...),
          ", D = ",
          format(x$D, ...),
          ">",
          sep="")
}

format.hyperbolic <- function(x, ...)
{
    paste("Arps hyperbolic decline: <qi = ",
          format(x$qi, ...),
          ", Di = ",
          format(x$Di, ...),
          ", b = ",
          format(x$b, ...),
          ">",
          sep="")
}

format.hyp2exp <- function(x, ...)
{
    paste("Arps hyperbolic-to-exponential decline: <qi = ",
          format(x$qi, ...),
          ", Di = ",
          format(x$Di, ...),
          ", b = ",
          format(x$b, ...),
          ", Df = ",
          format(x$Df, ...),
          ">",
          sep="")
}

print.arps <- function(x, ...)
{
    print(format(x, ...))
    invisible(x)
}

arps.with.buildup <- function(decl, initial.rate, time.to.peak)
{
    res <- list(decline=decl, 
                initial.rate=initial.rate,
                time.to.peak=time.to.peak,
                peak.rate=arps.q(decl, time.to.peak))
    class(res) <- c("buildup", "arps")
    res
}

arps.q.buildup <- function(decl, t)
{
    which.buildup <- which(t <= decl$time.to.peak)
    res <- numeric(length(t))
    buildup.m <- (decl$peak.rate - decl$initial.rate) / decl$time.to.peak
    res[which.buildup] <- buildup.m * t[which.buildup] + decl$initial.rate
    res[-which.buildup] <- arps.q(decl$decline, t[-which.buildup])
    res
}

arps.Np.buildup <- function(decl, t)
{
    which.buildup <- which(t <= decl$time.to.peak)
    res <- numeric(length(t))
    buildup.m <- (decl$peak.rate - decl$initial.rate) / decl$time.to.peak
    res[which.buildup] <- buildup.m * t[which.buildup]^2 * 0.5 +
      decl$initial.rate * t[which.buildup]
    res[-which.buildup] <- (arps.Np(decl$decline, t[-which.buildup])
      - arps.Np(decl$decline, decl$time.to.peak)
      + (buildup.m * decl$time.to.peak^2 * 0.5 +
         decl$initial.rate * decl$time.to.peak))
    res
}

arps.D.buildup <- function(decl, t)
{
    which.buildup <- which(t < decl$time.to.peak)
    res <- numeric(length(t))
    res[which.buildup] <- NA
    res[-which.buildup] <- arps.D(decl$decline, t[-which.buildup])
    res
}

format.buildup <- function(x, ...)
{
    paste(format(x$decline, ...),
          " with buildup: <initial rate = ",
          format(x$initial.rate, ...),
          ", time to peak = ",
          format(x$time.to.peak, ...),
          ", peak rate = ",
          format(x$peak.rate, ...),
          ">",
          sep="")
}
