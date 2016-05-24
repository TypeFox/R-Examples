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

HARMONIC_EPS <- 1e-10
EXPONENTIAL_EPS <- 1e-10

exponential.q <- function (qi, D, t)
{
    qi * exp(-D * t)
}

exponential.Np <- function (qi, D, t)
{
    if (D == 0)
        qi * t
    else
        qi / D * (1 - exp(-D * t))
}

harmonic.q <- function (qi, Di, t)
{
    qi / (1 + Di * t)
}

harmonic.Np <- function (qi, Di, t)
{
    if (Di == 0)
        qi * t
    else
        qi / Di * log(1 + Di * t)
}

harmonic.D <- function (Di, t)
{
    hyperbolic.D(Di, 1, t)
}

hyperbolic.q <- function (qi, Di, b, t)
{
    if (abs(b - 1) < HARMONIC_EPS)
        harmonic.q(qi, Di, t)
    else if (abs(b) < EXPONENTIAL_EPS)
        exponential.q(qi, Di, t)
    else
        qi * (1 + b * Di * t) ^ (-1/b)
}

hyperbolic.Np <- function (qi, Di, b, t)
{
    if (abs(b - 1) < HARMONIC_EPS)
        harmonic.Np(qi, Di, t)
    else if (abs(b) < EXPONENTIAL_EPS)
        exponential.Np(qi, Di, t)
    else if (Di == 0)
        qi * t
    else
        (qi / ((1 - b) * Di)) * (1 - (1 + b * Di * t) ^ (1 - (1/b)))
}

hyperbolic.D <- function (Di, b, t)
{
    Di / (1 + b * Di * t)
}

hyp2exp.transition <- function (Di, b, Df)
{
    if (Di < EXPONENTIAL_EPS || Df < EXPONENTIAL_EPS || b < HARMONIC_EPS)
        Inf
    else if (abs(Df - Di) < EXPONENTIAL_EPS)
        0
    else
        (Di / Df - 1) / (b * Di)
}

hyp2exp.q <- function (qi, Di, b, Df, t)
{
    t.trans <- hyp2exp.transition(Di, b, Df)
    q.trans <- hyperbolic.q(qi, Di, b, t.trans)

    q <- hyperbolic.q(qi, Di, b, t)
    q[t > t.trans] <- exponential.q(q.trans, Df, t[t > t.trans] - t.trans)

    q
}

hyp2exp.Np <- function (qi, Di, b, Df, t)
{
    t.trans <- hyp2exp.transition(Di, b, Df)
    q.trans <- hyperbolic.q(qi, Di, b, t.trans)
    Np.trans <- hyperbolic.Np(qi, Di, b, t.trans)

    Np <- hyperbolic.Np(qi, Di, b, t)
    Np[t > t.trans] <- Np.trans +
        exponential.Np(q.trans, Df, t[t > t.trans] - t.trans)

    Np
}

hyp2exp.D <- function (Di, b, Df, t)
{
    t.trans <- hyp2exp.transition(Di, b, Df)
    D <- hyperbolic.D(Di, b, t)
    D[t > t.trans] <- Df

    D
}

# from tangent effective
as.nominal <- function (D.eff,
                        from.period=c("year", "month", "day"),
                        to.period=c("year", "month", "day"))
{
    rescale.by.time(-log(1 - D.eff), from.period, to.period)
}

# to tangent effective
as.effective <- function (D.nom,
                          from.period=c("year", "month", "day"),
                          to.period=c("year", "month", "day"))
{
    1 - exp(-rescale.by.time(D.nom, from.period, to.period))
}

rescale.by.time <- function (value,
                             from.period=c("year", "month", "day"),
                             to.period=c("year", "month", "day"),
                             method=c("decline", "rate", "time"))
{
    from.period <- match.arg(from.period)
    to.period <- match.arg(to.period)
    method <- match.arg(method)

    if (method == "time") {
        tmp <- from.period
        from.period <- to.period
        to.period <- tmp
        stop.to <- "Invalid from.period."
        stop.from <- "Invalid to.period."
    } else if (method != "decline" && method != "rate") {
        stop("Invalid method.")
    } else {
        stop.to <- "Invalid to.period."
        stop.from <- "Invalid from.period."
    }

    if (from.period == to.period)
        value
    else if (from.period == "year") {
        if (to.period == "month")
            value / 12
        else if (to.period == "day")
            value / 365.25
        else
            stop(stop.to)
    } else if (from.period == "month") {
        if (to.period == "year")
            value * 12
        else if (to.period == "day")
            value / 30.4375
        else
            stop(stop.to)
    } else if (from.period == "day") {
        if (to.period == "year")
            value * 365.25
        else if (to.period == "month")
            value * 30.4375
        else
            stop(stop.to)
    } else
        stop(stop.from)
}
