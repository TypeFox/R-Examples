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

sse <- function(q, forecast)
{
    sum((q - forecast) ^ 2)
}

best.exponential <- function(q, t,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10) # = 0.99995 / [time] effective
  )
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1])
                     # Di = decline from first to second point
                   ),

                 # cost function
                 function (guess) sse(q, exponential.q(guess[1], guess[2], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2]),
         sse=res$objective)
}

best.hyperbolic <- function(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2)  # b <= 2.0
  )
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5),   # right-ish for a lot of wells currently coming on

                 # cost function
                 function (guess)
                     sse(q, hyperbolic.q(guess[1], guess[2], guess[3], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
         sse=res$objective)
}

best.hyp2exp <- function(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35) # Df <= 0.35
  )
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1), # Df = about 9% effective

                 # cost function
                 function (guess)
                     sse(q,
                         hyp2exp.q(guess[1], guess[2], guess[3], guess[4], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1],
                              Di=res$par[2],
                              b=res$par[3],
                              Df=res$par[4]),
         sse=res$objective)
}

best.exponential.curtailed <- function(q, t,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0  # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    t[length(t)])
  )
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   t[2]  # t.curtail = second t in vector
                   ),

                 # cost function
                 function (guess)
                     sse(q,
                         curtailed.q(arps.decline(guess[1], guess[2]),
                                     guess[3], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(arps.decline(qi=res$par[1], Di=res$par[2]),
                         res$par[3]),
         sse=res$objective)
}

best.hyperbolic.curtailed <- function(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0,  # b > 0
    0   # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    t[length(t)])
  )
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # right-ish for a lot of wells currently coming on
                   t[2]  # t.curtail = second t in vector
                   ),

                 # cost function
                 function (guess)
                     sse(q,
                         curtailed.q(
                                 arps.decline(guess[1], guess[2], guess[3]),
                                 guess[4], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(
             arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
             res$par[4]),
         sse=res$objective)
}

best.hyp2exp.curtailed <- function(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0,  # Df > 0
    0   # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35, # Df <= 0.35
    t[length(t)])
  )
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1,  # Df = about 9% effective
                   t[2]
                   ),

                 # cost function
                 function (guess)
                     sse(q,
                         curtailed.q(
                             arps.decline(guess[1], guess[2], guess[3], guess[4]),
                             guess[5], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(arps.decline(qi=res$par[1],
                                      Di=res$par[2],
                                      b=res$par[3],
                                      Df=res$par[4]),
                         res$par[5]),
         sse=res$objective)
}

best.fit <- function(q, t)
{
    exp <- best.exponential(q, t)
    hyp <- best.hyperbolic(q, t)
    h2e <- best.hyp2exp(q, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.curtailed.fit <- function(q, t)
{
    exp <- best.exponential.curtailed(q, t)
    hyp <- best.hyperbolic.curtailed(q, t)
    h2e <- best.hyp2exp.curtailed(q, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.exponential.from.Np <- function(Np, t,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10) # = 0.99995 / [time] effective)
  )
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1] / t[1], # qi = rate(t = first t in vector)
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1])
                     # Di = decline from first to second point
                   ),

                 # cost function
                 function (guess) sse(Np, exponential.Np(guess[1], guess[2], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2]),
         sse=res$objective)
}

best.exponential.from.interval <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10) # = 0.99995 / [time] effective)
  )
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1])
                     # Di = decline from first to second point
                   ),

                 # cost function
                 function (guess) sse(volume, diff(exponential.Np(
                   guess[1], guess[2], c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2]),
         sse=res$objective)
}

best.hyperbolic.from.Np <- function(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    2)  # b <= 2.0
  )
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1] / t[1], # qi = rate(t = first t in vector)
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5),   # right-ish for a lot of wells currently coming on

                 # cost function
                 function (guess)
                     sse(Np, hyperbolic.Np(guess[1], guess[2], guess[3], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
         sse=res$objective)
}

best.hyperbolic.from.interval <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    2)  # b <= 2.0
  )
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5),   # right-ish for a lot of wells currently coming on

                 # cost function
                 function (guess) sse(volume,
                   diff(hyperbolic.Np(guess[1], guess[2], guess[3], c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
         sse=res$objective)
}

best.hyp2exp.from.Np <- function(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35) # Df <= 0.35
  )
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1] / t[1], # qi = rate(t = first t in vector)
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1), # Df = about 9% effective

                 # cost function
                 function (guess)
                     sse(Np,
                         hyp2exp.Np(guess[1], guess[2], guess[3], guess[4], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1],
                              Di=res$par[2],
                              b=res$par[3],
                              Df=res$par[4]),
         sse=res$objective)
}

best.hyp2exp.from.interval <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35) # Df <= 0.35
  )
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1), # Df = about 9% effective

                 # cost function
                 function (guess) sse(volume,
                   diff(hyp2exp.Np(guess[1], guess[2], guess[3], guess[4],
                     c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.decline(qi=res$par[1],
                              Di=res$par[2],
                              b=res$par[3],
                              Df=res$par[4]),
         sse=res$objective)
}

best.exponential.curtailed.from.Np <- function(Np, t,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0  # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    t[length(t)])
  )
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1] / t[1], # qi = rate(t = first t in vector)
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   t[2]  # t.curtail = second t in vector
                   ),

                 # cost function
                 function (guess)
                     sse(Np,
                         curtailed.Np(arps.decline(guess[1], guess[2]),
                                     guess[3], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(arps.decline(qi=res$par[1], Di=res$par[2]),
                         res$par[3]),
         sse=res$objective)
}

best.exponential.curtailed.from.interval <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0, # qi > 0
    0, # D > 0
    0  # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    t[length(t)])
  )
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   t[2]  # t.curtail = second t in vector
                   ),

                 # cost function
                 function (guess) sse(volume, diff(curtailed.Np(
                   arps.decline(guess[1], guess[2]), guess[3], c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(arps.decline(qi=res$par[1], Di=res$par[2]),
                         res$par[3]),
         sse=res$objective)
}

best.hyperbolic.curtailed.from.Np <- function(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0,  # b > 0
    0   # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    t[length(t)])
  )
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1] / t[1], # qi = rate(t = first t in vector)
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # right-ish for a lot of wells currently coming on
                   t[2]  # t.curtail = second t in vector
                   ),

                 # cost function
                 function (guess)
                     sse(Np,
                         curtailed.Np(
                                 arps.decline(guess[1], guess[2], guess[3]),
                                 guess[4], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(
             arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
             res$par[4]),
         sse=res$objective)
}

best.hyperbolic.curtailed.from.interval <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0,  # b > 0
    0   # t.curtail > 0
  ),
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    t[length(t)])
  )
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # right-ish for a lot of wells currently coming on
                   t[2]  # t.curtail = second t in vector
                   ),

                 # cost function
                 function (guess) sse(volume, diff(curtailed.Np(
                   arps.decline(guess[1], guess[2], guess[3]), guess[4],
                     c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(
             arps.decline(qi=res$par[1], Di=res$par[2], b=res$par[3]),
             res$par[4]),
         sse=res$objective)
}

best.hyp2exp.curtailed.from.Np <- function(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0,  # Df > 0
    0
  ),
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35, # Df <= 0.35
    t[length(t)])
  )
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    res <- nlminb(c( # initial guesses
                   Np[1] / t[1], # qi = rate(t = first t in vector)
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1,  # Df = about 9% effective
                   t[2]
                   ),

                 # cost function
                 function (guess)
                     sse(Np,
                         curtailed.Np(
                             arps.decline(guess[1], guess[2], guess[3], guess[4]),
                             guess[5], t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(arps.decline(qi=res$par[1],
                                      Di=res$par[2],
                                      b=res$par[3],
                                      Df=res$par[4]),
                         res$par[5]),
         sse=res$objective)
}

best.hyp2exp.curtailed.from.interval <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0,  # Df > 0
    0
  ),
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35, # Df <= 0.35
    t[length(t)])
  )
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1,  # Df = about 9% effective
                   t[2]
                   ),

                 # cost function
                 function (guess) sse(volume, diff(curtailed.Np(
                   arps.decline(guess[1], guess[2], guess[3], guess[4]), guess[5],
                     c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=curtail(arps.decline(qi=res$par[1],
                                      Di=res$par[2],
                                      b=res$par[3],
                                      Df=res$par[4]),
                         res$par[5]),
         sse=res$objective)
}

best.fit.from.Np <- function(Np, t)
{
    exp <- best.exponential.from.Np(Np, t)
    hyp <- best.hyperbolic.from.Np(Np, t)
    h2e <- best.hyp2exp.from.Np(Np, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.fit.from.interval <- function(volume, t, t.begin=0.0)
{
    exp <- best.exponential.from.interval(volume, t, t.begin)
    hyp <- best.hyperbolic.from.interval(volume, t, t.begin)
    h2e <- best.hyp2exp.from.interval(volume, t, t.begin)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.curtailed.fit.from.Np <- function(Np, t)
{
    exp <- best.exponential.curtailed.from.Np(Np, t)
    hyp <- best.hyperbolic.curtailed.from.Np(Np, t)
    h2e <- best.hyp2exp.curtailed.from.Np(Np, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.curtailed.fit.from.interval <- function(volume, t, t.begin=0.0)
{
    exp <- best.exponential.curtailed.from.interval(volume, t, t.begin)
    hyp <- best.hyperbolic.curtailed.from.interval(volume, t, t.begin)
    h2e <- best.hyp2exp.curtailed.from.interval(volume, t, t.begin)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.exponential.with.buildup <- function(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10),        # = 0.99995 / [time] effective
  initial.rate=q[1], time.to.peak=t[which.max(q)])
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1])
                     # Di = decline from first to second point
                   ),

                 # cost function
                 function (guess) sse(q, arps.q(arps.with.buildup(
                   arps.decline(guess[1], guess[2]), initial.rate,
                   time.to.peak), t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2]),
           initial.rate, time.to.peak), sse=res$objective)
}

best.hyperbolic.with.buildup <- function(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2), # b <= 2.0
  initial.rate=q[1], time.to.peak=t[which.max(q)])
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5),   # right-ish for a lot of wells currently coming on

                 # cost function
                 function (guess) sse(q, arps.q(arps.with.buildup(
                   arps.decline(guess[1], guess[2], guess[3]),
                   initial.rate, time.to.peak), t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2],
           b=res$par[3]), initial.rate, time.to.peak), sse=res$objective)
}

best.hyp2exp.with.buildup <- function(q, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(q) * 5, # qi < qmax * 5
    10, # = 0.99995 / [time] effective
    2,  # b <= 2.0
    0.35), # Df <= 0.35
  initial.rate=q[1], time.to.peak=t[which.max(q)])
{
    if (length(q) != length(t) || length(q) <= 2)
        stop("Invalid lengths for q, t vectors.")

    res <- nlminb(c( # initial guesses
                   q[1], # qi = q(t = first t in vector)
                   -((q[2] - q[1]) / q[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1), # Df = about 9% effective

                 # cost function
                 function (guess) sse(q, arps.q(arps.with.buildup(
                   arps.decline(guess[1], guess[2], guess[3], guess[4]),
                   initial.rate, time.to.peak), t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2],
           b=res$par[3], Df=res$par[4]), initial.rate, time.to.peak),
         sse=res$objective)
}

best.fit.with.buildup <- function(q, t)
{
    exp <- best.exponential.with.buildup(q, t)
    hyp <- best.hyperbolic.with.buildup(q, t)
    h2e <- best.hyp2exp.with.buildup(q, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.exponential.from.Np.with.buildup <- function(Np, t,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10), # = 0.99995 / [time] effective
  initial.rate=Np[1] / t[1],
  time.to.peak=(t[which.max(diff(Np))] + t[which.max(diff(Np)) + 1]) / 2.0)
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    max.delta <- max(Np[1], diff(Np))
    t.max.delta <- ifelse(max.delta == Np[1], t[1],
      t[which.max(diff(Np)) + 1] - t[which.max(diff(Np))])

    res <- nlminb(c( # initial guesses
                   max.delta / t.max.delta, # qi = max Np delta / time delta
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1])
                     # Di = decline from first to second point
                   ),

                 # cost function
                 function (guess) sse(Np, arps.Np(arps.with.buildup(
                   arps.decline(guess[1], guess[2]), initial.rate,
                   time.to.peak), t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2]),
           initial.rate, time.to.peak), sse=res$objective)
}

best.exponential.from.interval.with.buildup <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0, # qi > 0
    0), # D > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10), # = 0.99995 / [time] effective
  initial.rate=volume[1] / (t[1] - t.begin),
  time.to.peak=(t - diff(c(t.begin, t)) / 2)[which.max(volume)])
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1])
                     # Di = decline from first to second point
                   ),

                 # cost function
                 function (guess) sse(volume, diff(arps.Np(arps.with.buildup(
                   arps.decline(guess[1], guess[2]), initial.rate,
                   time.to.peak), c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2]),
           initial.rate, time.to.peak), sse=res$objective)
}

best.hyperbolic.from.Np.with.buildup <- function(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    2), # b <= 2.0
  initial.rate=Np[1] / t[1],
  time.to.peak=(t[which.max(diff(Np))] + t[which.max(diff(Np)) + 1]) / 2.0)
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    max.delta <- max(Np[1], diff(Np))
    t.max.delta <- ifelse(max.delta == Np[1], t[1],
      t[which.max(diff(Np)) + 1] - t[which.max(diff(Np))])

    res <- nlminb(c( # initial guesses
                   max.delta / t.max.delta, # qi = max Np delta / time delta
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5),   # right-ish for a lot of wells currently coming on

                 # cost function
                 function (guess) sse(Np, arps.Np(arps.with.buildup(
                   arps.decline(guess[1], guess[2], guess[3]),
                   initial.rate, time.to.peak), t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2],
           b=res$par[3]), initial.rate, time.to.peak), sse=res$objective)
}

best.hyperbolic.from.interval.with.buildup <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0,  # Di > 0
    0), # b > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    2), # b <= 2.0
  initial.rate=volume[1] / (t[1] - t.begin),
  time.to.peak=(t - diff(c(t.begin, t)) / 2)[which.max(volume)])
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5),   # right-ish for a lot of wells currently coming on

                 # cost function
                 function (guess) sse(volume, diff(arps.Np(arps.with.buildup(
                   arps.decline(guess[1], guess[2], guess[3]),
                   initial.rate, time.to.peak), c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2],
           b=res$par[3]), initial.rate, time.to.peak), sse=res$objective)
}

best.hyp2exp.from.Np.with.buildup <- function(Np, t,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(c(Np[1], diff(Np)) / diff(c(0, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35), # Df <= 0.35
  initial.rate=Np[1] / t[1],
  time.to.peak=(t[which.max(diff(Np))] + t[which.max(diff(Np)) + 1]) / 2.0)
{
    # drop leading zero records
    which.nz <- Np != 0
    Np <- Np[which.nz]
    t <- t[which.nz]

    if (length(Np) != length(t) || length(Np) <= 2)
        stop("Invalid lengths for Np, t vectors.")

    max.delta <- max(Np[1], diff(Np))
    t.max.delta <- ifelse(max.delta == Np[1], t[1],
      t[which.max(diff(Np)) + 1] - t[which.max(diff(Np))])

    res <- nlminb(c( # initial guesses
                   max.delta / t.max.delta, # qi = max Np delta / time delta
                   -(((Np[2] - Np[1]) - Np[1]) / Np[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1), # Df = about 9% effective

                 # cost function
                 function (guess) sse(Np, arps.Np(arps.with.buildup(
                   arps.decline(guess[1], guess[2], guess[3], guess[4]),
                   initial.rate, time.to.peak), t)),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2],
           b=res$par[3], Df=res$par[4]), initial.rate, time.to.peak),
         sse=res$objective)
}

best.hyp2exp.from.interval.with.buildup <- function(volume, t, t.begin=0.0,
  lower=c( # lower bounds
    0,  # qi > 0
    0.35,  # Di > 0
    0,  # b > 0
    0), # Df > 0
  upper=c( # upper bounds
    max(volume / diff(c(t.begin, t))) * 5, # qi < max(rate) * 5
    10, # = 0.99995 / [time] effective
    5,  # b <= 2.0
    0.35), # Df <= 0.35
  initial.rate=volume[1] / (t[1] - t.begin),
  time.to.peak=(t - diff(c(t.begin, t)) / 2)[which.max(volume)])
{
    if (length(volume) != length(t) || length(volume) <= 2)
        stop("Invalid lengths for volume, t vectors.")

    res <- nlminb(c( # initial guesses
                   volume[1] / (t[1] - t.begin), # qi = rate(t = first t in vector)
                   -((volume[2] - volume[1]) / volume[1]) / (t[2] - t[1]),
                     # Di = decline from first to second point
                   1.5,  # b = right-ish for a lot of wells currently coming on
                   0.1), # Df = about 9% effective

                 # cost function
                 function (guess) sse(volume, diff(arps.Np(arps.with.buildup(
                   arps.decline(guess[1], guess[2], guess[3], guess[4]),
                   initial.rate, time.to.peak), c(t.begin, t)))),

                 # bounds
                 lower=lower, upper=upper)

    list(decline=arps.with.buildup(arps.decline(qi=res$par[1], Di=res$par[2],
           b=res$par[3], Df=res$par[4]), initial.rate, time.to.peak),
         sse=res$objective)
}

best.fit.from.Np.with.buildup <- function(Np, t)
{
    exp <- best.exponential.from.Np.with.buildup(Np, t)
    hyp <- best.hyperbolic.from.Np.with.buildup(Np, t)
    h2e <- best.hyp2exp.from.Np.with.buildup(Np, t)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}

best.fit.from.interval.with.buildup <- function(volume, t, t.begin=0.0)
{
    exp <- best.exponential.from.interval.with.buildup(volume, t, t.begin)
    hyp <- best.hyperbolic.from.interval.with.buildup(volume, t, t.begin)
    h2e <- best.hyp2exp.from.interval.with.buildup(volume, t, t.begin)

    if (exp$sse <= hyp$sse && exp$sse <= h2e$sse)
        exp
    else if (hyp$sse <= exp$sse && hyp$sse <= h2e$sse)
        hyp
    else
        h2e
}
