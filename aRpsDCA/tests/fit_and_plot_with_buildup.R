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

library(aRpsDCA)

fitme.exponential.t <- seq(0, 5, 1 / 12) # 5 years
fitme.exponential.q <- arps.q(arps.with.buildup(arps.decline(
    1000, # Bbl/d
    as.nominal(0.70) # / year
    ), 250, 1.5 / 12),
    fitme.exponential.t
) * rnorm(n=length(fitme.exponential.t), mean=1, sd=0.1) # perturb

exponential.fit <- best.exponential.with.buildup(fitme.exponential.q, fitme.exponential.t)
cat(paste("SSE:", exponential.fit$sse))
dev.new()
plot(fitme.exponential.q ~ fitme.exponential.t, main="Exponential Fit",
     col="blue", log="y", xlab="Time", ylab="Rate")
lines(arps.q(exponential.fit$decline, fitme.exponential.t) ~ fitme.exponential.t,
      col="red")
legend("topright", pch=c(1, NA), lty=c(NA, 1), col=c("blue", "red"), legend=c("Actual", "Fit"))

fitme.hyperbolic.t <- seq(0, 5, 1 / 12) # 5 years
fitme.hyperbolic.q <- arps.q(arps.with.buildup(arps.decline(
    1000, # Bbl/d
    as.nominal(0.70), # / year
    1.9
    ), 250, 1.5 / 12),
    fitme.hyperbolic.t
) * rnorm(n=length(fitme.hyperbolic.t), mean=1, sd=0.1) # perturb

hyperbolic.fit <- best.hyperbolic.with.buildup(fitme.hyperbolic.q, fitme.hyperbolic.t)
cat(paste("SSE:", hyperbolic.fit$sse))
dev.new()
plot(fitme.hyperbolic.q ~ fitme.hyperbolic.t, main="Hyperbolic Fit",
     col="blue", log="y", xlab="Time", ylab="Rate")
lines(arps.q(hyperbolic.fit$decline, fitme.hyperbolic.t) ~ fitme.hyperbolic.t,
      col="red")
legend("topright", pch=c(1, NA), lty=c(NA, 1), col=c("blue", "red"), legend=c("Actual", "Fit"))

fitme.hyp2exp.t <- seq(0, 5, 1 / 12) # 5 years
fitme.hyp2exp.q <- arps.q(arps.with.buildup(arps.decline(
    1000, # Bbl/d
    as.nominal(0.70), # / year
    1.9,
    as.nominal(0.15) # / year
    ), 250, 1.5 / 12),
    fitme.hyp2exp.t
) * rnorm(n=length(fitme.hyp2exp.t), mean=1, sd=0.1) # perturb

hyp2exp.fit <- best.hyp2exp.with.buildup(fitme.hyp2exp.q, fitme.hyp2exp.t)
cat(paste("SSE:", hyp2exp.fit$sse))
dev.new()
plot(fitme.hyp2exp.q ~ fitme.hyp2exp.t, main="Hyperbolic-to-Exponential Fit",
     col="blue", log="y", xlab="Time", ylab="Rate")
lines(arps.q(hyp2exp.fit$decline, fitme.hyp2exp.t) ~ fitme.hyp2exp.t,
      col="red")
legend("topright", pch=c(1, NA), lty=c(NA, 1), col=c("blue", "red"), legend=c("Actual", "Fit"))

overall.best <- best.fit.with.buildup(fitme.hyp2exp.q, fitme.hyp2exp.t)
cat(paste("SSE:", overall.best$sse))
dev.new()
plot(fitme.hyp2exp.q ~ fitme.hyp2exp.t, main="Overall Best Fit (h2e Data)",
     col="blue", log="y", xlab="Time", ylab="Rate")
lines(arps.q(overall.best$decline, fitme.hyp2exp.t) ~ fitme.hyp2exp.t,
      col="red")
legend("topright", pch=c(1, NA), lty=c(NA, 1), col=c("blue", "red"), legend=c("Actual", "Fit"))

