# Copyright 2014 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


firthglm.family <- function(family)
{
    if (!(family$family %in% c("binomial", "gaussian", "Gamma",
                               "inverse.gaussian", "poisson",
                               "quasibinomial", "quasipoisson"))) {
        print(family)
        stop("'family' not supported")
    }

    # make sure link is canonical
    if (!switch(family$family,
                binomial = family$link == "logit",
                gaussian = family$link == "identity",
                Gamma = family$link == "inverse",
                inverse.gaussian = family$link == "1/mu^2",
                poisson = family$link == "log",
                quasibinomial = family$link == "logit",
                quasipoisson = family$link == "log")) {
             print(family)
             warning("noncanonical link functions are not supported")
    }

    # add skewness, kurtosis, and entropy functions
    if (family$family == "binomial" || family$family == "quasibinomial") {
        if (family$link == "logit") {
            family$linkinv <- function(eta) {
                eta.min <- log(.Machine$double.xmin)
                eta.max <- log(2^(.Machine$double.digits) - 1)

                eta <- pmin(pmax(eta, eta.min), eta.max)
                tmp <- exp(eta)
                tmp / (1 + tmp)
            }
            family$mu.eta <- function(eta) {
                eta.min <- log(.Machine$double.xmin)
                eta.max <- log(2^(.Machine$double.digits) - 1)

                eta <- pmin(pmax(eta, eta.min), eta.max)
                opexp <- 1 + exp(eta)
                exp(eta) / (opexp * opexp)
            }
        }
        family$skewness <- function(mu) {
            (1 - 2 * mu) / sqrt(mu * (1 - mu))
        }
        family$kurtosis <- function(mu) {
            (1 - 6 * mu * (1 - mu)) / (mu * (1 - mu))
        }
        family$entropy <- function(mu) {
            ifelse(mu == 0 | mu == 1, 0, mu * log(mu) + (1 - mu) * log(1-mu))
        }
    } else if (family$family == "gaussian") {
        family$skewness <- function(mu) {
            (mu[!is.na(mu)] <- 0)
        }
        family$kurtosis <- function(mu) {
            (mu[!is.na(mu)] <- 0)
        }
        family$entropy <- function(mu) {
            mu^2 / 2
        }
    } else if (family$family == "Gamma") {
        family$skewness <- function(mu) {
            (mu[!is.na(mu)] <- 2)
        }
        family$kurtosis <- function(mu) {
            (mu[!is.na(mu)] <- 6)
        }
        family$entropy <- function(mu) {
            -1 - log(mu)
        }
    } else if (family$family == "inverse.gaussian") {
        family$skewness <- function(mu) {
            3 * sqrt(mu)
        }
        family$kurtosis <- function(mu) {
            15 * mu
        }
        family$entropy <- function(mu) {
            1 / (2 * mu)
        }
    } else if (family$family == "poisson" || family$family == "quasipoisson") {
        family$skewness <- function(mu) {
            1/sqrt(mu)
        }
        family$kurtosis <- function(mu) {
            1/mu
        }
        family$entropy <- function(mu) {
            ifelse(mu == 0, 0, mu * log(mu) - mu)
        }
    } else {
        print(family)
        stop("unsupported family")
    }

    family
}


