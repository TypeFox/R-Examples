# rand.R -- functions for generating random variates
# Copyright (C) 2015 Matthew Clegg

#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  http://www.r-project.org/Licenses/

rpar <- function (n, rho, sigma_M, sigma_R, M0=0, R0=0, include.state = FALSE, robust=FALSE, nu=par.nu.default()) {
	# Generates a random partially AR(1) instance with parameters rho, sigma_M and sigma_R
    # In other words, generates a random realization of the sequence
    #   X_t = M_t + R_t
    #   M_t = rho M_{t-1} + epsilon_{M,t}
    #   R_t = R_{t-1} + epsilon_{R,t}
    #   epsilon_{M,t} ~ N(0, sigma_M^2)
    #   epsilon_{R,t} ~ N(0, sigma_R^2)
    # If include.state is FALSE, returns a vector of length n consisting of the
    # randomly generated values x_t.  If include.state is TRUE, returns an n x 5
    # matrix, whose columns are x, m, r, eps_M, eps_R

    if (!robust) {
        eps_M <- rnorm(n, 0, sigma_M)
        eps_R <- rnorm(n, 0, sigma_R)
    } else {
        eps_M <- rt(n, nu) * sigma_M
        eps_R <- rt(n, nu) * sigma_R
    }        
    M <- numeric (n)
    M[1] <- rho * M0 + eps_M[1]
    for (i in 2:n) {
        M[i] <- rho * M[i-1] + eps_M[i]
    }
    R <- cumsum(eps_R) + R0
    X <- M + R
    if (include.state) {
        return(data.frame(X=X, M=M, R=R, eps_M=eps_M, eps_R=eps_R))
    } else {
        return(X)
    }
}
