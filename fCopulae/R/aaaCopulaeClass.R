
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################
# FUNCTION:                 COPULA SPECIFICATION:
#  fCOPULA                   S4 class representation
#  show                      S4 print method for copula specification
# FUNCTION:                 FRECHET COPULA:
#  pfrechetCopula            Computes Frechet copula probability
# FUNCTION:                 SPEARMAN'S RHO:
#  .copulaRho                Spearman's rho by integration for "ANY" copula
################################################################################


################################################################################



setClass("fCOPULA",
    # Description:
    #   Specifying and creating copula objects
    
    # Copula Representation:
    representation(
        call = "call",
        copula = "character",
        param = "list",
        title = "character",
        description = "character")
)


# ------------------------------------------------------------------------------


setMethod("show", "fCOPULA",
    function(object)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Print and Summary method for fCOPULA

    # Source:
    #   This function copies code from base:print.htest

    # FUNCTION:

    # Unlike print the argument for show is 'object'.
    x = object

    # Title:
    cat("\nTitle:\n ", x@title, "\n", sep = "")

    # Call:
    cat("\nCall:\n ")
    cat(paste(deparse(x@call), sep = "\n", collapse = "\n"), "\n", sep = "")

    # Copula Type:
    cat("\nCopula:\n ", x@copula, "\n", sep = "")

    # Model Parameter:
    if (length(x@param) != 0) {
        cat("\nModel Parameter(s):\n ")
        print(unlist(x@param), quote = FALSE)
    }

    # Description:
    cat("\nDescription:\n ", x@description, sep = "")
    cat("\n\n")

    # Return Value:
    invisible(object)
})


################################################################################
# Frechet Copulae:


pfrechetCopula <-
    function(u = 0.5, v = u, type = c("m", "pi", "w"),
    output = c("vector", "list"))
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes Frechet copula probability

    # Arguments:
    #   u, v - two numeric values or vectors of the same length at
    #       which the copula will be computed. If 'u' is a list then the
    #       the '$x' and '$y' elements will be used as 'u' and 'v'.
    #       If 'u' is a two column matrix then the first column will
    #       be used as 'u' and the the second as 'v'.
    #   type - the type of the Frechet copula. A character
    #       string selected from: "m", "pi", or "w".
    #   output - a character string specifying how the output should
    #       be formatted. By default a vector of the same length as
    #       'u' and 'v'. If specified as "list" then 'u' and 'v' are
    #       expected to span a two-dimensional grid as outputted by the
    #       function 'grid2d' and the function returns a list with
    #       elements '$x', 'y', and 'z' which can be directly used
    #       for example by 2D plotting functions.

    # Examples:
    #   persp(pfrechetCopula(u=grid2d(), output="list", type = "m"))
    #   persp(pfrechetCopula(u=grid2d(), output="list", type = "pi"))
    #   persp(pfrechetCopula(u=grid2d(), output="list", type = "w"))

    # FUNCTION:

    # Match Arguments:
    type = type[1] # Allow for "psp" ... # type = match.arg(type)
    output = match.arg(output)

    # Settings:
    if (is.list(u)) {
        v = u[[2]]
        u = u[[1]]
    }
    if (is.matrix(u)) {
        v = u[, 1]
        u = u[, 2]
    }

    # Compute Copula Probability:
    if (type == "m") {
        # C(u,v) = min(u,v)
        C.uv = apply(cbind(u, v), 1, min)
    }
    if (type == "pi") {
        # C(u, v) = u * v
        C.uv = u * v
    }
    if (type == "w") {
        # C(u,v) = max(u+v-1, 0)
        C.uv = apply(cbind(X = u+v-1, Y = rep(0, length = length(u))), 1, max)
    }
    if (type == "psp") {
        # C(u,v) = u*v/(u+v-u*v)
        C.uv = u*v/(u+v-u*v)
    }

    # Add Control:
    attr(C.uv, "control") <- unlist(list(type = type))

    # As List ?
    if (output == "list") {
        N = sqrt(length(u))
        x = u[1:N]
        y = matrix(v, ncol = N)[1, ]
        C.uv = list(x = x, y = y,  z = matrix(C.uv, ncol = N))
    }

    # Return Value:
    C.uv
}


################################################################################


.copulaRho =
    function(rho = NULL, alpha = NULL, param = NULL,
    family = c("elliptical", "archm", "ev", "archmax"),
    type = NULL, error = 1e-3, ...)
{   
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Spearman's rho by integration for "ANY" copula

    # Notes:
    #   pellipticalCopula(u, v, rho,   param, type, output, border)
    #   parchmCopula     (u, v, alpha,        type, output, alternative)
    #   pevCopula        (u, v,        param, type, output, alternative)
    #   parchmaxCopula   (u, v,        param, type, output )

    # Examples:
    #   .copulaRho(rho = 0.5, family = "elliptical", type = "norm")
    #   .copulaRho(alpha = 1, family = "archm", type = "1")
    #   .copulaRho(param = 2, family = "ev", type = "galambos")

    # FUNCTION:

    # Match Arguments:
    family = match.arg(family)

    # Type:
    if (is.null(type)) {
        family = "elliptical"
        type = "norm"
    } else {
        type = as.character(type)
    }

    # 2D Function to be integrated:
    rho <<- rho
    alpha <<- alpha
    param <<- param
    type <<- type
    if (family == "elliptical") {
        dCopulaRho <- function(x, y) {
            C = pellipticalCopula(x, y, rho = rho, param = param, type = type)
            12 * (C - x*y )
        }
    } else if (family == "archm") {
        if (is.null(alpha)) alpha <<- archmParam(type)$param
        check = archmCheck(alpha, type)
        dCopulaRho <- function(x, y) {
            C = parchmCopula(x, y, alpha = alpha, type = type)
            12 * (C - x*y )
        }
    } else if (family == "ev") {
        dCopulaRho <- function(x, y) {
            C = pevCopula(x, y, param = param, type = type)
            12 * (C - x*y )
        }
    }
    # else if (family == "archmax") {
    #     dCopulaRho <- function(x, y) {
    #         C = parchmaxCopula(x, y, param = param, type = type)
    #         12 * (C - x*y )
    #     }
    # }

    # Integrate:
    ans = integrate2d(dCopulaRho, error = error)
    Rho = ans$value
    error = ans$error

    # Result:
    control = list(rho = rho, alpha = alpha, param = param,
        family = family, type = type, error = signif(error, 3))
    attr(Rho, "control") <- unlist(control)

    # Return Value:
    Rho
}


################################################################################

