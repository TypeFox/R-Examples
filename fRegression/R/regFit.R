
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


###############################################################################
# FUNCTION:             PARAMETER ESTIMATION:
#  regFit                Wrapper Function for Regression Models
#  .lmFit                 Linear Regression Model
#  .rlmFit                Robust Linear Regression Model
#  .glmFit                Generalized Linear Model
#  .gamFit                Generalized Additive Model
#  .pprFit                Projection Pursuit Regression Model
#  .nnetFit               Feedforward Neural Network Model
#  .polymarsFit           Polytochomous MARS Model
###############################################################################


###############################################################################
# MODEL:        PACKAGE     print   plot   summary   print     predict
#                                    persp           summary
#   lm          stats       x       x      x         x         x
#   rlm         MASS
#   glm         stats       x       -      x         x         x
#   gam         mgcv        x       x      x         x         x
#   ppr         stats       x       x      x         x         x
#   nnet        nnet        x       -      x         x         x
#   polymars*   polspline   -       xx     x         -         x
###############################################################################


regFit <-
    function (formula, data, family = gaussian,
    use = c("lm", "rlm", "glm", "gam", "ppr", "nnet", "polymars"),
    title = NULL, description = NULL, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Common function call for several selected regression models.

    # Details:
    #   This is a wrapper function for the following regrssion models:
    #     LM          Linear Regression Modelling
    #     RLM         Robust Linear Regression Modelling
    #     GLM         Generalized Linear Modelling
    #     GAM         Generalized Additive Modelling
    #     PPR         Projection Pursuit Regression
    #     NNET        Feedforward Neural Net
    #     POLYMARS    Polytochomous MARS Modeling

    # Notes:
    #   Available Methods are
    #   "print", "plot", "summary", "predict"
    #   "coef", "formula", "residuals" "fitted", "vcov"

    # Example:
    #   regFit(Y ~ X1 + X2 + X3, regSim())

    # FUNCTION:

    # Match Arguments:
    use <- match.arg(use)
    if (missing(data)) data <- NULL

    # Transform data into a dataframe
    if (!is.null(data)) {
      Data <- if (inherits(data, "timeSeries")) data else as.timeSeries(data)
      data <- as.data.frame(data)
    } else {
      Data <- data <- NULL
    }

    # Function to be called:
    fun <- paste(".", match.arg(use), sep = "")

    # Title:
    if (is.null(title)) {
        if (use == "lm") title = "Linear Regression Modeling"
        if (use == "rlm") title = "Robust Linear Regression Modeling"
        if (use == "glm") title = "Generalized Linear Modeling"
        if (use == "gam") title = "Generalized Additive Modeling"
        if (use == "ppr") title = "Projection Pursuit Regression"
        if (use == "nnet") title = "Feedforward Neural Network Modeling"
        if (use == "polymars") title = "Polytochomous MARS Modeling"
    }

    # Description:
    if (is.null(description)) description = description()

    # Compose Command to be Called:
    cmd <- match.call()
    if (!is.null(cmd$use)) cmd = cmd[-match("use", names(cmd), 0)]
    cmd[[1]] <- as.name(fun)
    # Use this to access hidden functions in a parent frame:
    #cmd[[1]] <- substitute(fRegression:::f, list(f=as.name(fun)))

    # Ensure that data is a data.frame
    if (!is.null(cmd$data)) cmd$data <- as.name("data")
    # Use this to directly pass the argument from the parent frame:
    #if (!is.null(cmd$data)) cmd$data <- call("as.data.frame", cmd$data)

    # Fit Regression Model:
    fit <- eval(cmd)
    # Use this to evaluate in parent frame:
    #fit <- eval(cmd, parent.frame())

    # Add "cmd" to Fit:
    fit$cmd <- cmd

    # Add "xlevels" to Fit (if missing):
    if (is.null(fit$xlevels)) fit$xlevels = list()

    # Add "residuals" and "fitted" to Fit (to be sure ...):
    fit$residuals <- as.vector(fit$residuals)
    fit$fitted.values <- as.vector(fit$fitted.values)

    # Add "parameters" as Alternative:
    fit$parameters <- fit$coef

    # Extend to class "list":
    class(fit) <- c("list", class(fit))
    if (!inherits(fit, "lm")) class(fit) = c(class(fit), "lm")

    # Return Value:
    new("fREG",
        call = as.call(match.call()),
        formula = as.formula(formula),
        family = as.character(gaussian()),
        method = use,
        # data is as.data.frame(data), Data is as.timeSeries(data):
        data = list(data = data, Data = Data),
        fit = fit,
        residuals = fit$residuals,
        fitted = fit$fitted.values,
        title = as.character(title),
        description = as.character(description)
    )
}


###############################################################################


.lm <- 
  function(...)
{
  stats::lm(...)
}


# -----------------------------------------------------------------------------


.rlm <- 
  function(...)
{
  MASS::rlm(...)
}


# -----------------------------------------------------------------------------


.glm <- 
  function(...)
{
  stats::glm(...)
}


# -----------------------------------------------------------------------------


.gam <- 
  function(...)
{
  mgcv::gam(...)
}


# -----------------------------------------------------------------------------


.ppr <- 
  function(..., nterms = 2)
{
  stats::ppr(..., nterms = nterms)
}


# -----------------------------------------------------------------------------


.nnet <- 
  function(..., trace = FALSE, size = 2, linout = TRUE)
{
  nnet::nnet(..., trace = trace, size = size, linout = linout)
}


# -----------------------------------------------------------------------------


.polymars <- 
  function(...) 
{
  .polymarsFormula(...)
}


###############################################################################


