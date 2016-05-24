
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


################################################################################
# FUNCTION:             REGRESSION MODELLING DESCRIPTION:
#  regSim                Returns a regression example data set
#  regFit.dataframe
#  regFit.valueSlots
#  predict.fREG          Predicts values from a fitted regression model
#  regFit.nonDefaults
#  generalizedModels
################################################################################


test.regSim <-
    function()
{
    # Plot Parameters:
    par(ask = FALSE)
    par(mfrow = c(3, 1))

    # Simulate Artificial LM:
    X = regSim(model = "LM3", n = 365)
    head(X)
    plot(X[, "Y"], type = "l", main = "LM3", xlab = "1970", ylab = "Y")

    # Simulate Artificial LOGIT:
    X = regSim(model = "LOGIT3", n = 365)
    head(X)
    plot(X[, "Y"], type = "l", main = "LOGIT3", xlab = "1970", ylab = "Y")

    # Simulate Artificial GAM:
    X = regSim(model = "GAM3", n = 365)
    head(X)
    plot(X[, "Y"], type = "l", main = "GAM3", xlab = "1970", ylab = "Y")

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.regFit.dataframe <-
    function()
{
    # Working with timeSeries Objects ...
    DATA = regSim(model = "GAM3", n = 100)
    head(DATA)
    class(DATA)

    # Regression Fit:
    LM       = regFit(Y ~ X1 + X2, data = DATA, use = "lm")
    RLM      = regFit(Y ~ X1 + X2, data = DATA, use = "rlm")
    AM       = regFit(Y ~ X1 + X2, data = DATA, use = "gam")
    PPR      = regFit(Y ~ X1 + X2, data = DATA, use = "ppr")
    POLYMARS = regFit(Y ~ X1 + X2, data = DATA, use = "polymars")
    ## NNET     = regFit(Y ~ X1 + X2, data = DATA, use = "nnet")
    # ... a note on AM the smoothing functions are added by default!
    # this is different to gam()


    # Print Method:
    print(LM)
    print(RLM)
    print(AM)
    print(PPR)
    print(POLYMARS)
    ## print(NNET)

    # Plot Method:
    par(ask = FALSE)
    par(mfrow = c(1, 1))
    # plot(LM, which = "all")                             # CHECK which !!!
    # plot(RLM, which = "all")
    # plot(AM, which = "all")
    # plot(PPR, which = "all")
    # plot(POLYMARS, which = "all")
    # plot(NNET, which = "all")

    # Summary Method:
    summary(LM)
    summary(RLM)
    summary(AM)
    summary(PPR)
    summary(POLYMARS)
    ## summary(NNET)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.regFit.valueSlots <-
    function()
{

    # Working with timeSeries Objects ...
    DATA = regSim(model = "GAM3", n = 100)
    head(DATA)
    class(DATA)
    
    require(mgcv)

    # Modelling:
    LM    = regFit(Y ~ X1 + X2, data = DATA, use = "lm")
    RLM   = regFit(Y ~ X1 + X2, data = DATA, use = "rlm")
    AM    = regFit(Y ~ s(X1) + s(X2),  DATA, use = "gam")
    PPR   = regFit(Y ~ X1 + X2, data = DATA, use = "ppr")
    POLYMARS = regFit(Y ~ X1 + X2, data = DATA, use = "polymars")
    NNET  = regFit(Y ~ X1 + X2, data = DATA, use = "nnet")

    # Extract:

    # call = "call"
    # formula = "formula"
    # family = "character"
    # method = "character"
    # data = "data.frame"
    # fit = "list"
    # residuals = "timeSeries"
    # fitted.values = "timeSeries"
    # title = "character"
    # description = "character"

    LM@call
    RLM@call
    AM@call
    PPR@call
    POLYMARS@call
    NNET@call

    LM@formula
    RLM@formula
    AM@formula                                                       # CHECK !!!
    PPR@formula
    POLYMARS@formula
    NNET@formula

    LM@family[1:2]
    RLM@family[1:2]
    AM@family[1:2]
    PPR@family[1:2]
    POLYMARS@family[1:2]
    NNET@family[1:2]

    LM@method
    RLM@method
    AM@method
    PPR@method
    POLYMARS@method
    NNET@method

    # Note the residuals are time tmeSeries objects!
    print(LM@residuals[c(1,100)])
    print(RLM@residuals[c(1,100)])
    print(AM@residuals[c(1,100)])
    print(PPR@residuals[c(1,100)])
    print(POLYMARS@residuals[c(1,100)])
    print(NNET@residuals[c(1,100)])

    # Note the fitted values are time tmeSeries objects!
    print(LM@fitted[c(1,100)])
    print(RLM@fitted[c(1,100)])
    print(AM@fitted[c(1,100)])
    print(PPR@fitted[c(1,100)])
    print(POLYMARS@fitted[c(1,100)])
    print(NNET@fitted[c(1,100)])

    # Returns a Title, by default the name of the algorithm applied:
    LM@title
    RLM@title
    AM@title
    PPR@title
    POLYMARS@title
    NNET@title

    # Returns a Description, by default Date/Time and user:
    LM@description
    RLM@description
    AM@description
    PPR@description
    POLYMARS@description
    NNET@description

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.predict.fREG <-
    function()
{
    # Working with timeSeries Objects ...
    DATA <- regSim(model = "GAM3", n = 100)
    head(DATA)
    class(DATA)

   require(mgcv)
    
    # Regression Fit:
    LM    = regFit(Y ~ X1 + X2, data = DATA, use = "lm")
    RLM   = regFit(Y ~ X1 + X2, data = DATA, use = "rlm")
    AM    = regFit(Y ~ s(X1) + s(X2),  DATA, use = "gam")
    PPR   = regFit(Y ~ X1 + X2, data = DATA, use = "ppr")
    POLYMARS = regFit(Y ~ X1 + X2, data = DATA, use = "polymars")
    NNET  = regFit(Y ~ X1 + X2, data = DATA, use = "nnet")

    # Just to rmember - Predict:
    #   predict.fREG(object, newdata, se.fit = FALSE, type = "response", ...)

    # Selext some rows to predict:
    set.seed(4711)
    N <- round(runif(5, 1, 100), 0)

    # Predict Response:
    predict(LM,    DATA[N, ])
    predict(RLM,   DATA[N, ])
    predict(AM,    DATA[N, ])
    predict(PPR,   DATA[N, ])
    predict(POLYMARS, DATA[N, ])
    ## predict(NNET,  DATA[N, ])

    # Predict Response:
    predict(LM,    DATA[N, ], type = "response")
    predict(RLM,   DATA[N, ], type = "response")
    predict(AM,    DATA[N, ], type = "response")
    predict(PPR,   DATA[N, ], type = "response")
    predict(POLYMARS, DATA[N, ], type = "response")
    ## predict(NNET,  DATA[N, ], type = "response")

    # Predict Response with Standard Errors:
    predict(LM,    DATA[N, ], se.fit = TRUE)
    predict(RLM,   DATA[N, ], se.fit = TRUE)
    predict(AM,    DATA[N, ], se.fit = TRUE)
    predict(PPR,   DATA[N, ], se.fit = TRUE)
    predict(POLYMARS, DATA[N, ], se.fit = TRUE)
    ## predict(NNET,  DATA[N, ], se.fit = TRUE)

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.regFit.nonDefaults <-
    function()
{

    # Simulate Data - a data frame:
    DATA = regSim(model = "GAM3", n = 100)
    head(DATA)
    class(DATA)

    # Simulate Data - a timeSeries object:
    DATA = as.timeSeries(DATA)
    head(DATA)
    class(DATA)

    # LM:
    LM1 = regFit(Y ~ X1 + X2,      DATA, use = "lm")
    print(LM1)
    LM2 = regFit(Y ~ 1 + X1 + X2,  DATA)
    print(LM2)
    LM3 = regFit(Y ~ -1 + X1 + X2, DATA)
    print(LM3)
    LM4 = regFit(Y ~ X1 + log(X2), DATA)
    print(LM4)

    require(mgcv)
    
    # AM:
    AM1 = regFit(Y ~ s(X1) + s(X2), data = DATA, use = "gam")
    print(AM1)
    # AM2 = regFit(Y ~ s(X1) + s(X2), DATA, "gam",
    #   method = gam.method(pearson = TRUE))
    # print(AM2)

    # PPR:
    par(ask = FALSE)
    par(mfrow = c(1, 1))
    PPR1 = regFit(Y ~ sin(X1) + exp(X2), DATA, "ppr", nterms = 4,
        sm.method = "supsmu", use = "ppr")
    PPR2 = regFit(Y ~ sin(X1) + exp(X2), DATA, "ppr", nterms = 4,
        sm.method = "spline", use = "ppr")
    PPR3 = regFit(Y ~ sin(X1) + exp(X2), DATA, "ppr", nterms = 3,
        sm.method = "gcvspline", use = "ppr")
    ## termPlot(PPR1)
    ## termPlot(PPR2)
    ## termPlot(PPR3)

    # POLYMARS:
    POLYMARS <- regFit(Y ~ X1 + X2 + X3, DATA, use = "polymars")
    POLYMARS <- regFit(Y ~ X1*X2 + X2*X3 + X3*X1, DATA, use = "polymars")

    # NNET
    # todo ...

    # Return Value:
    return()
}


# ------------------------------------------------------------------------------


test.generalizedModels <-
    function()
{
    # Generalized * Models:

    M1 <- matrix(c(
       1, 0.80, 0.83, 0.66, 1.9, 1.100, 0.996,
       1, 0.90, 0.36, 0.32, 1.4, 0.740, 0.992,
       0, 0.80, 0.88, 0.70, 0.8, 0.176, 0.982,
       0, 1.00, 0.87, 0.87, 0.7, 1.053, 0.986,
       1, 0.90, 0.75, 0.68, 1.3, 0.519, 0.980,
       0, 1.00, 0.65, 0.65, 0.6, 0.519, 0.982,
       1, 0.95, 0.97, 0.92, 1.0, 1.230, 0.992,
       0, 0.95, 0.87, 0.83, 1.9, 1.354, 1.020,
       0, 1.00, 0.45, 0.45, 0.8, 0.322, 0.999,
       0, 0.95, 0.36, 0.34, 0.5, 0.000, 1.038,
       0, 0.85, 0.39, 0.33, 0.7, 0.279, 0.988,
       0, 0.70, 0.76, 0.53, 1.2, 0.146, 0.982,
       0, 0.80, 0.46, 0.37, 0.4, 0.380, 1.006,
       0, 0.20, 0.39, 0.08, 0.8, 0.114, 0.990,
       0, 1.00, 0.90, 0.90, 1.1, 1.037, 0.990,
       1, 1.00, 0.84, 0.84, 1.9, 2.064, 1.020,
       0, 0.65, 0.42, 0.27, 0.5, 0.114, 1.014,
       0, 1.00, 0.75, 0.75, 1.0, 1.322, 1.004,
       0, 0.50, 0.44, 0.22, 0.6, 0.114, 0.990,
       1, 1.00, 0.63, 0.63, 1.1, 1.072, 0.986,
       0, 1.00, 0.33, 0.33, 0.4, 0.176, 1.010,
       0, 0.90, 0.93, 0.84, 0.6, 1.591, 1.020,
       1, 1.00, 0.58, 0.58, 1.0, 0.531, 1.002,
       0, 0.95, 0.32, 0.30, 1.6, 0.886, 0.988,
       1, 1.00, 0.60, 0.60, 1.7, 0.964, 0.990,
       1, 1.00, 0.69, 0.69, 0.9, 0.398, 0.986,
       0, 1.00, 0.73, 0.73, 0.7, 0.398, 0.986),
       byrow = TRUE, ncol = 7)
    colnames(M1) = c("Y", "X1", "X2", "X3", "X4", "X5", "X6")
    D1 = data.frame(M1)
    D1

    # fit.glm = glm(Y ~ X1 + X2 + X3 + X4 + X5 + X6, data = D1,
    #     family = binomial("logit"))
    # fit.gam = gam(Y ~ s(X1) + s(X2) + s(X3) + s(X4) + s(X5) + s(X6),
    #     data = D1, family = binomial("logit"))


    M2 <- matrix(c(
        0,29,62,
        0,30,83,
        0,31,74,
        0,31,88,
        0,32,68,
        1,29,41,
        1,30,44,
        1,31,21,
        1,32,50,
        1,33,33),
       byrow = TRUE, ncol = 3)

    colnames(M2) = c("Y", "X1", "X2")
    D2 = data.frame(M2)
    D2

    plot  (D2[1:5, "X1"], D2[1:5, "X2"],
        xlim = range(D2[, "X1"]), ylim = range(D2[, "X2"]),
        pch = 19, col = "blue")
    points(D2[6:10,"X1"], D2[6:10,"X2"],
        pch = 19, col = "red")
    U = range(D2[, "X1"])
    V = 2*U - 6
    lines(U, V, lty = 3, col = "grey")

    fit.glm = glm(Y ~ X1 + X2, data = D2, family = binomial("logit"))
    print(fit.glm)

    # Return Value:
    return()
}


################################################################################
