
## Method: 'show':         Display the object, by printing, plotting or 
##                         whatever suits its class
## Method: 'print':        Prints its argument and returns it invisibly 
##                         (via invisible(x))
## Method: 'plot':         Generic function for plotting of R objects

## Method: 'summary':      Summaries of the results of various model 
##                         fitting functions - for LaTeX
## Method: 'coef':         Regression coefficients (Ph.Eur labels)
## Method: 'anova':        Anova table (Ph.Eur labels)
## Method: 'residuals':
## Method: 'fitted':
## Method: 'vcov':
## Method: 'predict':      for prediction, including confidence and 
##                         prediction intervals
## Method: 'confint':      for confidence intervals of parameters
## Method: 'lm.influence': for regression diagnostics

setMethod("show",
    signature(object = "plaFit"),
    function (object)
    {
        X <- object@inpArgs
        pla.fit(X$data,
                sampleLabels = X$sampleLabels,
                indexOfReference = X@indexOfReference,
                StdName = X@StdName,
                design = X$design,
                dfAdj = X$dfAdj,
                dr = X$dr,
                main = X$main,
                alpha = X$alpha,
                factor = X$factor,
                show = TRUE,
                returnPotencyEstimates = FALSE)
    }
)

setMethod("print",
    signature(x = "plaFit"),
    function (x, ...)
    {
        X <- x@inpArgs
        pla.fit(X$data,
                sampleLabels = X$sampleLabels,
                indexOfReference = X@indexOfReference,
                StdName = X@StdName,
                design = X$design,
                dfAdj = X$dfAdj,
                dr = X$dr,
                main = X$main,
                alpha = X$alpha,
                factor = X$factor,
             ## Sweave = TRUE,
                show = TRUE, ...,
                returnPotencyEstimates = TRUE)
    }
)

setMethod("plot",
    signature(x = "plaFit"),
    function (x, y, ...)
    {
        X <- x@inpArgs
        cbbPalette <- c("#000000", "#009E73", "#e79f00", "#9ad0f3",
                        "#0072B2", "#D55E00", "#CC79A7", "#F0E442")
        samples <- X$sampleLabels
        colors <- cbbPalette[2:length(samples)]
        tests <- x@tests
        pla.plots(X$data,
                  sampleLabels = X$sampleLabels,
                  indexOfReference = X@indexOfReference,
                # StdName = X@StdName,
                  design = X$design,
                  main = X$main,
                  tests = x@tests, colTst = colors, ...)
    }
)

##

setMethod("summary",
    signature(object = "plaFit"),
    function (object, ...)
    {
        summary(object@lm, ...)
    }
)
setMethod("coef",
    signature(object = "plaFit"),
    function (object, ...)
    {
        coef(object@lm, ...)
    }
)
setMethod("anova",
    signature(object = "plaFit"),
    function (object, ...)
    {
        anova(object@lm, ...)
    }
)
setMethod("residuals",
    signature(object = "plaFit"),
    function (object, ...)
    {
        residuals(object@lm, ...)
    }
)
setMethod("fitted",
    signature(object = "plaFit"),
    function (object, ...)
    {
        fitted(object@lm, ...)
    }
)
setMethod("vcov",
    signature(object = "plaFit"),
    function (object, ...)
    {
        vcov(object@lm, ...)
    }
)
setMethod("predict",
    signature(object = "plaFit"),
    function (object, ...)
    {
        predict(object@lm, ...)
    }
)
setMethod("confint",
    signature(object = "plaFit"),
    function (object, parm, level = 0.95, ...)
    {
        confint(object@lm, ...)
    }
)
setMethod("lm.influence",
    signature(model = "plaFit"),
    function (model, do.coef = TRUE)
    {
        lm.influence(model@lm)
    }
)

##
