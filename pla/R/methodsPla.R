
## Method: 'show':      Display the object, by printing, plotting or 
##                      whatever suits its class
## Method: 'print':     Prints its argument and returns it invisibly 
##                      (via invisible(x))
## Method: 'plot':      Generic function for plotting of R objects

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
    signature(object = "pla"),
    function (object)
    {
        pla.fit(object@data,
                sampleLabels = object@sampleLabels,
                indexOfReference = object@indexOfReference,
                StdName = object@StdName,
                design = object@design,
                dfAdj = object@dfAdjustment,
                dr = object@dilutionRatio,
                main = object@assayTitle,
                alpha = object@alpha,
                factor = object@factor,
                show = TRUE,
                returnPotencyEstimates = FALSE)
    }
)

setMethod("print",
    signature(x = "pla"),
    function (x, ...)
    {
        if (any(names(list(...)) == "show"))
            pla.fit(x@data,
                    sampleLabels = x@sampleLabels,
                    indexOfReference = x@indexOfReference,
                    StdName = x@StdName,
                    design = x@design,
                    dfAdj = x@dfAdjustment,
                    dr = x@dilutionRatio,
                    main = x@assayTitle,
                    alpha = x@alpha,
                    factor = x@factor,
                    ## Sweave = TRUE,
                    ...,
                    returnPotencyEstimates = TRUE)
        else
            pla.fit(x@data,
                    sampleLabels = x@sampleLabels,
                    indexOfReference = x@indexOfReference,
                    StdName = x@StdName,
                    design = x@design,
                    dfAdj = x@dfAdjustment,
                    dr = x@dilutionRatio,
                    main = x@assayTitle,
                    alpha = x@alpha,
                    factor = x@factor,
                    ## Sweave = TRUE,
                    show = TRUE, ...,
                    returnPotencyEstimates = TRUE)
    }
)

setMethod("plot",
    signature(x = "pla"),
    function (x, y, ...)
    {
        pla.plots(x@data,
                  sampleLabels = x@sampleLabels,
                  indexOfReference = x@indexOfReference,
                # StdName = x@StdName,
                  design = x@design,
                  main = x@assayTitle,
                  colRef = x@colors[1],
                  colTst = x@colors[-1], ...)
    }
)

setMethod("as.data.frame",
    signature(x = "pla"),
    function (x, ...)
    {
        return(x@data)
    }
)

setMethod("summary",
    signature(object = "pla"),
    function (object, ...)
    {
        fit <- pla.fit(object@data,
                       sampleLabels = object@sampleLabels,
                       indexOfReference = object@indexOfReference,
                       StdName = object@StdName,
                       design = object@design,
                       dfAdj = object@dfAdjustment,
                       dr = object@dilutionRatio,
                       main = object@assayTitle,
                       alpha = object@alpha,
                       factor = object@factor,
                       show = FALSE, ...,
                       returnPotencyEstimates = TRUE)
        summary(fit@lm, ...)

    }
)
setMethod("coef",
    signature(object = "pla"),
    function (object, ...)
    {
        fit <- pla.fit(object@data,
                       sampleLabels = object@sampleLabels,
                       indexOfReference = object@indexOfReference,
                       StdName = object@StdName,
                       design = object@design,
                       dfAdj = object@dfAdjustment,
                       dr = object@dilutionRatio,
                       main = object@assayTitle,
                       alpha = object@alpha,
                       factor = object@factor,
                       show = FALSE, ...,
                       returnPotencyEstimates = TRUE)
        coef(fit@lm, ...)
    }
)
setMethod("anova",
    signature(object = "pla"),
    function (object, ...)
    {
        fit <- pla.fit(object@data,
                       sampleLabels = object@sampleLabels,
                       indexOfReference = object@indexOfReference,
                       StdName = object@StdName,
                       design = object@design,
                       dfAdj = object@dfAdjustment,
                       dr = object@dilutionRatio,
                       main = object@assayTitle,
                       alpha = object@alpha,
                       factor = object@factor,
                       show = FALSE, ...,
                       returnPotencyEstimates = TRUE)
        anova(fit@lm, ...)
    }
)
setMethod("residuals",
    signature(object = "pla"),
    function (object, ...)
    {
        fit <- pla.fit(object@data,
                       sampleLabels = object@sampleLabels,
                       indexOfReference = object@indexOfReference,
                       StdName = object@StdName,
                       design = object@design,
                       dfAdj = object@dfAdjustment,
                       dr = object@dilutionRatio,
                       main = object@assayTitle,
                       alpha = object@alpha,
                       factor = object@factor,
                       show = FALSE, ...,
                       returnPotencyEstimates = TRUE)
        residuals(fit@lm, ...)
    }
)
setMethod("fitted",
    signature(object = "pla"),
    function (object, ...)
    {
        ## Not working!
        fit <- pla.fit(object@data,
                       sampleLabels = object@sampleLabels,
                       indexOfReference = object@indexOfReference,
                       StdName = object@StdName,
                       design = object@design,
                       dfAdj = object@dfAdjustment,
                       dr = object@dilutionRatio,
                       main = object@assayTitle,
                       alpha = object@alpha,
                       factor = object@factor,
                       show = FALSE, ...,
                       returnPotencyEstimates = TRUE)
        fitted(fit@lm, ...)
    }
)
setMethod("vcov",
    signature(object = "pla"),
    function (object, ...)
    {
        fit <- pla.fit(object@data,
                       sampleLabels = object@sampleLabels,
                       design = object@design,
                       dfAdj = object@dfAdjustment,
                       dr = object@dilutionRatio,
                       main = object@assayTitle,
                       alpha = object@alpha,
                       factor = object@factor,
                       show = FALSE, ...,
                       returnPotencyEstimates = TRUE)
        vcov(fit@lm, ...)
    }
)
setMethod("predict",
    signature(object = "pla"),
    function (object, ...)
    {
        fit <- pla.fit(object@data,
                       sampleLabels = object@sampleLabels,
                       indexOfReference = object@indexOfReference,
                       StdName = object@StdName,
                       design = object@design,
                       dfAdj = object@dfAdjustment,
                       dr = object@dilutionRatio,
                       main = object@assayTitle,
                       alpha = object@alpha,
                       factor = object@factor,
                       show = FALSE, ...,
                       returnPotencyEstimates = TRUE)
        predict(fit@lm, ...)
    }
)
setMethod("confint",
    signature(object = "pla"),
    function (object, parm, level = 0.95, ...)
    {
        fit <- pla.fit(object@data,
                       sampleLabels = object@sampleLabels,
                       indexOfReference = object@indexOfReference,
                       StdName = object@StdName,
                       design = object@design,
                       dfAdj = object@dfAdjustment,
                       dr = object@dilutionRatio,
                       main = object@assayTitle,
                       alpha = object@alpha,
                       factor = object@factor,
                       show = FALSE, ...,
                       returnPotencyEstimates = TRUE)
        confint(fit@lm, ...)
    }
)
setMethod("lm.influence",
    signature(model = "pla"),
    function (model, do.coef = TRUE)
    {
        fit <- pla.fit(model@data,
                       sampleLabels = model@sampleLabels,
                       indexOfReference = model@indexOfReference,
                       StdName = model@StdName,
                       design = model@model,
                       dfAdj = model@dfAdjustment,
                       dr = model@dilutionRatio,
                       main = model@assayTitle,
                       alpha = model@alpha,
                       factor = model@factor,
                       show = FALSE,
                       returnPotencyEstimates = TRUE)
        lm.influence(fit@lm)
    }
)
