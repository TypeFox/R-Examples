
## Function: readAssayTable (readAssayResult)

## Function: assayTable2frame (table2frame)
## Function: data2assayFrame (data2frame)

## Function: pheur325 (pheur325)
## Function: pla.fit (PLAfit)
## Function: pla.plot (PLAplots)
## Function: plotSamples (plotSamples)
## Function: jitterDilutions (shiftData)

## Class: assayData

setClass(Class = "assayData",
         representation =
         representation(internals     = "list",
                        status        = "character",
                        labels        = "character",
                        model         = "character",
                        design        = "character",
                        anova         = "matrix",
                        potency       = "matrix",
                        estimates     = "matrix",
                        tableRaw      = "data.frame",
                        table         = "matrix", ## Not used, but ...
                        array         = "array",  ## Not used, but ...
                        assay         = "array",  ## Not used, but ...
                        log           = "numeric",
                        fun           = "function",
                        inv           = "function",
                        factor        = "numeric",
                        dilutionRatio = "numeric",
                        dfAdjustment  = "numeric",
                        comment       = "character",
                        description   = "character",
                        resume        = "character",
                        operator      = "character",
                        date          = "character",
                        projectTitle  = "character",
                        assayTitle    = "character"),
         prototype =
         prototype(internals     = list(),
                   status        = "",
                   labels        = "",
                   model         = "",
                   design        = "",
                   anova         = matrix(nrow=0, ncol=0),
                   potency       = matrix(nrow=0, ncol=0),
                   estimates     = matrix(nrow=0, ncol=0),
                   tableRaw      = data.frame(),
                   log           = 0,
                   fun           = function (x) x,
                   inv           = function (x) x,
                   factor        = 1,
                   dilutionRatio = 2,
                   dfAdjustment  = 0,
                   table         = matrix(nrow=0, ncol=0),
                   array         = array(dim = c(0, 0, 0)),
                   assay         = array(dim = c(0, 0)),
                   comment       = "",
                   description   = "",
                   resume        = "",
                   operator      = "",
                   date          = "",
                   projectTitle  = "",
                   assayTitle    = "")
         )

setClass(Class = "assayFrame",
         representation =
         representation(),
         prototype =
         prototype(),
         contains = "assayData")

setClass(Class = "assayTable",
         representation =
         representation(),
         prototype =
         prototype(),
         contains = "assayData")

## setClass(Class = "assayDataMarked",
##          representation = representation(mark = "character"),
##          prototype = prototype(mark = ""),
##          contains = "assayData")

## Class: assayModel

setClass(Class = "assayModel",
         representation =
         representation(data = "data.frame",
                        alpha = "numeric",
                        dfAdjustment = "numeric",
                        dilutionRatio = "numeric",
                        selectFun = "function",
                        factor = "numeric",
                        sampleLabels = "character",
                        indexOfReference = "numeric",
                        StdName = "character",
                        colors = "character",
                        projectTitle = "character",
                        assayTitle = "character",
                        imputeMissing = "logical",
                        model = "character",
                        design = "character"),
         prototype =
         prototype(data = data.frame(),
                   alpha = 0.05,
                   dfAdjustment = 0,
                   dilutionRatio = 2,
                   selectFun = function (array) NULL,
                   factor = 1.0,
                   projectTitle = "",
                   assayTitle   = "",
                   sampleLabels = c("S", "T"),
                   indexOfReference = 1,
                   StdName = "S",
                   colors = c("black", "grey"),
                   imputeMissing = FALSE,
                   model = "",
                   design = ""))

## Class: fpl / ssa - S-share (Sigmoid):

## Class: sra - Slope Ratio:

## Class: pla - Parallel Line:

setClass(Class = "pla",
         representation =
         representation(),
         prototype =
         prototype(),
         contains = "assayModel")

setClass(Class = "plaCRD",
         representation =
         representation(),
         prototype =
         prototype(),
         contains = "pla")

setClass(Class = "plaLSD",
         representation =
         representation(),
         prototype =
         prototype(),
         contains = "pla")

setClass(Class = "plaRBD",
         representation =
         representation(),
         prototype =
         prototype(),
         contains = "pla")

## Class: assayFitUlisted:

## setClass(Class = "assayFitUnlisted",
##          representation =
##          representation(lm         = "lm",
##                         pheur      = "list",
##                         tests      = "list",
##                         anova      = "data.frame",
##                       # AT         = "data.frame",
##                       # SS         = "numeric",
##                       # reg        = "numeric",
##                       # K          = "numeric",
##                       # KP         = "numeric",
##                         relAnova   = "data.frame",
##                         relPotency = "data.frame"),
##          prototype =
##          prototype(lm         = list(),
##                    pheur      = list(),
##                    tests      = NULL,
##                    anova      = data.frame(),
##                  # AT         = NULL,
##                  # SS         = 0,
##                  # reg        = 0,
##                  # K          = matrix(nrow=0, ncol=0),
##                  # KP         = matrix(nrow=0, ncol=0),
##                    relAnova   = data.frame(),
##                    relPotency = data.frame()
##                    ),
##          contains = "assayModel")

## Class: plaFit:

setClass(Class = "plaFit",
         representation =
         representation(inpArgs    = "list",
                        model      = "character",
                        design     = "character",
                        lm         = "lm",
                        pheur      = "list",
                        tests      = "list",
                        anova      = "data.frame",
                      # AT         = "data.frame",
                      # SS         = "numeric",
                      # reg        = "numeric",
                      # K          = "numeric",
                      # KP         = "numeric",
                        relAnova   = "matrix",
                        relPotency = "matrix"
                      # relAnova   = "data.frame",
                      # relPotency = "data.frame"
                        ),
         prototype =
         prototype(model      = "",
                   design     = "",
                   lm         = list(),
                   pheur      = list(),
                   anova      = matrix(nrow=0, ncol=0),
                   tests      = list(),
                   anova      = data.frame(),
                 # AT         = NULL,
                 # SS         = 0,
                 # reg        = 0,
                 # K          = matrix(nrow=0, ncol=0),
                 # KP         = matrix(nrow=0, ncol=0),
                   relAnova   = matrix(nrow=0, ncol=0),
                   relPotency = matrix(nrow=0, ncol=0)
                 # relAnova   = data.frame(),
                 # relPotency = data.frame()
                   ),
         )
