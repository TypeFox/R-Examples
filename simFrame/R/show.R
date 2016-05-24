# ----------------------
# Author: Andreas Alfons
#         KU Leuven
# ----------------------

## data control
setMethod("show", "DataControl", 
    function(object) {
        cat("Number of observations to be generated:\n")
        print(getSize(object))
#        cat("\nGenerating (distribution) function:\n")
#        print(getDistribution(object))
        dots <- getDots(object)
        if(length(dots) > 0) {
            cat("\nAdditional parameters:\n")
            print(dots)
        }
        colnames <- getColnames(object)
        if(!is.null(colnames)) {
            cat("\nColumn names of the resulting data frame:\n")
            print(colnames)
        }
    })


## sample control
setMethod("show", "VirtualSampleControl", 
    function(object) {
        cat("Number of samples to be set up:\n")
        print(getK(object))
    })

# single-stage sampling
setMethod("show", "SampleControl", 
    function(object) {
        callNextMethod()
        design <- getDesign(object)
        haveDesign <- length(design) > 0
        if(haveDesign) {
            cat("\nDesign variable(s) for stratified sampling:\n")
            print(design)
        }
        group <- getGrouping(object)
        if(length(group) > 0) {
            cat("\nGrouping variable giving sampling units:\n")
            print(group)
        }
        size <- getSize(object)
        if(!is.null(size)) {
            tmp <- "\nSize of the resulting samples"
            if(haveDesign) {
                if(length(size) == 1) tmp <- paste(tmp, "for each stratum")
                else tmp <- paste(tmp, "by stratum")
            }
            cat(paste(tmp, ":\n", sep=""))
            print(size)
        }
    })

# two-stage sampling
setMethod("show", "TwoStageControl", 
    function(object) {
        callNextMethod()
        design <- getDesign(object)
        haveDesign <- length(design) > 0
        if(haveDesign) {
            cat("\nDesign variable(s) for stratified sampling:\n")
            print(design)
        }
        group <- getGrouping(object)
        haveSSUs <- length(group) > 1
        if(haveSSUs) {
            cat("\nGrouping variables giving PSUs and SSUs, respectively:\n")
        } else cat("\nGrouping variable giving PSUs:\n")
        print(group)
        size <- getSize(object, stage=1)
        if(!is.null(size)) {
            tmp <- "\nNumber of PSUs to be sampled"
            if(haveDesign) {
                if(length(size) == 1) tmp <- paste(tmp, "for each stratum")
                else tmp <- paste(tmp, "by stratum")
            }
            cat(paste(tmp, ":\n", sep=""))
            print(size)
        }
        size <- getSize(object, stage=2)
        if(!is.null(size)) {
            tmp <- if(haveSSUs) "SSUs" else "individuals"
            tmp <- paste("\nNumber of", tmp, "to be sampled")
            tmp <- paste(tmp, if(length(size) == 1) "for each " else "by ")
            cat(paste(tmp, "PSU:\n", sep=""))
            print(size)
        }
    })

## sample setup
setMethod("show", "SampleSetup", 
    function(object) {
        k <- length(object)
        if(k == 0) cat("No samples are set up\n")
        else {
            if(k == 1) cat("Indices of observations in the sample:\n\n")
            else {
                cat("Indices of observations for each of the", k, "samples:\n\n")
            }
            print(getIndices(object))
        }
    })

setMethod("show", "SummarySampleSetup", 
    function(object) {
        size <- getSize(object)
        k <- length(size)
        if(k == 0) cat("No samples are set up\n")
        else {
            if(length(unique(size)) == 1) {  # also TRUE if k == 1
                if(k == 1) msg <- "%i sample of size %i is set up\n"
                else msg <- "%i samples of size %i are set up\n"
                cat(sprintf(msg, k, size[1]))
            } else {
                cat(k, "samples are set up\n\n")
                cat("The sample sizes are:\n")
                print(size)
            }
        }
    })

## contamination control
setMethod("show", "VirtualContControl", 
    function(object) {
        target <- getTarget(object)
        if(is.null(target)) cat("All variables are target variables\n")
        else {
            if(length(target) == 1) cat("Target variable:\n")
            else cat("Target variables:\n")
            print(target)
        }
        cat("\nEpsilon:\n")
        print(getEpsilon(object))
    })

setMethod("show", "ContControl", 
    function(object) {
        callNextMethod()
        group <- getGrouping(object)
        if(length(group) > 0) {
            cat("\nGrouping variable for contaminating clusters:\n")
            print(group)
        }
        aux <- getAux(object)
        if(length(aux) > 0) {
            cat("\nVariable giving probability weights for selection:\n")
            print(aux)
        }
    })

## NA control
setMethod("show", "VirtualNAControl", 
    function(object) {
        target <- getTarget(object)
        if(is.null(target)) cat("All variables are target variables\n")
        else {
            if(length(target) == 1) cat("Target variable:\n")
            else cat("Target variables:\n")
            print(target)
        }
        cat("\nMissing value rates:\n")
        print(getNArate(object))
    })

setMethod("show", "NAControl", 
    function(object) {
        callNextMethod()
        group <- getGrouping(object)
        if(length(group) > 0) {
            cat("\nGrouping variable for setting clusters to NA:\n")
            print(group)
        }
        aux <- getAux(object)
        if(length(aux) > 0) {
            if(length(aux) == 1) {
                cat("\nVariable giving probability weights for selection:\n")
            } else {
                cat("\nVariables giving probability weights for selection:\n")
            }
            print(aux)
        }
    })

## strata information
setMethod("show", "Strata", 
    function(object) {
        cat("Strata:\n")
        print(getValues(object))
#        cat("\nLegend:\n")
#        print(getLegend(object))
        cat("\nStratum sizes:\n")
        show(summary(object))
    })

## simulation control
setMethod("show", "SimControl", 
    function(object) {
        contControl <- getContControl(object)
        if(is.null(contControl)) {
            cat("Contamination will not be added\n")
        } else {
            cat("Settings for contamination:\n")
            show(contControl)
        }
        NAControl <- getNAControl(object)
        blanks <- if(is.null(contControl) && is.null(NAControl)) "\n" else "\n\n"
        cat(blanks)
        if(is.null(NAControl)) {
            cat("NAs will not be inserted\n")
            blanks <- "\n"
        } else {
            cat("Settings for inserting NAs:\n")
            show(NAControl)
            blanks <- "\n\n"
        }
        design <- getDesign(object)
        if(length(design) > 0) {
            cat(blanks)
            cat("Variable(s) giving domains for splitting simulations:\n")
            print(design)
            blanks <- "\n"
        }
        SAE <- getSAE(object)
        if(isTRUE(SAE)) {
            cat(blanks)
            cat("Small area estimation is expected to be used\n")
        }
    })

## simulation results
setMethod("show", "SimResults", function(object) print(getValues(object)))
