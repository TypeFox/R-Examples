#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* campus *.* lmu *.* de]
## Time-stamp: <[extractSamples.R] by DSB Die 30/03/2010 10:09 (CEST)>
##
## Description:
## Extract samples and prediction info from BayesX results directory.
##
## History:
## 05/02/2009   file creation
## 06/02/2009   add na.string "." to import of predictions table
## 10/02/2009   add import of spatial estimates
## 23/02/2009   plug into BayesX package,
##              see SVN logs for the further history! 
#####################################################################################


## returns list of MCMC objects for all samples in BayesX results directory
extractSamples <- function(directoryWithBasename,
                           logfile= # log with the prior specifications, necessary for pspline
                                        # and MCMC option extraction
                           file.path(dirname(directoryWithBasename),
                                     "log.txt"))
{
    ## extracts
    resBasename <- basename(directoryWithBasename)
    resDirname <- dirname(directoryWithBasename)

    ## upgrade to full absolute path
    directoryWithBasename <- file.path(resDirname, resBasename)
    
    ## which files are to be processed?
    resFiles <- list.files(path=resDirname,
                           pattern=
                           paste(resBasename,
                                 "_.+_sample\\.raw",
                                 sep=""),
                           full.names=TRUE)

    ## Extract MCMC parameters from log file
    bayesxLog <- readLines(logfile)

    ## use internal helper function, defined in separate R code file
    ## to extract the numbers
    numbers <- getNumbers(beforeStringsList=
                          list(Iterations="Number of iterations:",
                               BurnIn="Burn-in period:",
                               Thin="Thinning parameter:"),
                          stringVector=bayesxLog)

    ## setup return list
    ret <- list()

    ## convenience functions:
    convert2Mcmc <- function(samples)   # just a shortcut
        coda::mcmc(data=samples,
                   start=numbers$BurnIn + 1,
                   end=numbers$Iterations,
                   thin=numbers$Thin)
    
    readData <- function(file)          # just a shortcut
        read.table(file,
                   header=TRUE,
                   row.names=1)

    getResFile <- function(sampleFile)  # just a shortcut
        sub(pattern="(.+)_sample\\.raw",
            replacement="\\1\\.res",
            x=sampleFile)
    
    readNamedSamples <-
        function(sampleFile, # .raw file
                 resFile=getResFile(sampleFile)) # corresponding .res file containing the varnames 
        {
            ## read unnamed data
            sampleData <- readData(sampleFile)
            ## and label correctly, if possible
            sampleNames <- readData(resFile)[, 1] # assumes that names are in first column!
            if(identical(ncol(sampleData),
                         length(sampleNames)))
            {
                colnames(sampleData) <- sampleNames
            }
            else
            {
                warning("For sample file ", sampleFile, " the names read from ", resFile,
                        " did not match the number of columns / parameters!")
            }
            
            return(sampleData)
        }

    ## process fixed effects
    fixedInds <- grep(pattern=".+FixedEffects[[:digit:]]+_sample\\.raw",
                      x=resFiles)
    if(length(fixedInds))
    {
        samplesMatrix <- matrix(nrow=
                                with(numbers,
                                     (Iterations - BurnIn) / Thin),
                                ncol=0)
        
        for(sampleFile in resFiles[fixedInds])
        {
            samplesMatrix <- cbind(samplesMatrix,
                                   as.matrix(readNamedSamples(sampleFile)))
        }
        ret$FixedEffects <- convert2Mcmc(samplesMatrix)
    }

    ## process random effects
    randomInds <- grep(pattern=".+_random_sample\\.raw",
                       x=resFiles)

    ## start random effects list if there is at least one such sample file
    if(length(randomInds))
        ret$RandomEffects <- list()  
    for(sampleFile in resFiles[randomInds])
    {
        ## read function samples from this file, first without names,
        ## because we cannot be sure that the res file has not too few name
        ## because there could be a fixed effect included (see below)
        functionSamples <- readData(sampleFile)

        ## and the variance samples from the corresponding file
        varianceSamples <- readData(sub(pattern="(.+)_sample\\.raw",
                                        replacement="\\1_variance_sample\\.raw",
                                        x=sampleFile))[, 1]
        
        filenameParts <- strsplit(x=basename(sampleFile),
                                  split="_f_",
                                  fixed=TRUE)[[1]]
       
        ## the first part before "_f_" and after the basename
        ## is the name of the covariate for which a random effect was specified
        covName <-
            if(identical(filenameParts[1],
                         resBasename))
                "const"
            else
                sub(pattern=
                    paste(resBasename, "_",
                          sep=""),
                    replacement="",
                    x=filenameParts[1])

        ## the name of the group id is before the tail of the second filename part 
        idName <- sub(pattern="_random_sample\\.raw",
                      replacement="",
                      x=filenameParts[2])

        ## check if fixed effects are included in random effects sample matrix,
        ## and move them to the fixed effects. This is e.g. the case when 
        ## "time * id(random)" is included in the formula,
        ## and not "time + time * id(random, nofixed)".

        if(file.exists(fixedPartResFile <- sub(pattern="sample.raw",
                                               replacement="fixed.res",
                                               x=sampleFile,
                                               fixed=TRUE)))
        {
            ## move this last column to the fixed effects:

            ## the new fixed effects mcmc matrix
            oldColnames <- colnames(ret$FixedEffects)
            ret$FixedEffects <- convert2Mcmc(cbind(ret$FixedEffects,
                                                   functionSamples[, ncol(functionSamples)]))
            colnames(ret$FixedEffects) <- c(oldColnames, covName)

            ## delete the column from the random effects samples
            functionSamples <- functionSamples[, - ncol(functionSamples)]
        }

        ## assign proper column names (the ID strings) to the random effect samples
        sampleNames <- readData(getResFile(sampleFile))[, 1] 
        colnames(functionSamples) <- sampleNames
        
        ## now insert the mcmc objects converted samples into the hierarchy
        theseSamples <- list(list(functionSamples=convert2Mcmc(functionSamples),
                                  varianceSamples=convert2Mcmc(varianceSamples)))
        names(theseSamples) <- covName
        
        ret$RandomEffects[[idName]] <- c(ret$RandomEffects[[idName]],
                                         theseSamples)
    }

   
    ## process nonlinear functions with rw or spatial priors
    rwInds <- grep(pattern=".+_(rw|spatial)_sample\\.raw",
                   x=resFiles)  
    for(sampleFile in resFiles[rwInds])
    {
        ## read function samples from this file
        functionSamples <- readNamedSamples(sampleFile)

        ## and the variance samples from the corresponding file
        varianceSamples <- readData(sub(pattern="(.+)_sample\\.raw",
                                        replacement="\\1_variance_sample\\.raw",
                                        x=sampleFile))[, 1]
        
        ## coerce to MCMC objects and insert into list with correct name
        functionName <- sub(pattern=
                            paste(directoryWithBasename,
                                  "(.+)_(rw|spatial)_sample\\.raw",
                                  sep="_"),
                            replacement="\\1",
                            x=sampleFile)
        ret[[functionName]] <- list(functionSamples=convert2Mcmc(functionSamples),
                                    varianceSamples=convert2Mcmc(varianceSamples))
    }   
    
    ## process nonlinear functions modelled as psplines
    psplineInds <- grep(pattern=".+_pspline_sample\\.raw",
                        x=resFiles)  
    for(sampleFile in resFiles[psplineInds])
    {       

        ## get corresponding covariate values (or gridpoints if it was restricted)
        covValues <- readData(getResFile(sampleFile))[, 1]

        ## get name of function
        functionName <- sub(pattern=
                            paste(directoryWithBasename,
                                  "(.+)_pspline_sample\\.raw",
                                  sep="_"),
                            replacement="\\1",
                            x=sampleFile)
        
        ## extract pspline parameters from log file
        optionsPart <- bayesxLog[grep(pattern=
                                      paste("[[:blank:]]*OPTIONS FOR P-SPLINE TERM:",
                                            functionName),
                                      x=bayesxLog)
                                 + (1:10)] # look for numbers in the next ten lines
        
        optionsNumbers <- getNumbers(beforeStringsList=
                                     list(knots="Number of knots:",
                                          degree="Degree of Splines:"),
                                     stringVector=optionsPart)

        ## derive knot locations
        eps <- 0.001
        minx <- min(covValues) - eps
        maxx <- max(covValues) + eps
        
        step <- (maxx - minx) / (optionsNumbers$knots - 1)
        knots <- seq(from=minx - optionsNumbers$degree * step,
                     to=maxx + optionsNumbers$degree * step,
                     by=step)

        ## read the coefficients samples
        coefSamples <- as.matrix(readData(sampleFile))
        
        ## create a function which returns function samples at given x values
        getFunctionSamples <- function(xValues)
        {
            ## build design matrix of basis function values
            design <- splines::spline.des(knots=knots,
                                          x=xValues,
                                          ord=optionsNumbers$degree + 1)$design
            
            ## and generate function samples
            ret <- tcrossprod(coefSamples, design)
            colnames(ret) <- xValues

            ## return
            return(convert2Mcmc(ret))
        }
         
        ## write function and variance samples into list,
        ## but also save the function
        ret[[functionName]] <-
            list(functionSamples=getFunctionSamples(covValues),
                 varianceSamples=
                 convert2Mcmc(readData(sub(pattern="(.+)_sample\\.raw", 
                                           replacement="\\1_variance_sample\\.raw",
                                           x=sampleFile))[, 1]),
                 getFunctionSamples=getFunctionSamples)
    }

    ## process deviance
    devianceInd <- grep(pattern=
                        paste(directoryWithBasename,
                              "deviance_sample\\.raw",
                              sep="_"),
                        x=resFiles)
    if(length(devianceInd))
    {
        ## read the data
        tmp <- readData(resFiles[devianceInd])

        ## get the pD
        ret$pD <- tmp["p_D", ]

        ## and the DIC
        ret$DIC <- tmp["DIC", ]
        
        # and the deviance samples (all but last two lines with pD and DIC)
        ret$Deviance <- convert2Mcmc(head(tmp, -2)) 
    }    
   
    ## process LASSO coefficients
    lassoInds <- grep(pattern=
                      paste(directoryWithBasename, 
                            "shrinkage_lasso",
                            sep="_"),
                      x=resFiles)
    if (length(lassoInds))
    {
        ret$lassoCoefficients <- convert2Mcmc(readData(resFiles[lassoInds]))
    }

    ## process Ridge coefficients
    ridgeInds <- grep(pattern=
                      paste(directoryWithBasename, 
                            "shrinkage_ridge",
                            sep="_"),
                      x=resFiles)
    if (length(ridgeInds))
    {
        ret$ridgeCoefficients <- convert2Mcmc(readData(resFiles[ridgeInds]))
    }
    ## todo: save lasso/ridge variances here as well?

    ## process scale parameter, if it exists (e.g. not for Poisson, but for Gaussian regression)
    scaleInd <- grep(pattern=
                     paste(directoryWithBasename,
                           "scale_sample\\.raw",
                           sep="_"),
                     x=resFiles)
    if(length(scaleInd))
    {
        ret$scale <- convert2Mcmc(readData(resFiles[scaleInd])[, 1])
    }
    
    ## process samples of means, if they exist
    meanInd <- grep(pattern=
                    paste(directoryWithBasename,
                          "predictmu_mean_sample\\.raw",
                          sep="_"),
                    x=resFiles)
    if(length(meanInd))
    {
        ret$means <- readData(resFiles[meanInd])
        
        ## the original names are b_1, b_2, ...
        colnames(ret$means) <- gsub(pattern="b_",
                                    replacement="",
                                    x=colnames(ret$means),
                                    fixed=TRUE)
        
        ret$means <- convert2Mcmc(ret$means)
    }

    ## which are the indexes of the samples files which have not been used yet?
    unusedInds <- setdiff(seq_along(resFiles),
                          c(devianceInd,
                            randomInds,
                            fixedInds,
                            psplineInds,
                            rwInds,
                            lassoInds,
                            ridgeInds,
                            scaleInd,
                            meanInd))

    ## if there are any unused files, extract the samples from them
    for (sampleFile in resFiles[unusedInds])
    {
        parName <- gsub(pattern=
                        paste(directoryWithBasename, 
                              "_(.+)_sample\\.raw",
                              sep = ""),
                        replacement="\\1",
                        x=sampleFile)
        
        ret[[parName]] <- convert2Mcmc(readData(sampleFile))
    }

    ## process prediction means, if they exist
    if(file.exists(predictMeanFile <- paste(directoryWithBasename,
                                            "predictmean.raw",
                                            sep="_")))
    {
        ret$PredictMeans <- read.table(predictMeanFile,
                                       header=TRUE,
                                       na.strings=c("NA", "."))
    }    
    
    ## finished!
    return(ret)
}
