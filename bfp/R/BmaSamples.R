#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[BmaSamples.R] by DSB Mit 25/04/2012 16:17 (CEST)>
##
## Description:
## Sample from models in a BayesMfp object using "BmaSamples" for MC model averaging
## the predictor functions. Also for prediction at new covariate points.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 04/09/2008   modified for new hyper-g methodology:
##              - add shrinkage slot to return list
##              - sample the shrinkage factors
##              - treat intercept
## 22/10/2008   fixed indexing the return bfp/uc lists: use names rather than integers,
## 24/10/2008   fixed bug in saving the bfp means of a model: transpose means matrix
## 05/11/2008   now an BayesMfp object of length 1 can be sampled as well,
##              add a newline after sampling the last model
## 09/11/2008   extend definition to enable sampling from the posterior predictive
##              distribution at new covariate points.
## 13/11/2008   use new internal function for creating model matrix for newdata
## 14/11/2008   save the random seed when entering the function
## 29/11/2008   remove unnecessary "pr" copy of priorSpecs
## 01/10/2009   coerce newdata to data.frame, so we don't need NROW here
##              (so a list can also be passed without breaking everything)
## 05/10/2009   some comments
## 05/11/2009   only construct new design matrix if there are any "newdata"!
## 16/06/2010   add predictMeans to return list, which is necessary for the log score
##              estimation (besides the already existing variance samples). The
##              actual predictive samples for "newdata" are now generated at the
##              end of the function.
## 17/06/2010   correct shifting of the new design matrix columns
## 25/04/2012   add option to include 0 samples (when the FP or linear covariate
##              is not included in the model), defaults to FALSE for backwards
##              compatibility
#####################################################################################

BmaSamples <-
    function(object,         # valid BayesMfp object including the models over which to average
             sampleSize = length (object) * 10, # sample size
             postProbs = posteriors (object),# vector of posterior probabilites that will be normalized within
             gridList = list (), # optional list of appropriately named grid vectors for fp evaluation
                                        # default is a length 201 grid per covariate additional to the observed values
                                        # (two are at the endpoints)
             gridSize = 203, # obvious. if there are many observed values, gridSize may be much
                                        # lower!
             newdata = NULL, # new covariate data with exactly the names (and preferably ranges) as before
             verbose = TRUE,  # should information on progress been printed?
             includeZeroSamples=FALSE)  # include zero samples?
{
    ## check the length of the object
    if(! (length(object) >= 1))
        stop(simpleError("There has to be at least one model in the BayesMfp object."))

    ## coerce newdata to data frame
    newdata <- as.data.frame(newdata)
    nNewObs <- nrow(newdata)
    
    ## start filling the return list
    ret <- list ()
    
    ret$priorSpecs <- attr (object, "priorSpecs")
    ret$termNames <- termNames <- attr (object, "termNames")
    ret$shiftScaleMax <- attr (object, "shiftScaleMax")
    ret$y <- attr (object, "y")
    ret$x <- attr (object, "x")
    ret$newdata <- newdata

    if(exists(".Random.seed", .GlobalEnv))
        ret$randomSeed <- get(".Random.seed", .GlobalEnv)
    
    ## extract handy things for local computations
    inds <- attr (object, "indices")
    fpIndsSeq <- seq_along (inds$bfp)
    
    colNames <- colnames (ret$x)
    nObs <- nrow(ret$x)
    yMean <- mean(ret$y)
    alpha <- ret$priorSpecs$a

    ## covariates matrix for newdata:
    if(nNewObs > 0L)
    {
        tempX <- constructNewdataMatrix(BayesMfpObject=object,
                                        newdata=newdata)
    }
    
    ## draw model indices
    objNames <- as.numeric (names (object))
    postProbs <- postProbs / sum (postProbs)
    ret$modelFreqs <-
        if(length(objNames) == 1)       # if there is only one model ...
        {
            rep(objNames, sampleSize)
        } else {                          # else more than one model ...
            sample (objNames,
                    size = sampleSize,
                    replace = TRUE,
                    prob = postProbs)
        }
    ret$modelFreqs <- table (ret$modelFreqs)
    nams <- names (ret$modelFreqs)

    ## save model summary
    ret$modelData <- as.data.frame (object)
    ret$modelData[, c ("bmaProb", "bmaFreq")] <- 0
    ret$modelData[, "bmaProb"] <- postProbs
    ret$modelData[nams, "bmaFreq"] <- ret$modelFreqs / sampleSize

    ret$sampleSize <- sampleSize

    ## reserve space for samples
    ret$sigma2 <- numeric (sampleSize)        # regression variance

    ret$shrinkage <- numeric(sampleSize) # shrinkage factor t=g/(1+g)
    
    nFix <- length (inds$fixed)
    ret$fixed <- matrix (nrow = sampleSize, ncol = nFix,
                         dimnames = list (NULL, colNames[inds$fixed])) # fixed coefficients

    ## samples of fractional polynomial function means evaluated at grids
    ## will be elements of this list:
    ret$bfp <- list ()                  
    for (i in fpIndsSeq){
        fpName <- termNames$bfp[i]

        ## determine number of samples for this FP
        m <-
            if(includeZeroSamples)
            {
                sampleSize
            } else {
                tmp <- sapply (nams,
                               function (one){
                                   as.logical (length (object[[one]]$powers[[i]]))})
        
                sum (tmp * ret$modelFreqs)
            }
        
        if (m > 0){
            ## determine additional grid values:
            obs <- ret$x[, inds$bfp[i], drop = FALSE]
            if (is.null (g <- gridList[[fpName]])){
                g <- seq (from = min (obs), to = max (obs), length = gridSize) # additional scaled grid
            }
            ## resulting total grid:
            g <- union (obs, g)
            gridSizeTotal <- length (g)

            orderG <- order(g)          # sort
            g <- g[orderG]

            mat <- matrix (nrow = m, ncol = gridSizeTotal) # for samples
            attr (mat, "whereObsVals") <- match (obs, g) # save position of observed values

            attr (mat, "scaledGrid") <- matrix (g, nrow = gridSizeTotal, ncol = 1, dimnames = list (NULL, fpName))
            attr (mat, "counter") <- 0

            ret$bfp[[fpName]] <- mat
        }
    }

    ## uncertain fixed form covariates coefficients samples
    ## will be elements of this list:
    ret$uc <- list ()                   
    for (i in seq_along (inds$ucList)){
        ucName <- termNames$uc[i]

        ## determine number of samples for this UC group
        m <-
            if(includeZeroSamples)
            {
                sampleSize
            } else {
                tmp <- sapply (nams,
                               function (one){
                                   any (object[[ one ]]$ucTerms == i)})
                sum (tmp * ret$modelFreqs)
            }
        
        if(m > 0)
        {
            mat <- matrix (nrow = m, ncol = length (inds$ucList[[i]]),
                           dimnames = list (NULL, colNames[inds$ucList[[i]]]))
            attr (mat, "counter") <- 0  # counts how many samples are already in the matrix

            ret$uc[[ucName]] <- mat
        }
    }

    ## here are the model-specific fits from all models in object:
    ret$fitted <- matrix (nrow = length(object), ncol = nObs)

    ## for samples from the posterior predictive means:
    if(nNewObs > 0L)
    {
        ret$predictMeans <-
            matrix(nrow=nNewObs,
                   ncol=sampleSize)
    }

    ## echo sampling start
    if (verbose)
        cat ("Starting sampling, current model is number ")

    ## now sample from each model as often as indicated by the modelFreqs
    ## and save fit and predictive samples
    sampleCounter <- 0      # invariant: already sampleCounter samples processed 
    ## process every model in object (to obtain fitted values along the way) 
    for (j in seq_along (object)) 
    {
        if (verbose)
            cat (j, " ")

        mod <- object[j]                # get model and

        design <- getDesignMatrix (mod) # its (centered) design matrix
                                        # (including intercept) and
        dim <- ncol(design)

        post <- getPosteriorParms (mod, design = design) # posterior parameters with posterior
                                        # expected shrinkage factor
        ret$fitted[j, ] <- fitted (mod, design = design, post = post) # to obtain fit

        if((modName <- names (mod)) %in% nams)     # sampling from this model?
        {
            m <- ret$modelFreqs[modName]
            oneInds <- seq_len (m)

            ## shrinkage factors
            shrinkage <- ret$shrinkage[sampleCounter + oneInds] <-
                if(dim > 1)
                    rshrinkage(n=m, R2=mod[[1]]$R2, nObs=nObs, p=ncol(design), alpha=alpha)
                else                    # if this is the null model: no shrinkage
                    1
                
            ## compute corresponding bStar parameters
            bStar <- (1 - shrinkage * mod[[1]]$R2) * attr(mod, "SST") / 2            

            ## regression variance
            ret$sigma2[sampleCounter + oneInds] <-
                theseVariances <-
                    rinvGamma(n=m, post$aStar, bStar)

            ## intercept
            ret$fixed[sampleCounter + oneInds, ] <-
                theseIntercepts <-
                    yMean + sqrt(bStar / post$aStar / nObs) * rt(n=m, df=nObs - 1) 
            
            ## begin prediction with these intercepts and the noise:
            if(nNewObs)
            {
                ## in each sample, all new obs have the same intercept
                ret$predictMeans[, sampleCounter + oneInds] <-
                    rep(theseIntercepts,
                        each=nNewObs)
            }
            
            ## if this is not the null model:
            if(dim > 1){

                ## then sample the effects 
                simCoefs <- matrix(data=rt(n=(dim-1) * m, df=nObs - 1),
                                   nrow=dim-1,
                                   ncol=m)
                                        # number of design columns (dim-1) x number of samples (m)
            
                simCoefs <- backsolve(r=post$XtXroot, x=simCoefs, k=dim-1)
                simCoefs <-
                    sapply(shrinkage, FUN="*", post$betaOLS) +
                        sweep(simCoefs,
                              MARGIN=2,
                              STATS=sqrt(shrinkage * bStar / post$aStar),
                              FUN="*")

                ## and sample from the likelihoods with these effects for the new covariate data
                if(nNewObs)
                {
                    ## copy model
                    tempMod <- mod
                    
                    ## correct model matrix in tempMod to new data matrix
                    attr(tempMod, "x") <- tempX

                    ## this is not necessary, because xCentered is not used by getDesignMatrix!
                    ## attr(tempMod, "xCentered") <- scale(tempX, center=TRUE, scale=FALSE)

                    ## get the correct design matrix (without the intercept column),
                    ## using the shifts of the original data!
                    newDesignNonFixed <- getDesignMatrix(tempMod,
                                                         center=FALSE)
                    newDesignNonFixed <- sweep(newDesignNonFixed,
                                               MARGIN=2L,
                                               attr(design, "shifts"))[, -1L] # intercept is discarded here 
                    
                    ## and add this to the predictive means
                    ret$predictMeans[, sampleCounter + oneInds] <-
                        ret$predictMeans[, sampleCounter + oneInds] +
                            newDesignNonFixed %*% simCoefs
                    
                }
            }

            ## sort samples into containers, from upper to lower coefficients
            rowCounter <- 0                 # invariant: already rowCounter rows of simCoefs processed         

            for (k in fpIndsSeq){           # next: fp functions means
                fpName <- termNames$bfp[k]
                at <- attributes (ret$bfp[[fpName]])                                  
                pi <- mod[[1]]$powers[[k]]
                
                if (len <- length (pi)) # if there is at least one power
                {
                    xMat <- getFpTransforms (at[["scaledGrid"]], pi, center=TRUE)
                    means <-  xMat %*% simCoefs[rowCounter + seq_len (len),, 
                                                drop = FALSE]
                    
                    rowCounter <- rowCounter + len
                } else if(includeZeroSamples) { # no power, but include zero samples

                    means <- matrix(data=0,
                                    nrow=nrow(at[["scaledGrid"]]),
                                    ncol=m)
                }

                if((len > 0L) | includeZeroSamples)
                {
                    count <- at[["counter"]]
                    ret$bfp[[fpName]][count + oneInds, ] <- t(means)
                    attr (ret$bfp[[fpName]], "counter") <- count + m
                }
            }

            for(ucInd in seq_along(termNames$uc)) # uncertain fixed form
            {
                ucName <- termNames$uc[ucInd]
                p <- ncol(ret$uc[[ucName]])

                if(ucInd %in% mod[[1]]$ucTerms) # if included
                {
                    coefSamples <- 
                        simCoefs[rowCounter + seq_len(p), ]
                    
                    rowCounter <- rowCounter + p   
                } else if(includeZeroSamples) { # not included, but include
                                        # zero samples
                    coefSamples <-
                        matrix(data=0,
                               nrow=p,
                               ncol=m)
                }

                if((ucInd %in% mod[[1]]$ucTerms) | includeZeroSamples)
                {
                    count <- attr (ret$uc[[ucName]], "counter")
                    ret$uc[[ucName]][count + oneInds, ] <- t(coefSamples)
                    attr (ret$uc[[ucName]], "counter") <- count + m
                }
            }

            sampleCounter <- sampleCounter + m # correct invariant
        }
    }
    

    ## add predictive samples (mainly for backwards compatibility)
    if(nNewObs > 0L)
    {
        ## in each sample, all new obs have the same sd, so we
        ## must replicate the corresponding sd's appropriately.
        ret$predictions <- matrix(rnorm(n=nNewObs * sampleSize,
                                        mean=ret$predictMeans,
                                        sd=
                                        rep(sqrt(ret$sigma2),
                                            each=nNewObs)),
                                  nrow=nNewObs,
                                  ncol=sampleSize)
    }
    
    ## be sure that we have a newline after the
    ## model numbers printed on the screen:
    if(verbose)
        cat("\n")

    ## attach S3 class
    class (ret) <- "BmaSamples"

    ## finally return the samples
    return (ret)
}
