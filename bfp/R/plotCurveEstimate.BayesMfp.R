#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[plotCurveEstimate.BayesMfp.R] by DSB Mit 26/01/2011 14:32 (CET)>
##
## Description:
## Plot predictor curve estimates based on a single model.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 03/09/2008   all HPDs will be conditional on a given shrinkage factor,
##              consistent with the summary.BayesMfp method.
##              Simulations marginalizing also g can only be done with BmaSamples.
##              Center the fp part stored in xMat.
## 04/09/2008   use the center argument of getFpTransforms
## 02/10/2009   also return the transform parameters, safer indexing of them
##              inside the function
## 03/03/2010   also accept a "main" plotting argument
## 20/09/2010   create matplotList$y in such a way that no R CMD check note is
##              triggered.
## 26/01/2011   add options "partialResids" and "hpd" 
#####################################################################################

`plotCurveEstimate.BayesMfp` <-
function (                          # plot fp estimate, optionally with credible intervals and / or bands
          model,                    # for this BayesMfp object, only first model will be recognized
          termName,                 # string denoting an fp term as written by summary
          plevel = 0.95,            # credible level for pointwise HPD, NULL means no pwise
                                        # HPD
          slevel = plevel,          # credible level for simultaneous credible band, NULL means no band
          plot = TRUE,              # if FALSE, only return values to produce plot
          legendPos = "topleft",    # position where to place the mode values for the coefficients,
                                    # NULL means no leg.
          rug=FALSE,                # add rug to plot?
          partialResids=TRUE,
          hpd=TRUE,
          grid = NULL,              # vector of unscaled abscissae, default is a
                                        # length 201 grid over observed range
          post = getPosteriorParms (model), # may be computed beforehand, the shrinkage factor can
                                        # be given here in the call to getPosteriorParms
          gridSize = 201,           # obvious
          numSim = 500,             # number of simulations for simultaneous credible band
          ...,                      # arguments for plotting with matplot
          main=NULL
          )
{
    model <- model[1]
    powers <- model[[1]]$powers

    ret <- list()

    ## determine position
    fpInd <- which (attr (model, "termNames")$bfp == termName)
    if (!length (powers[[fpInd]]))
        stop ("There is no term with ", termName, " included in this model.\n")
    inds <- attr (model, "indices")
    fpCol <- inds$bfp[fpInd]

    ## x values
    tr <- attr (model, "shiftScaleMax")[fpInd, c ("shift", "scale")]
    obsVals <- attr (model, "x")[, fpCol]

    if (is.null (grid)){                # default grid
        grid <- seq (from = min (obsVals), to = max (obsVals), length = gridSize)
        ret$original <- grid * tr["scale"] - tr["shift"]
    } else {                            # scale grid
        ret$original <- grid
        grid <- (grid + tr["shift"]) / tr["scale"]
    }
    ret$grid <- grid

    ## compute pwise data
    formerIndexDesign <- length (inds$fixed) + length (unlist (powers[seq_len (fpInd - 1)]))
    indPart <- formerIndexDesign + seq_along (powers[[fpInd]])

    mStarPart <- post$mStar[indPart, drop = FALSE]
    VStarPart <- post$VStar[indPart, indPart, drop = FALSE]

    xcol <- matrix (grid, nrow = length (grid), ncol = 1, dimnames = list (NULL, termName))
    xMat <- getFpTransforms (xcol, powers[[fpInd]], center=TRUE)
    ret$mode <- as.vector (xMat %*% mStarPart)

    if (!is.null (plevel)){
        sqrtScaleDiag <- sqrt (post$bStar / post$aStar * rowSums ((xMat %*% VStarPart) * xMat))
        quant <- qt ((1 + plevel) / 2, df = 2 * post$aStar)
        ret$plower <- ret$mode - sqrtScaleDiag * quant
        ret$pupper <- ret$mode + sqrtScaleDiag * quant
    }

    ## simulate from coefficients marginal to obtain simultaneous credible band
    if (!is.null (slevel)){
        simVals <- rmvt (n = numSim, mu = mStarPart, sigma = post$bStar / post$aStar * VStarPart,
                         df = 2 * post$aStar) # coefs in rows
        simVals <- tcrossprod (simVals, xMat) # respective simulated means in cols
        bandData <-
            if(hpd)
                scrHpd(simVals,
                       level = slevel,
                       mode = ret$mode)
            else
                scrBesag(simVals,
                         level=slevel)
        
        ret$slower <- bandData[1, ]
        ret$supper <- bandData[2, ]
    }

    ## partial residuals
    design <- getDesignMatrix (model)
    resids <- residuals (model, design = design, post = post)
    parResids <- as.vector (design[, indPart, drop = FALSE] %*% mStarPart) + resids

    if (plot){
        ## determine plotting arguments
        matplotList <- list (...)
        if (is.null (matplotList$xlab))
            matplotList$xlab <- termName
        if (is.null (matplotList$ylab)){
            front <- paste ("alpha[", seq_len (ncol (xMat)), "] * ", colnames (xMat), collapse = " + ", sep = "")
            if (any (tr != c (0, 1))){
                middle <- " after the transform "
                back <-
                    if (tr["shift"] != 0){
                        if (tr["scale"] != 1)
                            paste(termName, "%<-% (", termName, " + ", tr["shift"], ") %/% ", tr["scale"])
                        else
                            paste(termName, "%<-%", termName, " + ", tr["shift"])
                    } else {
                        paste(termName, "%<-%", termName, "%/%", tr["scale"])
                    }
                annotation <- substitute (expression (paste (f, m, b)),
                                          list (f = parse (text = front)[[1]],
                                                m = middle,
                                                b = parse (text = back)[[1]]))
            } else {
                annotation <- substitute (expression (f),
                                          list (f = parse (text = front)[[1]])
                                          )
            }
            matplotList$ylab <- eval (annotation)
        }
        if (is.null (matplotList$lty))
            matplotList$lty <- 1
        if (is.null (matplotList$col))
            matplotList$col <- c ("black", "blue", "blue", "green", "green")
        if (is.null (matplotList$type))
            matplotList$type <- "l"
        matplotList$x <- ret$original

        matplotList$y <- as.data.frame(ret)
        notCols <- which(names(matplotList$y) %in% c("original", "grid"))
        matplotList$y <- matplotList$y[, - notCols]
        
        if (is.null (matplotList$ylim))
            matplotList$ylim <- range (c (parResids, matplotList$y))

        ## and plot:

        ret$obsVals <- obsVals * tr["scale"] - tr["shift"]
        
        ## first the points (the partial residuals)
        plot(ret$obsVals, parResids,
             type=if(partialResids) "p" else "n",
             xlab=matplotList$xlab,
             ylab=matplotList$ylab,
             ylim=matplotList$ylim,
             cex = 0.5,
             col = "gray",
             main=main)

        ## possibly the rug
        rug <- as.logical(rug)
        if(isTRUE(rug))
        {
            rug (jitter (ret$obsVals), col = "gray")
        }

        ## then last the curves, so that they are not over painted over by points/rug
        matplotList$add <- TRUE
        do.call (matplot, matplotList)

        if (!is.null (legendPos)){      # legend with mode of coefficients
            legend <- paste ("hat(alpha)[", seq_len (ncol (xMat)), "] == ", round (mStarPart, 5), sep = "")
            legend <- as.expression ( parse (text = legend))
            legend (legendPos, legend = legend, bty = "n")
        }
    }

    ret$partialResids <- parResids
    
    ## also save the transform parameters in the return value
    ret$transform <- tr

    ## then invisibly return 
    invisible (ret)
}

