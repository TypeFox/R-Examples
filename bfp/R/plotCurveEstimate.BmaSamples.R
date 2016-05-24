#####################################################################################
## Author: Daniel Sabanes Bove [daniel *.* sabanesbove *a*t* ifspm *.* uzh *.* ch]
## Project: Bayesian FPs
## 
## Time-stamp: <[plotCurveEstimate.BmaSamples.R] by DSB Mit 25/04/2012 15:54 (CEST)>
##
## Description:
## Plot predictor curve estimates based on (MC) Bayesian model average.
##
## History:
## 04/07/2008   copy from thesis function collection.
## 22/10/2008   fix indexing bug with creation of partial residuals
## 25/02/2009   don't estimate/plot the median curve, as it complicates the matter
##              without real benefit.
## 16/09/2009   add rug argument, and plot curves last to obtain clearer graphics
## 02/10/2009   adapt return list to the analogue function for BayesMfp objects
##              ("original" instead of "originalGrid")
## 03/03/2010   also accept a "main" plotting argument
## 20/09/2010   create matplotList$y in such a way that no R CMD check note is
##              triggered.
## 26/01/2011   add options "partialResids" and "hpd"
## 25/04/2012   add "median" as output
#####################################################################################

## todo: options partialResids and hpd as in hypergsplines!
plotCurveEstimate.BmaSamples <-
    function (                          # plot fp estimate, optionally with credible intervals and / or bands
              model,                    # for this BmaSamples object
              termName,                 # string denoting an fp term as written by summary
              plevel = 0.95,            # credible level for pointwise HPD, NULL means no pwise HPD
              slevel = plevel,          # credible level for simultaneous credible band, NULL means no band
              plot = TRUE,              # if FALSE, only return values to produce plot
              legendPos = "topleft",    # where to place sample size information? if NULL it is not
                                        # printed
              rug=FALSE,                # add rug to plot?
              partialResids=TRUE,
              hpd=TRUE,
              ...,                       # arguments for plotting with matplot
              main=NULL
              )
{
    ret <- list()

    if (is.null (mat <- model$bfp[[termName]]))
        stop ("There were no samples which include ", termName, " in this model average sample!\n")

    ## x values
    ret$grid <- g <- as.vector (attr (mat, "scaledGrid"))
    tr <- model$shiftScaleMax[termName, c ("shift", "scale")]
    ret$original <- g * tr[2] - tr[1]

    ## account for non-identifiability: subtract function means
    ## mat <- sweep (mat, 1, rowMeans (mat))

    ## compute pwise data
    ret$mean <- colMeans(mat, na.rm=TRUE)
    ret$median <- apply(mat, 2L, median, na.rm=TRUE)

    if (!is.null (plevel)){
        plowerUpper <-
            if(hpd)
                apply(mat, 2, empiricalHpd, level = plevel)
            else
                apply(mat, 2, quantile,
                      p=c((1 - plevel) / 2, (1 + plevel) / 2))
        
        ret$plower <- plowerUpper[1, ]
        ret$pupper <- plowerUpper[2, ]
    }

    ## simultaneous credible band around the mean
    if (!is.null (slevel)){
        bandData <-
            if(hpd)
                scrHpd(mat,
                       level = slevel,
                       mode = ret$mean)
            else
                scrBesag(mat,
                         level=slevel)

        ret$slower <- bandData[1, ]
        ret$supper <- bandData[2, ]
    }

    ## partial residuals, attention because of possible ties between observed grid values in data!
    resids <- residuals (model)
    pos <- attr (mat, "whereObsVals")

    ## obsScaledVals <- g[pos]
    ## dataSetPos <- match(model$x[, termName], obsScaledVals)
    parResids <- ret$mean[pos] + resids

    if (plot){
        ## determine plotting arguments for matlines
        matplotList <- list (...)
        if (is.null (matplotList$xlab))
            matplotList$xlab <- termName
        if (is.null (matplotList$ylab)){
            front <- paste ("Average partial predictor g(", termName, ")", sep = "")
            if (any (tr != c (0, 1))){
                middle <- " after the transform "
                back <-
                    if (tr[1] != 0){
                        if (tr[2] != 1)
                            paste(termName, "%<-% (", termName, " + ", tr[1], ") %/% ", tr[2])
                        else
                            paste(termName, "%<-%", termName, " + ", tr[1])
                    } else {
                        paste(termName, "%<-%", termName, "%/%", tr[2])
                    }
                annotation <- substitute (expression (paste (f, m, b)),
                                          list (f = front,
                                                m = middle,
                                                b = parse (text = back)[[1]]))
            } else {
                annotation <- front
            }
            matplotList$ylab <- eval (annotation)
        }
        if (is.null (matplotList$lty))
            matplotList$lty <- 1
        if (is.null (matplotList$col))
            matplotList$col <- c ("black", "gray", "blue", "blue", "green", "green")
        if (is.null (matplotList$type))
            matplotList$type <- "l"
        matplotList$x <- ret$original

        matplotList$y <- as.data.frame(ret)
        notCols <- which(names(matplotList$y) %in% c("original", "grid"))
        matplotList$y <- matplotList$y[, - notCols]
        
        if (is.null (matplotList$ylim))
            matplotList$ylim <- range (c (parResids, matplotList$y))

        ## and plot:

        ## first the points
        ret$obsVals <- ret$original[pos]
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

        ret$sampleSize <- attr (mat, "counter")
        if (!is.null (legendPos)){
            legendText <- paste ("sample size:", format (ret$sampleSize, big.mark = " "))
            legend (legendPos, legend = legendText, bty = "n")
        }
    }
    ret$partialResids <- parResids
    ret$transform <- tr

    invisible (ret)
}
