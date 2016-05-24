## This function plots all graphs
##' @importFrom 'graphics' 'par' 'points' 'legend'
##' @importFrom 'stats' 'as.formula' 'predict'
##' @importFrom 'grDevices' 'dev.off'

plotAll <- function(x, dirPath = file.path(".", "figs"),
                    figArgs = list(res = 150, units = "in", height = 8, width = 8)){
    dir.create(dirPath)
    ## ---------------
    ## Set colors
    ## ---------------
    ## Colors from "set1" of RColorBrewer.
    mycols <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
                "#FF7F00", "#FFFF33", "#A65628", "#F781BF", "#999999")

    ## ---------------
    ## Plots
    ## ---------------
    allres <- x$all
    ## predictors
    prds <- rownames(allres[,,1])
    prds <- prds[prds != x$response]
    ## Responses
    rsps <- colnames(allres[,,1])
    rsps <- rsps[rsps != x$stressor]

    ## Go through pair by pair
    for(i in 1:length(prds)){
        for(j in 1:length(rsps)){
            if(i == j) next
            ## list of results of function forms on this pair
            allForms <- allres[prds[i], rsps[j],]

            ## initiate a png file
            pngArgs <- list(filename = file.path(dirPath,
                                paste0(rsps[j], "~", prds[i], ".png")))
            pngArgs <- c(pngArgs, figArgs)
            do.call(what = "png", args = pngArgs)
            
            ## plot original data
            par(mar = c(5.1, 4.1, 4.1, 11.1))
            cat("Plotting", rsps[j], "~", prds[i], "\n")
            plotargs <- list(formula = as.formula(paste0(rsps[j], "~", prds[i])),
                             data = x$data,
                             main = paste0(rsps[j], "~", prds[i],
                                 paste0("\nBest fit: ", x$best[prds[i], rsps[j]])))

            do.call(what = "plot", args = plotargs)

            ## Generate x values for fitted lines
            xs.dense <- seq(min(x$data[prds[i]], na.rm = TRUE),
                            max(x$data[prds[i]], na.rm = TRUE),
                            length = ifelse(nrow(x$data > 100), nrow(x$data), 100))
            newdata <- list(xs.dense)
            names(newdata) <- prds[i]

            ## Form 1: simple linear
            ## if(i == 1 & j == 2) browser()
            lgd.label <- NULL
            lgd.col <- NA
            nlines <- 0
            if(inherits(allForms[["SL"]], "lm")){
                ys <- predict(allForms[["SL"]], newdata = newdata)
                points(xs.dense, ys, type = "l", col = mycols[1])
                nlines <- nlines + 1
                lgd.label[nlines] <- "Simple Linear"
                lgd.col[nlines] <- mycols[1]
            }
            ## Form 2: Quadratic
            if(inherits(allForms[["Quad"]], "lm")){
                ys <- predict(allForms[["Quad"]], newdata = newdata)
                points(xs.dense, ys, type = "l", col = mycols[2])
                nlines <- nlines + 1
                lgd.label[nlines] <- "Quadratic"
                lgd.col[nlines] <- mycols[2]
            }
            ## Form 3: Simple Quadratic
            if(inherits(allForms[["SQuad"]], "lm")){
                ys <- predict(allForms[["SQuad"]], newdata = newdata)
                points(xs.dense, ys, type = "l", col = mycols[3])
                nlines <- nlines + 1
                lgd.label[nlines] <- "Simple Quadratic"
                lgd.col[nlines] <- mycols[3]
            }
            ## Form 4: Exponential
            if(inherits(allForms[["Exp"]], "lm")){
                ys <- predict(allForms[["Exp"]], newdata = newdata)
                points(xs.dense, ys, type = "l", col = mycols[4])
                nlines <- nlines + 1
                lgd.label[nlines] <- "Exponential"
                lgd.col[nlines] <- mycols[4]                
            }
            ## Form 5: log
            if(inherits(allForms[["Log"]], "lm")){
                ys <- predict(allForms[["Log"]], newdata = newdata)
                points(xs.dense, ys, type = "l", col = mycols[5])
                nlines <- nlines + 1
                lgd.label[nlines] <- "Logarithm"
                lgd.col[nlines] <- mycols[5]
            }
            ## Form 5: nls
            if(inherits(allForms[["nls"]], "nls")){
                ys <- predict(allForms[["nls"]], newdata = newdata)
                points(xs.dense, ys, type = "l", col = mycols[5])
                nlines <- nlines + 1
                lgd.label[nlines] <- "a + b * exp(c * x)"
                lgd.col[nlines] <- mycols[6]
            }
            if(length(lgd.label) > 0)
                legend("right", inset=c(-0.4,0),
                       legend = lgd.label, col = lgd.col, lty = 1, lwd = 1.5,
                       title = "Fittings", xpd = TRUE)
            dev.off()
        }
    }
    invisible(NULL)
}
