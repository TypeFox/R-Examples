################################################################################
################################################################################
# R-Funktionen zum rriskDistributions Paket (rrisk Projekt)
# Auftraggeber: Bundesinstitut für Risikobewertung, Berlin
# Auftragnehmer: Stat-Up, München
# ------------------------------------------------------------------------------
# Autor: Natalia Belgorodski, Stat-Up
# Anpassungen von Lutz Göhring, Lutz Göhring Consulting
#   Update: 05.05.2015 (Anpassungen an CRAN-Forderungen)
# Weitere Anpassungen von Matthias Flor, BfR
#   Update: 16.03.2016 (erneute Anpassungen an CRAN-Forderungen)
################################################################################
################################################################################

is.error <- function(x) inherits(x, "try-error")

#*******************************************************************************
#*******************************************************************************
# FUNCTIONS FOR FITTING BY PERCENTILES
#*******************************************************************************
#*******************************************************************************


################################################################################
################################################################################
#' This function provides a GUI for choosing a most appropriate continuous
#' distribution for known quantiles.
#'
#' The argument \code{tolPlot} defines a tolerance for plotting graphical diagnostics.
#' If the sums of the differences between the percentiles of the estimated distribution
#' and the given percentiles are smaller than this value, the distribution will be plotted.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#'
#' @name fit.perc
#' @aliases fit.perc
#' @title Choosing distribution by given quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage fit.perc(p = c(0.025, 0.5, 0.975), q = stats::qnorm(p), show.output = FALSE,
#'   tolPlot = 0.1, tolConv = 0.001, fit.weights = rep(1, length(p)))
#' @param p numerical vector of probabilities.
#' @param q numerical vector of quantiles.
#' @param show.output logical, if \code{TRUE} the output of the fitting functions \code{get.distribution.par} will be shown.
#' @param tolPlot single positive numerical value giving a tolerance for plotting graphical diagnostics. 
#'    If the sums of the differences between the distribution percentiles and
#'    the given percentiles are smaller than this value, the distribution will be plotted.
#' @param tolConv positive numerical value, the absolute convergence tolerance for reaching zero.
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @return Returns a named list containing a chosen distribution, its estimated parameters
#' and the data on which the estimation is based.
#' @note This function is used for defining a Monte-Carlo random variate item
#' (\code{mcrv}) in the \code{rrisk} project.
#' @keywords gui
#' @export
#' @examples
#' \dontrun{
#'     chosenDistr1 <- fit.perc()
#'     chosenDistr1
#'     
#'     chosenDistr2 <- fit.perc(tolPlot = 5)
#'     chosenDistr2
#'     
#'     chosenDistr3 <- fit.perc(p = c(0.3, 0.8, 0.9), q = c(10, 20, 40))
#'     chosenDistr3
#'     
#'     chosenDistr4 <- fit.perc(p = c(0.3, 0.8, 0.9), q = c(10, 30, 40))
#'     chosenDistr4
#'     
#'     chosenDistr5 <- fit.perc(p = c(0.3, 0.8, 0.9), q = c(10, 30, 40), tolPlot = 10)
#'     chosenDistr5
#'
#'     ## Fitting a PERT distribution
#'     p <- c(0.025, 0.5, 0.6, 0.975)
#'     q <- round(mc2d::qpert(p = p, min = 0, mode = 3, max = 10, shape = 5), digits = 2)
#'     chosenDistr6 <- fit.perc(p = p, q = q, tolPlot = 10)
#'     chosenDistr6
#' }
fit.perc <- function(p = c(0.025, 0.5, 0.975), 
                     q = stats::qnorm(p), 
                     show.output = FALSE, 
                     tolPlot = 0.1, tolConv = 0.001, 
                     fit.weights = rep(1, length(p))) {
    #-----------------------------------------------------------------------------
    # check consistency of the input data
    #-----------------------------------------------------------------------------
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights)) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the probability, percentiles and weights vectors are not of the same length!", call. = FALSE)
    }
    if (length(p) == 0 | length(q) == 0 | length(fit.weights) == 0) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, either the vector of probabilities, the vector of quantiles or the weights vector is empty!", call. = FALSE)
    }
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, one or more item(s) is (are) not numeric!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tolPlot) | length(tolPlot) != 1 | tolPlot < 0) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tolPlot' should be a single positive numerical value!", call. = FALSE)
    }
    if (!is.numeric(tolConv) | length(tolConv) != 1 | tolConv < 0) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tolConv' should be a single positive numerical value!", call. = FALSE)
    }
    if (any(fit.weights <= 0)) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, all items of the argument 'fit.weights' should be positive!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # load required libraries
    #-----------------------------------------------------------------------------
    #suppressWarnings(tcltk::tclRequire("Tktable"))
    if ( class(tcltk::tclRequire("Tktable")) != "tclObj" )
        stop("Tcl package \"Tktable\" required. Please install it.")
    #-----------------------------------------------------------------------------
    # define help variables for output
    #-----------------------------------------------------------------------------
    assign("tempEnvir", value = new.env())  # use default environment instead of envir = .GlobalEnv
    assign("comboDistributions", value = c(""), envir = tempEnvir)
    assign("chosenD", value = NA, envir = tempEnvir)
    assign("allParameters", value = NA, envir = tempEnvir)
    assign("fittedParams", value = NA, envir = tempEnvir)
    
    #-----------------------------------------------------------------------------
    # define help variables for tk commands and objects
    #-----------------------------------------------------------------------------
    pLabel <- "Enter probabilities, sep. by blank"
    qLabel <- "Enter percentiles, sep. by blank"
    wLabel <- "Enter weights, sep. by blank"
    tolPlotLabel <- "Plotting tolerance"
    tolConvLabel <- "Fitting tolerance"
    textWidth <- max(nchar(pLabel), nchar(qLabel), nchar(tolPlotLabel))
    headingFont1 <- tcltk::tkfont.create(size = 12, weight = "bold")
    headingFont2 <- tcltk::tkfont.create(weight = "bold", size = 10)
    pDefault <- round(p, digits = 2)
    qDefault <- round(q, digits = 2)
    wDefault <- round(fit.weights, digits = 2)
    tolPlotDefault <- tolPlot
    tolConvDefault <- tolConv
    p.tclVar <- tcltk::tclVar(pDefault)
    q.tclVar <- tcltk::tclVar(qDefault)
    w.tclVar <- tcltk::tclVar(wDefault)
    tolPlot.tclVar <- tcltk::tclVar(tolPlot)
    tolConv.tclVar <- tcltk::tclVar(tolConv)
    tclarray <- tcltk::tclArray()
    
    #-----------------------------------------------------------------------------
    # what happends by pressing "cancel" button
    #-----------------------------------------------------------------------------
    onCancel <- function(...) {
        #assign("chosenD", value = tcltk::tclvalue(tcltk::tkget(chooseCombobox)), envir = tempEnvir)
        tcltk::tkdestroy(fitpercWindow)
    } # end of onCancel()
    
    #-----------------------------------------------------------------------------
    # what happends by pressing "cancel" button
    #-----------------------------------------------------------------------------
    onOk <- function(...) {
        fittedParams.temp <- get("allParameters", envir = tempEnvir)
        if (!prod(is.na(fittedParams.temp))) {
            assign("chosenD", 
                   value = tcltk::tclvalue(tcltk::tkget(chooseCombobox)), 
                   envir = tempEnvir)
            chosenD <- get("chosenD", envir = tempEnvir)
            if (nchar(chosenD) > 0) {
                fittedParams <- fittedParams.temp[colnames(fittedParams.temp) == chosenD]
                
                if (colnames(fittedParams) == "norm") {
                    fittedParams <- c(mean = fittedParams[1, 1], 
                                      sd = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "beta") {
                    fittedParams <- c(shape1 = fittedParams[1, 1], 
                                      shape2 = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "cauchy") {
                    fittedParams <- c(location = fittedParams[1, 1], 
                                      scale = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "logis") {
                    fittedParams <- c(location = fittedParams[1, 1], 
                                      scale = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "t") {
                    fittedParams <- c(df = fittedParams[1, 1])
                } else if (colnames(fittedParams) == "chisq") {
                    fittedParams <- c(df = fittedParams[1, 1]) 
                } else if (colnames(fittedParams) == "chisqnc") {
                    fittedParams <- c(df = fittedParams[1, 1], 
                                      ncp = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "exp") {
                    fittedParams <- c(rate = fittedParams[1, 1])
                } else if (colnames(fittedParams) == "f") {
                    fittedParams <- c(df1 = fittedParams[1, 1], 
                                      df2 = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "gamma") {
                    fittedParams <- c(shape = fittedParams[1, 1], 
                                      rate = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "lnorm") {
                    fittedParams <- c(meanlog = fittedParams[1, 1], 
                                      sdlog = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "unif") {
                    fittedParams <- c(min = fittedParams[1, 1], 
                                      max = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "weibull") {
                    fittedParams <- c(shape = fittedParams[1, 1], 
                                      scale = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "triang") {
                    fittedParams <- c(min = fittedParams[1, 1],
                                      mode = fittedParams[2, 1], 
                                      max = fittedParams[3, 1])
                } else if (colnames(fittedParams) == "gompertz") {
                    fittedParams <- c(shape = fittedParams[1, 1],
                                      scale = fittedParams[2, 1])
                } else if (colnames(fittedParams) == "pert") {
                    fittedParams <- c(min = fittedParams[1, 1], 
                                      mode = fittedParams[2, 1], 
                                      max = fittedParams[3, 1], 
                                      shape = fittedParams[3, 1])
                } else if (colnames(fittedParams) == "tnorm") {
                    fittedParams <- c(mean = fittedParams[1, 1], 
                                      sd = fittedParams[2, 1], 
                                      lower = fittedParams[3, 1], 
                                      upper = fittedParams[4, 1])
                }
                assign("fittedParams", 
                       value = fittedParams, 
                       envir = tempEnvir)
            } else chosenD <- "NA"
        }
        tcltk::tkdestroy(fitpercWindow)
    } # end of onOK()
    
    #-----------------------------------------------------------------------------
    # what happends by pressing "reset" button
    #-----------------------------------------------------------------------------
    onReset <- function(...) {
        p.tclVar <- tcltk::tclVar(pDefault)
        q.tclVar <- tcltk::tclVar(qDefault)
        w.tclVar <- tcltk::tclVar(wDefault)
        tolPlot.tclVar <- tcltk::tclVar(tolPlotDefault)
        tolConv.tclVar <- tcltk::tclVar(tolConvDefault)
        tcltk::tkconfigure(pEntry, text = p.tclVar)
        tcltk::tkconfigure(qEntry, text = q.tclVar)
        tcltk::tkconfigure(wEntry, text = w.tclVar)
        tcltk::tkconfigure(tolPlotEntry, text = tolPlot.tclVar)
        tcltk::tkconfigure(tolConvEntry, text = tolConv.tclVar)
    }  # end of onReset()
    
    #-----------------------------------------------------------------------------
    # what happends by pressing "clear" button
    #-----------------------------------------------------------------------------
    onClear <- function(...) {
        p.tclVar <- tcltk::tclVar("")
        q.tclVar <- tcltk::tclVar("")
        w.tclVar <- tcltk::tclVar("")
        tolPlot.tclVar <- tcltk::tclVar("")
        tolConv.tclVar <- tcltk::tclVar("")
        tcltk::tkconfigure(pEntry, text = p.tclVar)
        tcltk::tkconfigure(qEntry, text = q.tclVar)
        tcltk::tkconfigure(wEntry, text = q.tclVar)
        tcltk::tkconfigure(tolPlotEntry, text = tolPlot.tclVar)
        tcltk::tkconfigure(tolConvEntry, text = tolConv.tclVar)
    } # end of onClear()
    
    #-----------------------------------------------------------------------------
    # what happends by pressing "fit" button
    #-----------------------------------------------------------------------------
    onFit <- function(...) {
        # define help variable for p
        pTemp <- tcltk::tclvalue(tcltk::tkget(pEntry))
        pTemp <- strsplit(pTemp, " ")[[1]]
        pTemp <- c(apply(matrix(pTemp, nrow = 1), 1, 
                         function(x) sub(x, pattern = " ", replacement = "")))
        toRemove <- which(pTemp == "")
        if (length(toRemove > 0)) pTemp <- pTemp[-toRemove]
        p <- as.numeric(pTemp)
        
        # define help variables for q
        qTemp <- tcltk::tclvalue(tcltk::tkget(qEntry))
        qTemp <- strsplit(qTemp, " ")[[1]]
        qTemp <- c(apply(matrix(qTemp, nrow = 1), 1, 
                         function(x) sub(x, pattern = " ", replacement = "")))
        toRemove <- which(qTemp == "")
        if (length(toRemove > 0)) qTemp <- qTemp[-toRemove]
        q <- as.numeric(qTemp)
        
        # define help variables for w(fit.weights)
        wTemp <- tcltk::tclvalue(tcltk::tkget(wEntry))
        wTemp <- strsplit(wTemp, " ")[[1]]
        wTemp <- c(apply(matrix(wTemp, nrow = 1), 1,
                         function(x) sub(x, pattern = " ", replacement = "")))
        toRemove <- which(wTemp == "")
        if (length(toRemove > 0)) wTemp <- wTemp[-toRemove]
        fit.weights <- as.numeric(wTemp)
        
        # define help variable for tolPlot
        tolPlotTemp <- tcltk::tclvalue(tcltk::tkget(tolPlotEntry))
        tolPlotTemp <- strsplit(tolPlotTemp, " ")[[1]]
        tolPlotTemp <- c(apply(matrix(tolPlotTemp, nrow = 1), 1,
                               function(x) sub(x, pattern = " ", replacement = "")))
        toRemove <- which(tolPlotTemp == "")
        if (length(toRemove > 0)) tolPlotTemp <- tolPlotTemp[-toRemove]
        tolPlot <- as.numeric(tolPlotTemp)
        
        # define help variable for tolConv
        tolConvTemp <- tcltk::tclvalue(tcltk::tkget(tolConvEntry))
        tolConvTemp <- strsplit(tolConvTemp, " ")[[1]]
        tolConvTemp <- c(apply(matrix(tolConvTemp, nrow = 1), 1, 
                               function(x) sub(x, pattern = " ", replacement = "")))
        toRemove <- which(tolConvTemp == "")
        if (length(toRemove > 0)) tolConvTemp <- tolConvTemp[-toRemove]
        tolConv <- as.numeric(tolConvTemp)
        
        # check input data consistency
        if (length(p) == 0) {
            tcltk::tkmessageBox(message = "The inputs are empty! Please correct your input!", icon = "error")
            tcltk::tkfocus(pEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if (length(q) == 0) {
            tcltk::tkmessageBox(message = "The inputs are empty! Please correct your input!", icon = "error")
            tcltk::tkfocus(qEntry)
            stop("")
        }
        if (length(fit.weights) == 0) {
            tcltk::tkmessageBox(message = "The inputs are empty! Please correct your input!", icon = "error")
            tcltk::tkfocus(wEntry)
            stop("")
        }
        if (any(is.na(p))) {
            tcltk::tkmessageBox(message = "One or more item(s) is not numeric, please correct your input!", icon = "error")
            tcltk::tkfocus(pEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if ( any(is.na(q))) {
            tcltk::tkmessageBox(message = "One or more item(s) is not numeric, please correct your input!", icon = "error")
            tcltk::tkfocus(qEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if ( any(is.na(fit.weights))) {
            tcltk::tkmessageBox(message = "One or more item(s) is not numeric, please correct your input!", icon = "error")
            tcltk::tkfocus(wEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if ( any((fit.weights) <= 0)) {
            tcltk::tkmessageBox(message = "All weights should be positive, please correct your input!", icon = "error")
            tcltk::tkfocus(wEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights)) {
            tcltk::tkmessageBox(message = "The vectors of percentiles, probabilities and weights are not of the same lengths! Please correct your input!", icon = "error")
            tcltk::tkfocus(pEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if (prod(order(p) == seq(1:length(p))) == 0) {
            tcltk::tkmessageBox(message = "The vector of probabilities ist not ordered! Please correct your input!", icon = "error")
            tcltk::tkfocus(pEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if (prod(order(q) == seq(1:length(q))) == 0) {
            tcltk::tkmessageBox(message = "The vector of quantiles ist not ordered! Please correct your input!", icon = "error")
            tcltk::tkfocus(qEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if (min(p) < 0 | max(p) > 1) {
            tcltk::tkmessageBox(message = "Items of the probability vector should lie between 0 and 1! Please correct your input!", icon = "error")
            tcltk::tkfocus(pEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if (!prod(is.numeric(tolPlot)) | length(tolPlot) != 1 | any(tolPlot < 0)) {
            tcltk::tkmessageBox(message = "The tolerance for diagnostic plotting should be a single positive numerical value! Please correct your input!", icon = "error")
            tcltk::tkfocus(tolPlotEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        if (!prod(is.numeric(tolConv)) | length(tolConv) != 1 | any(tolConv < 0)) {
            tcltk::tkmessageBox(message = "The tolerance for distributions fitting should be a single positive numerical value! Please correct your input!", icon = "error")
            tcltk::tkfocus(tolConvEntry)
            stop("INVALID INPUT", call. = FALSE)
        }
        
        # clear image, table and combobox
        tcltk::tkconfigure(chooseCombobox, values = c(""))
        tcltk::tkset(chooseCombobox, c(""))
        tcltk::tkconfigure(fitResultTable, variable = tcltk::tclArray())
        tkrplot::tkrreplot(imgPlot, hscale = 1.4, vscale = 1.2,
                           fun = function() {
                               graphics::plot(stats::rnorm(20), 
                                              col = "white", 
                                              xlab = "Percentile", ylab = "Percent", 
                                              main = "Graphical diagnostics")
                           }
        )
        
        # calculate results matrix
        fit.results <- rriskFitdist.perc(p, q,
                                         show.output = FALSE, 
                                         tolConv = tolConv, fit.weights)
        
        if (!prod(is.na(fit.results))) { # if res.matrix is not empty
            res.matrix <- fit.results$results
            assign("allParameters", value = res.matrix[1:4, -c(1, 2)], envir = tempEnvir)
            comboDistributions <- colnames(res.matrix)[!apply(res.matrix, 2, 
                                                              function(x) ifelse(all(is.na(x)), TRUE, FALSE))]
            comboDistributions <- comboDistributions[-c(1, 2)]
            assign("comboDistributions", value = comboDistributions, envir = tempEnvir)
            
            # formatting results
            res.matrix <- as.matrix(
                apply(res.matrix, 
                      c(1, 2), 
                      function(x) {
                          if (!is.na(x)) {
                              if (x > 10000) {
                                  x <- "+infty"
                              } else if (x < (-10000)) {
                                  x <- "-infty"
                              }
                              else x <- x
                          } 
                      }
                )
            )
            res.matrix <- cbind(rownames(res.matrix), res.matrix)
            res.matrix <- rbind(colnames(res.matrix), res.matrix)
            res.matrix[1, 1] <- "Percent"
            
            # define data for tktable
            for (i in 0:(nrow(res.matrix) - 1)) {
                for (j in 0:(ncol(res.matrix) - 1)) {
                    temp <- unlist(res.matrix[i + 1, j + 1])
                    tclarray[[i, j]] <- temp
                }
            }
            # changed by LG: cf. comment in method plotDiagnostics.perc() where the message is produced:
            withCallingHandlers(
                tkrplot::tkrreplot(imgPlot,
                                   fun = function() plotDiagnostics.perc(fit.results, tolPlot = tolPlot),
                                   hscale = 1.4, vscale = 1.2),
                message = function(m) tcltk::tkmessageBox(message = m$message, icon = "error")
            )
            tcltk::tkconfigure(fitResultTable, variable = tclarray, rows = nrow(res.matrix))
            tcltk::tkconfigure(chooseCombobox, values = get("comboDistributions", envir = tempEnvir))
            tcltk::tkset(chooseCombobox, get("comboDistributions", envir = tempEnvir)[1])
        } else {# if results matrix is empty
            allParameters.temp <- get("allParameters", envir = tempEnvir)
            allParameters.temp <- apply(allParameters.temp, c(1, 2), function(x) x <- NA)
            assign("allParameters", value = allParameters.temp, envir = tempEnvir)
            tcltk::tkraise(fitpercWindow)
            tcltk::tkmessageBox(message = "Sorry, no pdfs found that match with the specified percentiles!", icon = "error")
            tcltk::tkfocus(pEntry)
            stop("no pdfs found that match with the specified percentiles")
        }
        
        #to.remove <- as.matrix(apply(res.mat, 2, function(x) ifelse(all(is.na(x)), TRUE, FALSE)))
        #res.mat <- res.mat[,which(to.remove == FALSE)]
        tcltk::tkraise(fitpercWindow)
    } # end of onFit()
    
    #-----------------------------------------------------------------------------
    # create GUI window and frames
    #-----------------------------------------------------------------------------
    fitpercWindow <- tcltk::tktoplevel(width = 860, height = 580)
    tcltk::tkwm.title(fitpercWindow, "Fitting continuous distributions to given percentiles")
    tcltk::tkwm.resizable(fitpercWindow, FALSE, FALSE)  # fixed size, not resizeable
    #tcltk::tkwm.maxsize(fitpercWindow, 880, 580)
    #tcltk::tkwm.minsize(fitpercWindow, 880, 580)
    allFrame <- tcltk::tkframe(fitpercWindow)
    InputImageFrame <- tcltk::tkframe(allFrame)
    
    #-----------------------------------------------------------------------------
    # create contents of input frame
    #-----------------------------------------------------------------------------
    fitFrame <- tcltk::tkframe(InputImageFrame, relief = "groove", borderwidth = 4)
    inputFrame <- tcltk::tkframe(fitFrame)
    
    tcltk::tkpack(tcltk::tklabel(inputFrame, text = "Input frame", font = headingFont1), pady = c(0, 10))
    tcltk::tkpack(tcltk::tklabel(inputFrame, text = pLabel, width = 30), pady = c(0, 10))
    pEntry <- tcltk::tkentry(inputFrame, width = 30, text = p.tclVar, textvariable = p.tclVar)
    tcltk::tkpack(pEntry, pady = c(0, 10))
    
    tcltk::tkpack(tcltk::tklabel(inputFrame, text = qLabel, width = 30))
    qEntry <- tcltk::tkentry(inputFrame, width = 30, text = q.tclVar, textvariable = q.tclVar)
    tcltk::tkpack(qEntry, pady = c(0, 10))
    
    tcltk::tkpack(tcltk::tklabel(inputFrame, text = wLabel, width = 30))
    wEntry <- tcltk::tkentry(inputFrame, width = 30, text = w.tclVar, textvariable = w.tclVar)
    tcltk::tkpack(wEntry, pady = c(0, 10))
    
    tcltk::tkpack(tcltk::tklabel(inputFrame, text = tolPlotLabel, width = 30))
    tolPlotEntry <- tcltk::tkentry(inputFrame, width = 30, text = tolPlot.tclVar, textvariable = tolPlot.tclVar)
    tcltk::tkpack(tolPlotEntry)
    tcltk::tkpack(inputFrame, side = "top", padx = c(5, 5), pady = c(5, 5))
    tcltk::tkpack(fitFrame)
    
    tcltk::tkpack(tcltk::tklabel(inputFrame, text = tolConvLabel, width = 30))
    tolConvEntry <- tcltk::tkentry(inputFrame, width = 30, text = tolConv.tclVar, textvariable = tolConv.tclVar)
    tcltk::tkpack(tolConvEntry)
    tcltk::tkpack(inputFrame, side = "top", padx = c(5, 5), pady = c(5, 5))
    tcltk::tkpack(fitFrame)
    
    buttonsFrame1 <- tcltk::tkframe(fitFrame)
    fitButton <- tcltk::ttkbutton(buttonsFrame1, width = 8, text = "Fit", command = onFit)
    clearButton <- tcltk::ttkbutton(buttonsFrame1, width = 8, text = "Clear", command = onClear)
    resetButton <- tcltk::ttkbutton(buttonsFrame1, width = 8, text = "Reset", command = onReset)
    tcltk::tkpack(fitButton, side = "left", padx = c(0, 10))
    tcltk::tkpack(clearButton, side = "left", padx = c(0, 10))
    tcltk::tkpack(resetButton, side = "left", padx = c(0, 0))
    tcltk::tkpack(buttonsFrame1, side = "bottom", pady = c(0, 15))
    tcltk::tkpack(fitFrame, side = "right")
    
    #-----------------------------------------------------------------------------
    # create contents of image frame
    #-----------------------------------------------------------------------------
    imgPlot <- tkrplot::tkrplot(InputImageFrame, hscale = 1.4, vscale = 1.2,
                                fun = function() {
                                    graphics::plot(stats::rnorm(20), 
                                                   col = "white", 
                                                   xlab = "Percentile", ylab = "Percent", 
                                                   main = "Graphical diagnostics")
                                })
    tcltk::tkpack(imgPlot, side = "left")
    tcltk::tkpack(InputImageFrame)
    
    #-----------------------------------------------------------------------------
    # create table frame
    #-----------------------------------------------------------------------------
    TableFrame <- tcltk::tkframe(allFrame)
    fitResultTable <- tcltk::tkwidget(TableFrame, "table", variable = tcltk::tclArray(),
                                      height = 10, rows = 12, cols = 20, 
                                      background = "white", borderwidth = 2, state = "disabled", 
                                      titlerows = 1, titlecols = 1, 
                                      resizeborders = "none", colwidth = 6, maxwidth = 1200,
                                      yscrollcommand = function(...) tcltk::tkset(yscr, ...), 
                                      selectmode = "extended")
    yscr <- tcltk::tkscrollbar(TableFrame, 
                               command = function(...) tcltk::tkyview(fitResultTable, ...))
    tcltk::tkgrid(fitResultTable, yscr, sticky = "ns")
    tcltk::tkpack(TableFrame)
    
    #-----------------------------------------------------------------------------
    # create buttons and combobox frame
    #-----------------------------------------------------------------------------
    buttonsFrame2 <- tcltk::tkframe(allFrame)
    okButton <- tcltk::ttkbutton(buttonsFrame2, 
                                 width = 10, text = "Ok",
                                 command = onOk)
    cancelButton <- tcltk::ttkbutton(buttonsFrame2, width = 10, text = "Cancel", 
                                     command = onCancel)
    chooseCombobox <- tcltk::ttkcombobox(buttonsFrame2, 
                                         values = get("comboDistributions", 
                                                      envir = tempEnvir), 
                                         width = 10, state = "readonly")
    tcltk::tkpack(tcltk::tklabel(buttonsFrame2, 
                                 text = "chosen distribution", 
                                 font = headingFont2, width = 15), 
                  side = "left", padx = c(0, 15))
    tcltk::tkpack(chooseCombobox, side = "left", padx = c(0, 200))
    tcltk::tkpack(okButton, side = "left", padx = c(0, 15))
    tcltk::tkpack(cancelButton, side = "left", padx = c(0, 15))
    tcltk::tkpack(buttonsFrame2, pady = c(10, 0))
    
    #-----------------------------------------------------------------------------
    # place allFrame on GUI window
    #-----------------------------------------------------------------------------
    tcltk::tkpack(allFrame, padx = c(15, 15), pady = c(0, 10))
    
    #-----------------------------------------------------------------------------
    # create matrix with fitting results
    #-----------------------------------------------------------------------------
    fit.results <- rriskFitdist.perc(pDefault, qDefault, 
                                     show.output = FALSE, tolConv = tolConv, 
                                     fit.weights = wDefault)
    
    #-----------------------------------------------------------------------------
    # fill image and table with fitting results
    #-----------------------------------------------------------------------------
    if (!prod(is.na(fit.results))) { # if res.matrix is not empty
        res.matrix <- fit.results$results
        assign("allParameters", value = res.matrix[1:4,-c(1, 2)], envir = tempEnvir)
        comboDistributions <- colnames(res.matrix)[!apply(res.matrix, 2, 
                                                          function(x) ifelse(all(is.na(x)), TRUE, FALSE))]
        comboDistributions <- comboDistributions[-c(1, 2)]
        assign("comboDistributions", value = comboDistributions, envir = tempEnvir)
        tcltk::tkconfigure(chooseCombobox, values = get("comboDistributions", envir = tempEnvir))
        tcltk::tkset(chooseCombobox, get("comboDistributions", envir = tempEnvir)[1])
        #-----------------------------------------------------------------------------
        # formatting results
        #-----------------------------------------------------------------------------
        res.matrix <- as.matrix(apply(res.matrix, c(1, 2),
                                      function(x) {
                                          if (!is.na(x)) {
                                              if (x > 10000) {
                                                  x <- "+infty"
                                              } else if (x < (-10000)) {
                                                  x <- "-infty"
                                              } else x <- x
                                          } 
                                      }
        ))
        res.matrix <- cbind(rownames(res.matrix), res.matrix)
        res.matrix <- rbind(colnames(res.matrix), res.matrix)
        res.matrix[1, 1] <- "Percent"
        tableRows <- nrow(res.matrix)
        
        #-----------------------------------------------------------------------------
        # create data (tclarray) for tktable
        #-----------------------------------------------------------------------------
        for (i in 0:(nrow(res.matrix) - 1)) {
            for (j in 0:(ncol(res.matrix) - 1)) {
                temp <- unlist(res.matrix[i + 1, j + 1])
                tclarray[[i, j]] <- temp
            }
        }
        # changed by LG: cf. comment in method plotDiagnostics.perc() where the message is produced:
        withCallingHandlers(
            tkrplot::tkrreplot(
                imgPlot,
                fun = function() plotDiagnostics.perc(fit.results, tolPlot = tolPlot),
                hscale = 1.4, vscale = 1.2),
            message = function(m) tcltk::tkmessageBox(message = m$message, icon = "error")
        )
        tcltk::tkconfigure(fitResultTable, variable = tclarray, rows = tableRows)
    } else {
        tcltk::tkmessageBox(
            message = "Sorry, no pdfs found that match with the specified percentiles!", 
            icon = "error")
        tcltk::tkfocus(pEntry)
    }
    
    tcltk::tkfocus(fitpercWindow)
    tcltk::tkwait.window(fitpercWindow)
    
    #-----------------------------------------------------------------------------
    # generate output
    #-----------------------------------------------------------------------------
    chosenD <- get("chosenD", envir = tempEnvir)
    fittedParams <- get("fittedParams", envir = tempEnvir)
    output <- list(data.frame(p, q), chosenD, fittedParams)
    names(output) <- c("p/q", "chosenDistr", "fittedParams")
    
    #-----------------------------------------------------------------------------
    # output message
    #-----------------------------------------------------------------------------
    print.on.exit <- function(chosenD) {
        if (!is.na(chosenD)) {
            exitMessage <- ""
            if (chosenD == "norm")           exitMessage <- "Normal (norm)"
            else if (chosenD == "beta")      exitMessage <- "Beta (beta)"
            else if (chosenD == "cauchy")    exitMessage <- "Cauchy (cauchy)"
            else if (chosenD == "logis")     exitMessage <- "Logistic (logis)"
            else if (chosenD == "t")         exitMessage <- "Student's t (t)"
            else if (chosenD == "chisq")     exitMessage <- "Chi-squared (chisq)"
            else if (chosenD == "chisqnc")   exitMessage <- "Non-central chi-squared (chisqnc)"
            else if (chosenD == "exp")       exitMessage <- "Exponential (exp)"
            else if (chosenD == "f")         exitMessage <- "F (f)"
            else if (chosenD == "gamma")     exitMessage <- "Gamma (gamma)"
            else if (chosenD == "lnorm")     exitMessage <- "Log-normal (lnorm)"
            else if (chosenD == "unif")      exitMessage <- "Uniform (unif)"
            else if (chosenD == "weibull")   exitMessage <- "Weibull (weibull)"
            else if (chosenD == "triang")    exitMessage <- "Triangular (triang)"
            else if (chosenD == "gompertz")  exitMessage <- "Gompertz (gompertz)"
            else if (chosenD == "pert")      exitMessage <- "Beta-pert (pert)"
            else if (chosenD == "tnorm")     exitMessage <- "Truncated normal (tnorm)"
        } else exitMessage <- chosenD
        cat(paste("Chosen continuous distribution is:", exitMessage))
        cat("\nFitted parameters are: \n")
        print(fittedParams)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    on.exit(print.on.exit(chosenD))
    return(invisible(output))
} # end of function fit.perc()



################################################################################
#' Diagnostic plot for choosing a most appropriate continuous probability for known quantiles
#'
#' This function plots distribution whose percentiles go through the given percentiles
#' \code{q}. The argument \code{tolPlot} controls this match.
#'
#' @name plotDiagnostics.perc
#' @aliases plotDiagnostics.perc
#' @title Graphical tools for choosing distribution by given quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage plotDiagnostics.perc(fit.results, tolPlot = 0.1)
#' @param fit.results a list containing fitting results as an output of the function \code{rriskFitdist.perc}.
#' @param tolPlot numerical value, if the sums of the differences between the distribution
#' percentiles and the given percentiles are smaller than this value, the distribution
#' will be plotted.
#' @return Only graphical output.
#' @keywords others
#' @export
#' @importFrom mc2d qtriang
#' @importFrom mc2d ptriang
#' @importFrom mc2d qpert
#' @importFrom mc2d ppert
#' @importFrom eha qgompertz
#' @importFrom eha pgompertz
#' @importFrom msm qtnorm
#' @importFrom msm ptnorm
#' @examples
#' p <- c(0.025, 0.5, 0.975)
#' q <- c(9.68, 29.20, 50.98)
#' fit.results1 <- rriskFitdist.perc(p = p, q = q, show.output = FALSE, tolConv = 0.5)
#' old.par <- graphics::par(mfrow = c(1, 2))
#' plotDiagnostics.perc(fit.results1)
#' plotDiagnostics.perc(fit.results1, tolPlot = 5)
#' graphics::par(old.par)
#'
#' p <- c(0.2, 0.7)
#' q <- c(2.6, 19.1)
#' fit.results2 <- rriskFitdist.perc(p = p, q = q, show.output = FALSE)
#' plotDiagnostics.perc(fit.results2)
#'
#' p <- c(0.3, 0.8, 0.9)
#' q <- c(10, 20, 40)
#' fit.results3 <- rriskFitdist.perc(p = p, q = q, show.output = FALSE)
#' plotDiagnostics.perc(fit.results3)
#'
#' p <- c(0.3, 0.8, 0.9)
#' q <- c(10, 30, 40)
#' fit.results4 <- rriskFitdist.perc(p = p, q = q, show.output = FALSE)
#' plotDiagnostics.perc(fit.results4)
#'
#' ## Example with fitted beta pert distribution
#' p <- c(0.025, 0.5, 0.6, 0.975)
#' q <- mc2d::qpert(p = p, min = 0, mode = 3, max = 10, shape = 5)
#' fit.results5 <- rriskFitdist.perc(p = p, q = q, show.output = FALSE)
#' plotDiagnostics.perc(fit.results5)
#' 
plotDiagnostics.perc <- function(fit.results, tolPlot = 0.1) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (prod(is.na(fit.results))) {
        stop("INVALID INPUT, the argument 'fit.results' is empty!", call. = FALSE)
    }
    if (!is.numeric(tolPlot) | length(tolPlot) != 1 | tolPlot < 0) {
        stop("INVALID INPUT, the argument 'tolPlot' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # creating help variables
    #-----------------------------------------------------------------------------
    res.mat <- fit.results$results
    p <- fit.results$`p/q`$p
    q <- fit.results$`p/q`$q
    len <- max(q) - min(q)
    min <- min(q) - len * 1.5
    max <- max(q) + len * 1.5
    
    #-----------------------------------------------------------------------------
    # plotting quantile points
    #-----------------------------------------------------------------------------
    xx <- seq(from = min, to = max, length = 300)
    leg.col <- NULL
    leg.txt <- NULL
    graphics::plot(p ~ q, type = "p", cex = 1.5, pch = 19, lwd = 2, xlim = c(min, max),
                   ylim = c(0, 1), xlab = "Percentile", ylab = "Percent", col = "lightblue",
                   main = "Graphical diagnostics")
    
    #-----------------------------------------------------------------------------
    # defining color for each distribution
    #-----------------------------------------------------------------------------
    norm.color <- "blue"
    beta.color <- "red"
    cauchy.color <- "black"
    chisq.color <- "magenta"
    logis.color <- "orange"
    t.color <- "lightpink4"
    f.color <- "lightgoldenrod4"
    gamma.color <- "olivedrab"
    lnorm.color <- "navajowhite3"
    weibull.color <- "darkorchid4"
    exp.color <- "yellow3"
    triang.color <- "green4"
    gompertz.color <- "peru"
    unif.color <- "darkturquoise"
    pert.color <- "hotpink4"
    tnorm.color <- "gold2"
    chisqnc.color <- "gold4"
    
    #-----------------------------------------------------------------------------
    # plotting chisqnc distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$chisqnc[1])) {
        test.quantiles <- suppressWarnings(stats::qchisq(p, df = res.mat$chisqnc[1], 
                                                         ncp = res.mat$chisqnc[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pchisq(xx, df = res.mat$chisqnc[1], 
                                          ncp = res.mat$chisqnc[2]) ~ xx, 
                            col = chisqnc.color, lwd = 2)
            leg.txt <- c(leg.txt, "n.-c. chi-square")
            leg.col <- c(leg.col, chisqnc.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # plotting tnorm distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$tnorm[1])) {
        test.quantiles <- suppressWarnings(msm::qtnorm(p, 
                                                       mean = res.mat$tnorm[1], 
                                                       sd = res.mat$tnorm[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(msm::ptnorm(xx, 
                                        mean = res.mat$tnorm[1], 
                                        sd = res.mat$tnorm[2],
                                        lower = res.mat$tnorm[3], 
                                        upper = res.mat$tnorm[4]) ~ xx, 
                            col = tnorm.color, lwd = 2)
            leg.txt <- c(leg.txt, "trunc. normal")
            leg.col <- c(leg.col, tnorm.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # plotting pert distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$pert[1])) {
        test.quantiles <- suppressWarnings(mc2d::qpert(p, 
                                                       min = res.mat$pert[1], 
                                                       mode = res.mat$pert[2], 
                                                       max = res.mat$pert[3], 
                                                       shape = res.mat$pert[4]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(mc2d::ppert(xx, 
                                        min = res.mat$pert[1], 
                                        mode = res.mat$pert[2],
                                        max = res.mat$pert[3], 
                                        shape = res.mat$pert[4]) ~ xx, 
                            col = pert.color, lwd = 2)
            leg.txt <- c(leg.txt, "Beta pert")
            leg.col <- c(leg.col, pert.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting triang distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$triang[1])) {
        test.quantiles <- suppressWarnings(mc2d::qtriang(p, 
                                                         min = res.mat$triang[1], 
                                                         mode = res.mat$triang[2], 
                                                         max = res.mat$triang[3]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(mc2d::ptriang(xx, 
                                          min = res.mat$triang[1], 
                                          mode = res.mat$triang[2],
                                          max = res.mat$triang[3]) ~ xx, 
                            col = triang.color, lwd = 2)
            leg.txt <- c(leg.txt, "Triangular")
            leg.col <- c(leg.col, triang.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting gompertz distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$gompertz[1])) {
        test.quantiles <- suppressWarnings(eha::qgompertz(p, 
                                                          shape = res.mat$gompertz[1] + 1/10000, 
                                                          scale = res.mat$gompertz[2] + 1/10000))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(eha::pgompertz(xx, 
                                           shape = res.mat$gompertz[1] + 1/10000, 
                                           scale = res.mat$gompertz[2] + 1/10000) ~ xx, 
                            col = gompertz.color, lwd = 2)
            leg.txt <- c(leg.txt, "Gompertz")
            leg.col <- c(leg.col, gompertz.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting normal distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$norm[1])) {
        test.quantiles <- suppressWarnings(stats::qnorm(p, 
                                                        mean = res.mat$norm[1], 
                                                        sd = res.mat$norm[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pnorm(xx, 
                                         mean = res.mat$norm[1], 
                                         sd = res.mat$norm[2]) ~ xx, 
                            col = norm.color, lwd = 2)
            leg.txt <- c(leg.txt, "Normal")
            leg.col <- c(leg.col, norm.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting beta ditribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$beta[1])) {
        test.quantiles <- suppressWarnings(stats::qbeta(p, 
                                                        shape1 = res.mat$beta[1], 
                                                        shape2 = res.mat$beta[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pbeta(xx, 
                                         shape1 = res.mat$beta[1], 
                                         shape2 = res.mat$beta[2]) ~ xx, 
                            col = beta.color, lwd = 2)
            leg.txt <- c(leg.txt, "Beta")
            leg.col <- c(leg.col, beta.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting cauchy distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$cauchy[1])) {
        test.quantiles <- suppressWarnings(stats::qcauchy(p, 
                                                          location = res.mat$cauchy[1], 
                                                          scale = res.mat$cauchy[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pcauchy(xx, 
                                           location = res.mat$cauchy[1], 
                                           scale = res.mat$cauchy[2]) ~ xx, 
                            col = cauchy.color, lwd = 2)
            leg.txt <- c(leg.txt, "Cauchy")
            leg.col <- c(leg.col, cauchy.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting chisq distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$chisq[1])) {
        test.quantiles <- suppressWarnings(stats::qchisq(p, df = res.mat$chisq[1]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pchisq(xx, df = res.mat$chisq[1]) ~ xx, 
                            col = chisq.color, lwd = 2)
            leg.txt <- c(leg.txt, "Chi-square")
            leg.col <- c(leg.col, chisq.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting logis ditribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$logis[1])) {
        test.quantiles <- suppressWarnings(stats::qlogis(p, 
                                                         location = res.mat$logis[1], 
                                                         scale = res.mat$logis[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::plogis(xx, 
                                          location = res.mat$logis[1], 
                                          scale = res.mat$logis[2]) ~ xx, 
                            col = logis.color, lwd = 2)
            leg.txt <- c(leg.txt, "Logistic")
            leg.col <- c(leg.col, logis.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting t distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$t[1])) {
        test.quantiles <- suppressWarnings(stats::qt(p, df = res.mat$t[1]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pt(xx, df = res.mat$t[1]) ~ xx, 
                            col = t.color, lwd = 2)
            leg.txt <- c(leg.txt, "Student")
            leg.col <- c(leg.col, t.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting exp ditribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$exp[1])) {
        test.quantiles <- suppressWarnings(stats::qexp(p, rate = res.mat$exp[1]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pexp(xx, rate = res.mat$exp[1]) ~ xx, 
                            col = exp.color, lwd = 2)
            leg.txt <- c(leg.txt, "Exponential")
            leg.col <- c(leg.col, exp.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting F distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$f[1])) {
        test.quantiles <- suppressWarnings(stats::qf(p, 
                                                     df1 = res.mat$f[1], 
                                                     df2 = res.mat$f[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pf(xx, 
                                      df1 = res.mat$f[1], 
                                      df2 = res.mat$f[2]) ~ xx, 
                            col = f.color, lwd = 2)
            leg.txt <- c(leg.txt, "F")
            leg.col <- c(leg.col, f.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting gamma distibution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$gamma[1])) {
        test.quantiles <- suppressWarnings(stats::qgamma(p, 
                                                         shape = res.mat$gamma[1], 
                                                         rate = res.mat$gamma[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pgamma(xx, 
                                          shape = res.mat$gamma[1], 
                                          rate = res.mat$gamma[2]) ~ xx, 
                            col = gamma.color, lwd = 2)
            leg.txt <- c(leg.txt, "Gamma")
            leg.col <- c(leg.col, gamma.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting lnorm ditribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$lnorm[1])) {
        test.quantiles <- suppressWarnings(stats::qlnorm(p, 
                                                         meanlog = res.mat$lnorm[1], 
                                                         sdlog = res.mat$lnorm[2] + 1/10000))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::plnorm(xx, 
                                          meanlog = res.mat$lnorm[1], 
                                          sdlog = res.mat$lnorm[2] + 1/10000) ~ xx, 
                            col = lnorm.color, lwd = 2)
            leg.txt <- c(leg.txt, "Lognormal")
            leg.col <- c(leg.col, lnorm.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting Weibull distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$weibull[1])) {
        test.quantiles <- suppressWarnings(stats::qweibull(p, 
                                                           shape = res.mat$weibull[1], 
                                                           scale = res.mat$weibull[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::pweibull(xx, 
                                            shape = res.mat$weibull[1], 
                                            scale = res.mat$weibull[2]) ~ xx, 
                            col = weibull.color, lwd = 2)
            leg.txt <- c(leg.txt, "Weibull")
            leg.col <- c(leg.col, weibull.color)
        }, silent = TRUE)
    }
    #-----------------------------------------------------------------------------
    # fitting unif distribution
    #-----------------------------------------------------------------------------
    if (!is.na(res.mat$unif[1])) {
        test.quantiles <- suppressWarnings(stats::qunif(p, 
                                                        min = res.mat$unif[1], 
                                                        max = res.mat$unif[2]))
        try(if (sum(abs(test.quantiles - q)) < tolPlot) {
            graphics::lines(stats::punif(xx, 
                                         min = res.mat$unif[1], 
                                         max = res.mat$unif[2]) ~ xx, 
                            col = unif.color, lwd = 2)
            leg.txt <- c(leg.txt, "Uniform")
            leg.col <- c(leg.col, unif.color)
        }, silent = TRUE)
    }
    
    #-----------------------------------------------------------------------------
    # creating legend
    #-----------------------------------------------------------------------------
    if (length(leg.txt) > 0) {
        graphics::legend("bottomright", legend = leg.txt, col = leg.col, lty = 1, bty = "n", lwd = 2)
    } else {
        # change by LG: Now, the tkmessageBox is produced in method fit.perc(): there the call of 
        # plotDiagnostics.perc() is checked if the message produced below is thrown. Given that case 
        # the message is caught and the tkmessageBox with this message produced:
        # tcltk::tkmessageBox(message = "Neither of the fitted distributions satisfies the tolerance constraint for plotting diagnostics!", icon = "error")
        message("Warning: Neither of the fitted distributions satisfies the tolerance constraint for plotting diagnostics!\n")
    }
} # end of plotDiagnostics.perc()


################################################################################

#' This function fits the amount of distribution families to given quantiles and returns
#' diagnostics that allow user to choose a most appropriate probability.
#'
#' Both inputs \code{p} and \code{q} should be of the same length. The items of
#' the probability vector \code{p} should lie between 0 and 1.
#'
#' @name rriskFitdist.perc
#' @aliases rriskFitdist.perc
#' @title Fitting an amount of distribution families by given quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage rriskFitdist.perc(p = c(0.025, 0.5, 0.975), q = c(9.68, 29.20, 50.98),
#'    show.output = TRUE, tolConv = 0.001, fit.weights = rep(1, length(p)))
#' @param p numerical vector giving probabilities.
#' @param q numerical vector giving percentiles.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param tolConv positive numerical value, the absolute convergence tolerance for reaching zero by fitting distributions
#' \code{get.norm.par} will be shown.
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @return Returns a list containing the data frame with the input vectors \code{p}
#' and \code{q} and the results matrix giving fitted distributions, estimated
#' parameters and a vector of theoretical percentiles calculated based on the
#' estimated parameters. If the consistency check of input parameters fails
#' the function returns \code{NA}.
#' @keywords fitdistrplus
#' @export
#' @importFrom mc2d qtriang
#' @importFrom mc2d qpert
#' @importFrom eha qgompertz
#' @importFrom msm qtnorm
#' @examples
#' fit.results1 <- rriskFitdist.perc(show.output = FALSE)
#' fit.results1
#'
#' fit.results2 <- rriskFitdist.perc(show.output = FALSE, tolConv = 0.6)
#' fit.results2
#'
#' p <- c(0.2, 0.7)
#' q <- c(2.6, 19.1)
#' fit.results3 <- rriskFitdist.perc(p = p, q = q, show.output = FALSE)
#' fit.results3
#'
#' p <- c(0.3, 0.8, 0.9)
#' q <- c(10, 20, 40)
#' fit.results4 <- rriskFitdist.perc(p = p, q = q, show.output = FALSE)
#' fit.results4
#'
#' ## Example with fitted pert distribution
#' p <- c(0.025, 0.5, 0.6, 0.975)
#' q <- mc2d::qpert(p = p, min = 0, mode = 3, max = 10, shape = 5)
#' fit.results5 <- rriskFitdist.perc(p = p, q = q, show.output = FALSE)
#' fit.results5
#' 
rriskFitdist.perc <- function(p = c(0.025, 0.5, 0.975), 
                              q = c(9.68, 29.20, 50.98), 
                              show.output = TRUE, 
                              tolConv = 0.001, 
                              fit.weights = rep(1, length(p))) {
    #-----------------------------------------------------------------------------
    # check general consistency of the input data
    #-----------------------------------------------------------------------------
    if (length(p) != length(q)) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vectors of probabilities, percentiles or/and weights are not of the same length!", call. = FALSE)
    }
    if (length(p) == 0 | length(q) == 0) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, either the vector of probabilities or the vector of quantiles is empty!", call. = FALSE)
    }
    if (length(fit.weights) == 0) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of weights is empty!", call. = FALSE)
    }
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, one of the following vectors is not numeric: probabilities, percentiles, weights!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (any(fit.weights <= 0)) {
        on.exit(return(invisible(NA)))
        stop("INVALID INPUT, all items of the argument 'fit.weights' should be positive!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # defining help variables
    #-----------------------------------------------------------------------------
    Perc <- c(p, 0.0001, 0.001, 0.01, 0.025, 0.05, 0.1, 0.5, 0.9, 0.95, 0.975, 0.99, 0.999, 0.9999)
    Quantiles <- c(q, rep(NA, 13))
    Quantiles <- Quantiles[order(Perc)]
    Perc <- Perc[order(Perc)]
    Quantiles <- Quantiles[!duplicated(Perc)]
    Perc <- Perc[!duplicated(Perc)]
    Perc <- round(Perc, digits = 4)
    res.mat <- data.frame(weight = rep(0, length(Perc) + 4), 
                          Quantiles = c(rep(NA, 4), Quantiles),
                          norm = NA, beta = NA, cauchy = NA, logis = NA, 
                          t = NA, chisq = NA, chisqnc = NA, exp = NA, 
                          f = NA, gamma = NA, lnorm = NA, unif = NA, 
                          weibull = NA, triang = NA, gompertz = NA, 
                          pert = NA, tnorm = NA)
    rownames(res.mat) <- c(paste("Para", 1:4, sep = ""), Perc * 100)
    res.mat[as.character(p * 100), "weight"] <- fit.weights
    
    message("\nBegin fitting distributions ---------------------------------------")
    
    ## tnorm
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.tnorm.par(p = p, q = q, 
                          show.output = show.output, 
                          plot = FALSE, 
                          tol = tolConv, 
                          fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$tnorm[1:4] <- parameters
        res.mat$tnorm[-c(1:4)] <- msm::qtnorm(p = Perc, 
                                              mean = parameters["mean"], 
                                              sd = parameters["sd"], 
                                              lower = parameters["lower"], 
                                              upper = parameters["upper"])
    }
    message("* fitting truncated normal distribution ... ", res)
    
    ## chisqnc
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.chisqnc.par(p = p, q = q, 
                            show.output = show.output, 
                            plot = FALSE, 
                            tol = tolConv, 
                            fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$chisqnc[1:2] <- parameters
        res.mat$chisqnc[-c(1:4)] <- stats::qchisq(p = Perc, 
                                                  df = parameters["df"], 
                                                  ncp = parameters["ncp"])
    } 
    message("* fitting non-central chi-square distribution ... ", res)
    
    ## pert
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.pert.par(p = p, q = q, 
                         show.output = show.output, 
                         plot = FALSE, 
                         tol = tolConv, 
                         fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$pert[1:4] <- parameters
        res.mat$pert[-c(1:4)] <- mc2d::qpert(p = Perc, 
                                             min = parameters["min"], 
                                             mode = parameters["mode"], 
                                             max = parameters["max"], 
                                             shape = parameters["shape"])
    }
    message("* fitting PERT distribution ... ", res)
    
    ## triang
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.triang.par(p = p, q = q, 
                           show.output = show.output, 
                           plot = FALSE, 
                           tol = tolConv, 
                           fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$triang[1:3] <- parameters
        res.mat$triang[-c(1:4)] <- mc2d::qtriang(p = Perc, 
                                                 min = parameters["min"], 
                                                 mode = parameters["mode"], 
                                                 max = parameters["max"])
    }
    message("* fitting triangular distribution ... ", res)
    
    ## gompertz
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.gompertz.par(p = p, q = q, 
                             show.output = show.output, 
                             plot = FALSE, 
                             tol = tolConv, 
                             fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$gompertz[1:2] <- parameters
        res.mat$gompertz[-c(1:4)] <- eha::qgompertz(p = Perc, 
                                                    shape = parameters["shape"], 
                                                    scale = parameters["scale"])
    }
    message("* fitting Gompertz distribution ... ", res)
    
    ## normal
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.norm.par(p = p, q = q, 
                         show.output = show.output, 
                         plot = FALSE, 
                         tol = tolConv, 
                         fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$norm[1:2] <- parameters
        res.mat$norm[-c(1:4)] <- stats::qnorm(p = Perc, 
                                              mean = parameters["mean"], 
                                              sd = parameters["sd"])
    }
    message("* fitting normal distribution ... ", res)
    
    ## beta
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.beta.par(p = p, q = q, 
                         show.output = show.output, 
                         plot = FALSE, 
                         tol = tolConv, 
                         fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$beta[1:2] <- parameters
        res.mat$beta[-c(1:4)] <- stats::qbeta(p = Perc, 
                                              shape1 = parameters["shape1"], 
                                              shape2 = parameters["shape2"])
    }
    message("* fitting beta distribution ... ", res)
    
    ## cauchy
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.cauchy.par(p = p, q = q, 
                           show.output = show.output, 
                           plot = FALSE, 
                           tol = tolConv, 
                           fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$cauchy[1:2] <- parameters
        res.mat$cauchy[-c(1:4)] <- stats::qcauchy(p = Perc, 
                                                  location = parameters["location"], 
                                                  scale = parameters["scale"])
    }
    message("* fitting Cauchy distribution ... ", res)
    
    ## chisq
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.chisq.par(p = p, q = q, 
                          show.output = show.output, 
                          plot = FALSE, 
                          tol = tolConv, 
                          fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$chisq[1] <- parameters
        res.mat$chisq[-c(1:4)] <- stats::qchisq(p = Perc, 
                                                df = parameters["df"])
    }
    message("* fitting chi-square distribution ... ", res)
    
    ## logis
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.logis.par(p = p, q = q, 
                          show.output = show.output, 
                          plot = FALSE, 
                          tol = tolConv, 
                          fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$logis[1:2] <- parameters
        res.mat$logis[-c(1:4)] <- stats::qlogis(p = Perc, 
                                                location = parameters["location"], 
                                                scale = parameters["scale"])
    }
    message("* fitting logistic distribution ... ", res)
    
    ## t
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.t.par(p = p, q = q, 
                      show.output = show.output, 
                      plot = FALSE, 
                      tol = tolConv, 
                      fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$t[1] <- parameters
        res.mat$t[-c(1:4)] <- stats::qt(p = Perc, 
                                        df = parameters["df"])
    }
    message("* fitting Student's t-distribution ... ", res)
    
    ## exp
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.exp.par(p = p, q = q, 
                        show.output = show.output, 
                        plot = FALSE, 
                        tol = tolConv, 
                        fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$exp[1] <- parameters
        res.mat$exp[-c(1:4)] <- stats::qexp(p = Perc, 
                                            rate = parameters["rate"])
    }
    message("* fitting exponential distribution ... ", res)
    
    ## F
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.f.par(p = p, q = q, 
                      show.output = show.output, 
                      plot = FALSE, 
                      tol = tolConv, 
                      fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$f[1:2] <- parameters
        res.mat$f[-c(1:4)] <- stats::qf(p = Perc, 
                                        df1 = parameters["df1"], 
                                        df2 = parameters["df2"])
    }
    message("* fitting F-distribution ... ", res)
    
    ## gamma
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.gamma.par(p = p, q = q, 
                          show.output = show.output, 
                          plot = FALSE, 
                          tol = tolConv, 
                          fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$gamma[1:2] <- parameters
        res.mat$gamma[-c(1:4)] <- stats::qgamma(p = Perc, 
                                                shape = parameters["shape"], 
                                                rate = parameters["rate"])
    }
    message("* fitting gamma distribution ... ", res)
    
    ## Weibull
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.weibull.par(p = p, q = q, 
                            show.output = show.output, 
                            plot = FALSE, 
                            tol = tolConv, 
                            fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$weibull[1:2] <- parameters
        res.mat$weibull[-c(1:4)] <- stats::qweibull(p = Perc, 
                                                    shape = parameters["shape"], 
                                                    scale = parameters["scale"])
    }
    message("* fitting Weibull distribution ... ", res)
    
    ## lognormal
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.lnorm.par(p = p, q = q, 
                          show.output = show.output, 
                          plot = FALSE, 
                          tol = tolConv, 
                          fit.weights = fit.weights)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$lnorm[1:2] <- parameters
        res.mat$lnorm[-c(1:4)] <- stats::qlnorm(p = Perc, 
                                                meanlog = parameters["meanlog"], 
                                                sdlog = parameters["sdlog"])
    }
    message("* fitting lognormal distribution ... ", res)
    
    ## uniform
    parameters <- NA
    try({
        parameters <- suppressMessages(suppressWarnings({
            get.unif.par(p = p, q = q, plot = FALSE)
        }))},
        silent = TRUE
    )
    res <- ifelse(any(is.na(parameters)), "failed", "OK")
    if (res == "OK") {
        res.mat$unif[1:2] <- parameters
        res.mat$unif[-c(1:4)] <- stats::qunif(p = Perc, 
                                              min = parameters["min"], 
                                              max = parameters["max"])
    }
    message("* fitting uniform distribution ... ", res)
    
    message("End fitting distributions -----------------------------------------\n")
    
    if (all(is.na(res.mat[1:4, -1]))) {
        if (is.element("package:rrisk", search())) { # wenn "rrisk" vorhanden, mache weiter. sonst breche ab.
            #on.exit(.generate.newitem())
            stop("\n\nSorry, no pdfs found that match with the specified percentiles\n", call. = FALSE)
        } else {
            on.exit(return(invisible(NA)))
            stop("\n\nSorry, no pdfs found that match with the specified percentiles\n", call. = FALSE)
        }
    }
    res.mat <- data.frame(apply(res.mat, c(1, 2), function(x) round(x, digits = 2)))
    output <- list(data.frame(p, q), res.mat)
    names(output) <- c("p/q", "results")
    return(output)
} # end of rriskFitdist.perc()


################################################################################
#' \code{get.beta.par} returns the parameters of estimated beta distribution where the
#' \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the numbers of weightings
#' must be identical and should be at least two. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.beta.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{CG} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.beta.par
#' @aliases get.beta.par
#' @title Fitting parameters of a Beta distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.beta.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE, plot = TRUE,
#'      tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a Beta distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the eror massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good til very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pbeta} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qbeta(p = c(0.025, 0.5, 0.975), shape1 = 2, shape2 = 3)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.beta.par(q = q)
#' get.beta.par(q = q, scaleX = c(0.001, 0.999))
#' get.beta.par(q = q, fit.weights = c(10, 1, 10))
#' get.beta.par(q = q, fit.weights = c(1, 10, 1))
#' get.beta.par(q = q, fit.weights = c(100, 1, 100))
#' get.beta.par(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qbeta(p = c(0.025, 0.5, 0.975), shape1 = 1, shape2 = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.beta.par(q = q)
#' get.beta.par(q = q, fit.weights = c(10, 1, 10))
#' get.beta.par(q = q, fit.weights = c(1, 10, 1))
#' get.beta.par(q = q, fit.weights = c(100, 1, 100))
#' get.beta.par(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qbeta(p = c(0.025, 0.5, 0.975), shape1 = 0.3, shape2 = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.beta.par(q = q)
#' get.beta.par(q = q, fit.weights = c(10, 1, 10))
#' get.beta.par(q = q, fit.weights = c(1, 10, 1))
#' get.beta.par(q = q, fit.weights = c(100, 1, 100))
#' get.beta.par(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qbeta(p = c(0.025, 0.975), shape1 = 2, shape2 = 3)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.beta.par(p = c(0.025, 0.975), q = q)
#' get.beta.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' get.beta.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 10))
#' get.beta.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1))
#' get.beta.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 100))
#' graphics::par(old.par)
#' 
get.beta.par <- function(p = c(0.025, 0.5, 0.975), q, 
                         show.output = TRUE, plot = TRUE, 
                         tol = 0.001, fit.weights = rep(1, length(p)), 
                         scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) < 0 | max(q) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain (0, 1) => beta distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights / sum(fit.weights)
    minimize <- function(shape) {
        summand <- suppressWarnings(stats::pbeta(q = q, 
                                                 shape1 = shape[1], 
                                                 shape2 = shape[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(0.1, 0.1), 
                            minimize, method = "L-BFGS-B", 
                            lower = 0.001, upper = 10000), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking the output
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(minimize, 
                                method = "CG"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'CG' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'CG' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("shape1", "shape2")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("shape1", "shape2")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # creating graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("shape1 = ", round(Par["shape1"], digits = 2))
        main2 <- paste("shape2 = ", round(Par["shape2"], digits = 2))
        main <- paste("Beta (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qbeta(p = min(p) * scaleX[1], 
                                      shape1 = Par["shape1"], 
                                      shape2 = Par["shape2"]),
                         stats::qbeta(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                      shape1 = Par["shape1"], 
                                      shape2 = Par["shape2"]))
        #Support <- seq(Support.lim[1], Support.lim[2], length = 200)
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pbeta(Support, Par["shape1"], Par["shape2"])
        graphics::plot(Support, Probability, 
                       type = "l", xlim = range(Support.lim, q), 
                       main = main, xlab = "Quantiles", 
                       sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
} # end of get.beta.par()



################################################################################
#' \code{get.cauchy.par} returns the parameters of a Cauchy distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least two. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the
#' 97.5th percentile, respectively. \code{get.cauchy.par} uses the R function
#' \code{optim} with the method \code{L-BFGS-B}.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.cauchy.par
#' @aliases get.cauchy.par
#' @title Fitting parameters of a Cauchy distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.cauchy.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'    plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a Cauchy distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the eror massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good til very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pcauchy} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qcauchy(p = c(0.025, 0.5, 0.975), location = 0, scale = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.cauchy.par(q = q)
#' get.cauchy.par(q = q, scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.5, 0.5))
#' graphics::par(old.par)
#'
#' q <- stats::qcauchy(p = c(0.025, 0.5, 0.975), location = 3, scale = 5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.cauchy.par(q = q, scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.5, 0.5))
#' graphics::par(old.par)
#'
#' q <- stats::qcauchy(p = c(0.025, 0.5, 0.975), location = 0.1, scale = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.cauchy.par(q = q, scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.5, 0.5))
#' get.cauchy.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.5, 0.5))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qcauchy(p = c(0.025, 0.975), location = 0.1, scale = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.cauchy.par(p = c(0.025, 0.975), q = q, scaleX = c(0.5, 0.5))
#' get.cauchy.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1), scaleX = c(0.5, 0.5))
#' get.cauchy.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1), scaleX = c(0.5, 0.5))
#' get.cauchy.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 10), scaleX = c(0.5, 0.5))
#' get.cauchy.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 100), scaleX = c(0.5, 0.5))
#' graphics::par(old.par)
#' 
get.cauchy.par <- function(p = c(0.025, 0.5, 0.975), q, 
                           show.output = TRUE, plot = TRUE, 
                           tol = 0.001, fit.weights = rep(1, length(p)), 
                           scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights / sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::pcauchy(q = q, 
                                                   location = theta[1], 
                                                   scale = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(1, 1), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(-10000, 0.001), 
                            upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(1, 1),  # formerly 'par = start' which was not defined!
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("location", "scale")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("location", "scale")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plot graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("location = ", round(Par["location"], digits = 2))
        main2 <- paste("scale = ", round(Par["scale"], digits = 2))
        main <- paste("Cauchy (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qcauchy(p = min(p) * scaleX[1], 
                                        location = Par["location"], 
                                        scale = Par["scale"]),
                         stats::qcauchy(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                        location = Par["location"], 
                                        scale = Par["scale"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pcauchy(Support, Par["location"], Par["scale"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.chisq.par} returns the parameters of a chi-square distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least one. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the
#' 97.5th percentile, respectively. \code{get.chisq.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted, if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.chisq.par
#' @aliases get.chisq.par
#' @title Fitting parameter of a chi-square distribution from one or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.chisq.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'    plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1)).
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE})
#' @return Returns fitted parameters of a chi-square distribution or missing
#' values (\code{NA}'s), if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pchisq} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qchisq(p = c(0.025, 0.5, 0.975), df = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.chisq.par(q = q)
#' get.chisq.par(q = q, fit.weights = c(10, 1, 10))
#' get.chisq.par(q = q, fit.weights = c(100, 1, 100))
#' get.chisq.par(q = q, fit.weights = c(1, 10, 1))
#' get.chisq.par(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qchisq(p = c(0.025, 0.5, 0.975), df = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.chisq.par(q = q, scaleX = c(0.1, 0.1))
#' get.chisq.par(q = q, fit.weights = c(10, 1, 10))
#' get.chisq.par(q = q, fit.weights = c(100, 1, 100))
#' get.chisq.par(q = q, fit.weights = c(1, 10, 1))
#' get.chisq.par(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qchisq(p = c(0.025, 0.5, 0.975), df = 20)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.chisq.par(q = q)
#' get.chisq.par(q = q, fit.weights = c(10, 1, 10))
#' get.chisq.par(q = q, fit.weights = c(100, 1, 100))
#' get.chisq.par(q = q, fit.weights =c(1, 10, 1))
#' get.chisq.par(q = q, fit.weights =c(1, 100, 1))
#' graphics::par(old.par)
#'
#' ## example with only one quantile
#' q <- stats::qchisq(p = c(0.025), df = 20)
#' old.par <- graphics::par(mfrow = c(1, 3))
#' get.chisq.par(p = c(0.025), q = q)
#' get.chisq.par(p = c(0.025), q = q, fit.weights = 10)
#' get.chisq.par(p = c(0.025), q = q, fit.weights = 100)
#' graphics::par(old.par)
#'
get.chisq.par <- function(p = c(0.025, 0.5, 0.975), q, 
                          show.output = TRUE, plot = TRUE, 
                          tol = 0.001, fit.weights = rep(1, length(p)), 
                          scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, The vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, Items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain [0, inf) => chi-square distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least one quantile must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights / sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::pchisq(q = q, df = theta) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    suppressWarnings(lm <- stats::lm(q ~ p)    )
    suppressWarnings(m <- stats::predict(lm, newdata = list(p = 0.5))[[1]])
    # using the approximation of the chisq median
    start <- ((3 * m^3 + 6 * m^2 + 2 * m)/81 + 2 * m * sqrt(2 * m + 3)/(81 * sqrt(3)))^(1 / 3) + (3 * m^2 + 4 * m)/(27 * ((3 * m^3 + 6 * m^2 + 2 * m)/81 + 2 * m * sqrt(2 * m + 3)/(81 * sqrt(3)))^(1/3)) + (3 * m + 2)/9
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = start, 
                            minimize, method = "L-BFGS-B", 
                            lower = 0.001, upper = 10000), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = start,
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("df")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("df")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plot graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("df = ", round(Par["df"], digits = 2))
        main <- paste("Chi-square (", main1, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qchisq(p = min(p) * scaleX[1], df = Par["df"]), 
                         stats::qchisq(p = (max(p) + (1 - max(p)) * scaleX[2]), df = Par["df"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pchisq(Support, Par["df"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.chisqnc.par} returns the parameters of a non-central chi-square
#' distribution where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weights
#' must be identical and should be at least one. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.chisqnc.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.chisqnc.par
#' @aliases get.chisqnc.par
#' @title Fitting parameters of a non-central chi-square distribution from one or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.chisqnc.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'    plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a non-central chi-square distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the eror massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pchisq} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qchisq(p = c(0.025, 0.5, 0.975), df = 2, ncp = 4)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.chisqnc.par(q = q)
#' get.chisqnc.par(q = q, scaleX = c(0.1, 0.9999999))
#' get.chisqnc.par(q = q, fit.weights = c(100, 1, 100))
#' get.chisqnc.par(q = q, fit.weights = c(10, 1, 10))
#' get.chisqnc.par(q = q, fit.weights = c(1, 100, 1))
#' get.chisqnc.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qchisq(p = c(0.025, 0.5, 0.975), df = 0.1, ncp = 0.4)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.chisqnc.par(q = q)
#' get.chisqnc.par(q = q, fit.weights = c(100, 1, 100))
#' get.chisqnc.par(q = q, fit.weights = c(10, 1, 10))
#' get.chisqnc.par(q = q, fit.weights = c(1, 100, 1))
#' get.chisqnc.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qchisq(p = c(0.025, 0.5, 0.975), df = 1, ncp = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.chisqnc.par(q = q)
#' get.chisqnc.par(q = q, fit.weights = c(100, 1, 100))
#' get.chisqnc.par(q = q, fit.weights = c(10, 1, 10))
#' get.chisqnc.par(q = q, fit.weights = c(1, 100, 1))
#' get.chisqnc.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantile
#' q <- stats::qchisq(p = c(0.025, 0.95), df = 20, ncp = 20)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.chisqnc.par(p = c(0.025, 0.975), q = q)
#' get.chisqnc.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1))
#' get.chisqnc.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 100))
#' get.chisqnc.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' get.chisqnc.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 10))
#' graphics::par(old.par)
#'
get.chisqnc.par <- function(p = c(0.025, 0.5, 0.975), q,
                            show.output = TRUE, plot = TRUE, 
                            tol = 0.001, fit.weights = rep(1, length(p)), 
                            scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, The vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, Items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain [0, inf) => chi-square distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least one quantile must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::pchisq(q = q, 
                                                  df = theta[1], 
                                                  ncp = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(1, 1), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(0.001, 0.00001), 
                            upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(1, 1),
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("df", "ncp")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("df", "ncp")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plot graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("df = ", round(Par["df"], digits = 2))
        main2 <- paste("ncp = ", round(Par["ncp"], digits = 2))
        main <- paste("non-central chi-square (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qchisq(p = min(p) * scaleX[1], df = Par["df"]),
                         stats::qchisq(p = (max(p) + (1 - max(p)) * scaleX[2]), df = Par["df"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pchisq(Support, Par["df"], Par["ncp"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.exp.par} returns the parameters of an exponential distribution
#' where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least one. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the
#' 97.5th percentile, respectively. \code{get.exp.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method \code{BFGS}
#' will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.exp.par
#' @aliases get.exp.par
#' @title Fitting parameters of an exponential distribution from one or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR),  \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.exp.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of an exponential distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the eror massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good til very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pexp} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qexp(p = c(0.025, 0.5, 0.975), rate = 2)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.exp.par(q = q)
#' get.exp.par(q = q, fit.weights = c(100, 1, 100))
#' get.exp.par(q = q, fit.weights = c(10, 1, 10))
#' get.exp.par(q = q, fit.weights = c(1, 100, 1))
#' get.exp.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qexp(p = c(0.025, 0.5, 0.975), rate = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.exp.par(q = q)
#' get.exp.par(q = q, fit.weights = c(100, 1, 100))
#' get.exp.par(q = q, fit.weights = c(10, 1, 10))
#' get.exp.par(q = q, fit.weights = c(1, 100, 1))
#' get.exp.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qexp(p = c(0.025, 0.5, 0.975), rate = 0.001)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.exp.par(q = q)
#' get.exp.par(q = q, fit.weights = c(100, 1, 100))
#' get.exp.par(q = q, fit.weights = c(10, 1, 10))
#' get.exp.par(q = q, fit.weights = c(1, 100, 1))
#' get.exp.par(q = q, tol = 0.2, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qexp(p = c(0.025, 0.5, 0.975), rate = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.exp.par(q = q)
#' get.exp.par(q = q, fit.weights = c(100, 1, 100))
#' get.exp.par(q = q, fit.weights = c(10, 1, 10))
#' get.exp.par(q = q, fit.weights = c(1, 100, 1))
#' get.exp.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## example with only one quantile
#' q <- stats::qexp(p = c(0.025), rate = 2)
#' old.par <- graphics::par(mfrow = c(1, 3))
#' get.exp.par(p = c(0.025), q = q)
#' get.exp.par(p = c(0.025), q = q, fit.weights = 10)
#' get.exp.par(p = c(0.025), q = q, fit.weights = 100)
#' graphics::par(old.par)
#'
get.exp.par <- function(p = c(0.025, 0.50,.975), q,
                        show.output = TRUE, plot = TRUE, 
                        tol = 0.001, fit.weights = rep(1, length(p)), 
                        scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(1) < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain [0, inf) => exponential distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least one quantile must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizinmg procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(par) {
        summand <- suppressWarnings(stats::pexp(q = q, rate = par) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    Start <- 1/mean(q)
    try1 <- try(
        fit <- stats::optim(par = Start, 
                            minimize, method = "L-BFGS-B", 
                            lower = 0.001, upper = 10000), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = Start,
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("rate")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("rate")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("rate = ", round(Par["rate"], digits = 2))
        main <- paste("Exponential (", main1, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qexp(p = min(p) * scaleX[1], 
                                     rate = Par["rate"]),
                         stats::qexp(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                     rate = Par["rate"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pexp(Support, Par["rate"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), 
                       main = main, xlab = "Quantiles", 
                       sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}



################################################################################
#' \code{get.f.par} returns the parameters of a F distribution where the
#' \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least two. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.f.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{BFGS} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.f.par
#' @aliases get.f.par
#' @title Fitting parameters of a F distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.f.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'    plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a F distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{stats::pf} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qf(p = c(0.025, 0.5, 0.975), df1 = 2, df2 = 10)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.f.par(q = q, scaleX = c(0.1, 0.5))
#' get.f.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.1, 0.5))
#' get.f.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.1, 0.5))
#' get.f.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.1, 0.5))
#' get.f.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.1, 0.5))
#' graphics::par(old.par)
#'
#' q <- stats::qf(p = c(0.025, 0.5, 0.975), df1 = 0.2, df2 = 0.3)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.f.par(q = q, scaleX = c(0.1, 0.2))
#' get.f.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.1, 0.999))
#' get.f.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.1, 0.2))
#' get.f.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.1, 0.9999))
#' get.f.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.1, 0.9999))
#' graphics::par(old.par)
#'
#' q <- stats::qf(p = c(0.025, 0.5, 0.975), df1 = 1, df2 = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.f.par(q = q, scaleX = c(0.1, 0.2))
#' get.f.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.1, 0.2))
#' get.f.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.1, 0.2))
#' get.f.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.1, 0.2))
#' get.f.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.1, 0.2))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qf(p = c(0.025, 0.975), df1 = 2, df2 = 3)
#' old.par <- graphics::par(mfrow = c(1, 3))
#' get.f.par(p = c(0.025, 0.975), q = q)
#' get.f.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1))
#' get.f.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' graphics::par(old.par)
#'
get.f.par <- function(p = c(0.025, 0.5, 0.975), q, 
                      show.output = TRUE, plot = TRUE,
                      tol = 0.001, fit.weights = rep(1, length(p)),
                      scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) <= 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain [0, inf) => F distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::pf(q = q, 
                                              df1 = theta[1], 
                                              df2 = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(1, 1), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(0.001, 0.001), upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(1, 1),
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("df1", "df2")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("df1", "df2")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("df1 = ", round(Par["df1"], digits = 2))
        main2 <- paste("df2 = ", round(Par["df2"], digits = 2))
        main <- paste("F (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qf(p = min(p) * scaleX[1], 
                                   df1 = Par["df1"], 
                                   df2 = Par["df2"]),
                         stats::qf(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                   df1 = Par["df1"], 
                                   df2 = Par["df2"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pf(Support, Par["df1"], Par["df2"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), 
                       main = main, xlab = "Quantiles", 
                       sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.gamma.par} returns the parameters of a gamma distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least two. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.gamma.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{BFGS} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.gamma.par
#' @aliases get.gamma.par
#' @title Fitting parameters of a gamma distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting),  \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.gamma.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1)).
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a gamma distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pgamma} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qgamma(p = c(0.025, 0.5, 0.975), shape = 10, rate = 10)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.gamma.par(q = q)
#' get.gamma.par(q = q, scaleX = c(0.00001, 0.9999))
#' get.gamma.par(q = q, fit.weights = c(100, 1, 100))
#' get.gamma.par(q = q, fit.weights = c(10, 1, 10))
#' get.gamma.par(q = q, fit.weights = c(1, 100, 1))
#' get.gamma.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qgamma(p = c(0.025, 0.5, 0.975), shape = 0.1, rate = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.gamma.par(q = q)
#' get.gamma.par(q = q, fit.weights = c(100, 1, 100))
#' get.gamma.par(q = q, fit.weights = c(10, 1, 10))
#' get.gamma.par(q = q, fit.weights = c(1, 100, 1))
#' get.gamma.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qgamma(p = c(0.025, 0.5, 0.975), shape = 1, rate = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.gamma.par(q = q)
#' get.gamma.par(q = q, fit.weights = c(100, 1, 100))
#' get.gamma.par(q = q, fit.weights = c(10, 1, 10))
#' get.gamma.par(q = q, fit.weights = c(1, 100, 1))
#' get.gamma.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qgamma(p = c(0.025, 0.975), shape = 10, rate = 10)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.gamma.par(p = c(0.025, 0.975), q = q)
#' get.gamma.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1))
#' get.gamma.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 100))
#' get.gamma.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' get.gamma.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 10))
#' graphics::par(old.par)
#'
get.gamma.par <- function(p = c(0.025, 0.5, 0.975), q, 
                          show.output = TRUE, plot = TRUE, 
                          tol = 0.001, fit.weights = rep(1, length(p)), 
                          scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain [0, inf) => Gamma distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::pgamma(q = q, 
                                                  shape = theta[1], 
                                                  rate = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(1, 1), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(0.001, 0.001), 
                            upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(1, 1),
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("shape", "rate")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("shape", "rate")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("shape = ", round(Par["shape"], digits = 2), sep = "")
        main2 <- paste("rate = ", round(Par["rate"], digits = 2), sep = "")
        main <- paste("Gamma (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qgamma(p = min(p) * scaleX[1], 
                                       shape = Par["shape"], 
                                       rate = Par["rate"]),
                         stats::qgamma(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                       shape = Par["shape"], 
                                       rate = Par["rate"]))
        Support <- seq(min(min(q), Support.lim[1]), max(max(q), Support.lim[2]), length = 200)
        Probability <- stats::pgamma(Support, Par["shape"], Par["rate"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.gompertz.par} returns the parameters of a Gompertz distribution
#' where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least two. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.gompertz.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{BFGS} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.gompertz.par
#' @aliases get.gompertz.par
#' @title Fitting parameters of a Gompertz distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.gompertz.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a Gompertz distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note Comply with a parametrization of this distribution. The definition of this
#' distribution in the literature is not unique.
#' \cr \cr
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{\link[eha:pgompertz]{pgompertz}} from the package \pkg{eha} for distribution 
#' implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- eha::qgompertz(p = c(0.025, 0.5, 0.975), shape = 2, scale = 5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.gompertz.par(q = q)
#' get.gompertz.par(q = q, fit.weights = c(100, 1, 100))
#' get.gompertz.par(q = q, fit.weights = c(10, 1, 10))
#' get.gompertz.par(q = q, fit.weights = c(1, 100, 1))
#' get.gompertz.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- eha::qgompertz(p = c(0.025, 0.5, 0.975), shape = 0.2, scale = 0.5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.gompertz.par(q = q)
#' get.gompertz.par(q = q, fit.weights = c(100, 1, 100))
#' get.gompertz.par(q = q, fit.weights = c(10, 1, 10))
#' get.gompertz.par(q = q, fit.weights = c(1, 100, 1))
#' get.gompertz.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- eha::qgompertz(p = c(0.025, 0.5, 0.975), shape = 1, scale = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.gompertz.par(q = q)
#' get.gompertz.par(q = q, fit.weights = c(100, 1, 100))
#' get.gompertz.par(q = q, fit.weights = c(10, 1, 10))
#' get.gompertz.par(q = q, fit.weights = c(1, 100, 1))
#' get.gompertz.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- eha::qgompertz(p = c(0.025, 0.975), shape = 2, scale = 5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.gompertz.par(p = c(0.025, 0.975), q = q)
#' get.gompertz.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1), scaleX = c(0.0001, 0.9999))
#' get.gompertz.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 100))
#' get.gompertz.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' get.gompertz.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 10))
#' graphics::par(old.par)
#'
get.gompertz.par <- function(p = c(0.025, 0.5, 0.975), q,
                             show.output = TRUE, plot = TRUE, 
                             tol = 0.001, fit.weights = rep(1, length(p)), 
                             scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) <= 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain (0, inf) => Gompertz distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(eha::pgompertz(q = q, 
                                                   shape = theta[1], 
                                                   scale = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(1, 1), 
                            minimize, method = "L-BFGS-B",
                            lower = c(0.001, 0.001), upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(1, 1),
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("shape", "scale")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("shape", "scale")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plot graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("shape = ", round(Par["shape"], digits = 2))
        main2 <- paste("scale = ", round(Par["scale"], digits = 2))
        main <- paste("Gompertz (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(eha::qgompertz(p = min(p) * scaleX[1], 
                                        shape = Par["shape"], 
                                        scale = Par["scale"]),
                         eha::qgompertz(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                        shape = Par["shape"], 
                                        scale = Par["scale"]))
        Support <- seq(min(min(q), Support.lim[1]), max(max(q), Support.lim[2]), length = 200)
        Probability <- eha::pgompertz(Support, Par["shape"], Par["scale"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.hyper.par} returns the parameters of a hypergeometric distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least three. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.hyper.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.hyper.par
#' @aliases get.hyper.par
#' @title Fitting parameters of a hypergeometric  distribution from three or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.hyper.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a hypergeometric distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{phyper} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qhyper(p = c(0.025, 0.5, 0.975), m = 5, n = 3, k = 3)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.hyper.par(q = q)
#' get.hyper.par(q = q, tol = 1)
#' get.hyper.par(q = q, fit.weights = c(100, 1, 100))
#' get.hyper.par(q = q, fit.weights = c(10, 1, 10))
#' get.hyper.par(q = q, fit.weights = c(1, 100, 1))
#' get.hyper.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qhyper(p = c(0.025, 0.5, 0.975), m = 10, n = 5, k = 4)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.hyper.par(q = q)
#' get.hyper.par(q = q, fit.weights = c(100, 1, 100))
#' get.hyper.par(q = q, fit.weights = c(10, 1, 10))
#' get.hyper.par(q = q, fit.weights = c(1, 100, 1))
#' get.hyper.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
get.hyper.par <- function(p = c(0.025, 0.5, 0.975), q, 
                          show.output = TRUE, plot = TRUE, 
                          tol = 0.001, fit.weights = rep(1, length(p)), 
                          scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) < 0 | any(abs(q - trunc(q)) != 0)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles must be integer and positive => Hypergeometric distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 3) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least three quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::phyper(q = q, m = theta[1], n = theta[2], k = theta[3]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(9, 6, 7), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(0.001, 0.001, 0.001), 
                            upper = c(10000, 10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(9, 6, 7),
                                minimize, 
                                method = "SANN"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'SANN' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'SANN' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("m", "n", "k")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("m", "n", "k")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("m = ", round(Par["m"], digits = 2), sep = "")
        main2 <- paste("n = ", round(Par["n"], digits = 2), sep = "")
        main3 <- paste("k = ", round(Par["k"], digits = 2), sep = "")
        main <- paste("Hypergeo (", main1, ", ", main2, ", ", main3, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qhyper(p = min(p) * scaleX[1], 
                                       m = Par["m"], n = Par["n"], k = Par["k"]),
                         stats::qhyper(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                       m = Par["m"], n = Par["n"], k = Par["k"]))
        Support <- seq(min(min(q), Support.lim[1]), max(max(q), Support.lim[2]), length = 200)
        Probability <- stats::phyper(Support, Par["m"], Par["n"], Par["k"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}



################################################################################
#' \code{get.lnorm.par} returns the parameters of a lognormal distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least two. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.lnorm.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{Nelder-Mead} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.lnorm.par
#' @aliases get.lnorm.par
#' @title Fitting parameters of a lognormal distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.lnorm.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a lognormal distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note Comply with a parametrization of this distribution. The definition of this
#' distribution in the literature is not unique.
#' \cr \cr
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{plnorm} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qlnorm(p = c(0.025, 0.5, 0.975), meanlog = 4, sdlog = 0.8)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.lnorm.par(q = q)
#' get.lnorm.par(q = q, fit.weights = c(100, 1, 100))
#' get.lnorm.par(q = q, fit.weights = c(10, 1, 10))
#' get.lnorm.par(q = q, fit.weights = c(1, 100, 1))
#' get.lnorm.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qlnorm(p = c(0.025, 0.5, 0.975), meanlog=-4, sdlog = 0.8)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.lnorm.par(q = q)
#' get.lnorm.par(q = q, fit.weights = c(100, 1, 100))
#' get.lnorm.par(q = q, fit.weights = c(10, 1, 10))
#' get.lnorm.par(q = q, fit.weights = c(1, 100, 1))
#' get.lnorm.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qlnorm(p = c(0.025, 0.5, 0.975), meanlog = 1, sdlog = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.lnorm.par(q = q)
#' get.lnorm.par(q = q, fit.weights = c(100, 1, 100))
#' get.lnorm.par(q = q, fit.weights = c(10, 1, 10))
#' get.lnorm.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.000001, 0.99999999))
#' get.lnorm.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qlnorm(p = c(0.025, 0.5, 0.975), meanlog = 0.1, sdlog = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.lnorm.par(q = q)
#' get.lnorm.par(q = q, fit.weights = c(100, 1, 100))
#' get.lnorm.par(q = q, fit.weights = c(10, 1, 10))
#' get.lnorm.par(q = q, fit.weights = c(1, 100, 1))
#' get.lnorm.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qlnorm(p = c(0.025, 0.975), meanlog = 4, sdlog = 0.8)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.lnorm.par(p = c(0.025, 0.975), q = q)
#' get.lnorm.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1), scaleX = c(0.1, 0.001))
#' get.lnorm.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 100), scaleX = c(0.1, 0.001))
#' get.lnorm.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' get.lnorm.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 10))
#' graphics::par(old.par)
#'
get.lnorm.par <- function(p = c(0.025, 0.5, 0.975), q, 
                          show.output = TRUE, plot = TRUE, 
                          tol = 0.001, fit.weights = rep(1, length(p)), 
                          scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) <= 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain (0, inf) => Lognormal distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::plnorm(q = q, meanlog = theta[1], sdlog = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(1, 0.35), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(-10000, 0.001), 
                            upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(1, 3),
                                minimize, 
                                method = "Nelder-Mead"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'Nelder-Mead' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'Nelder-Mead' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("meanlog", "sdlog")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("meanlog", "sdlog")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("meanlog = ", round(Par["meanlog"], digits = 2))
        main2 <- paste("sdlog = ", round(Par["sdlog"], digits = 2))
        main <- paste("Lognormal (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qlnorm(p = min(p) * scaleX[1], 
                                       meanlog = Par["meanlog"], 
                                       sdlog = Par["sdlog"]),
                         stats::qlnorm(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                       meanlog = Par["meanlog"], 
                                       sdlog = Par["sdlog"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::plnorm(Support, Par["meanlog"], Par["sdlog"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.logis.par} returns the parameters of a logistic distribution
#' where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least two. Using the default \code{p},
#' the three corresponding quantiles are the 2.5th percentile, the median and the
#' 97.5th percentile, respectively. \code{get.logis.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{BFGS} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.logis.par
#' @aliases get.logis.par
#' @title Fitting parameters of a logistic distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.logis.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'    plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a logistic distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{plogis} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qlogis(p = c(0.025, 0.5, 0.975), location = 0, scale = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.logis.par(q = q)
#' get.logis.par(q = q, scaleX = c(0.5, 0.5))
#' get.logis.par(q = q, fit.weights = c(100, 1, 100))
#' get.logis.par(q = q, fit.weights = c(10, 1, 10))
#' get.logis.par(q = q, fit.weights = c(1, 100, 1))
#' get.logis.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qlogis(p = c(0.025, 0.5, 0.975), location = 0, scale = 3)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.logis.par(q = q)
#' get.logis.par(q = q, fit.weights = c(100, 1, 100))
#' get.logis.par(q = q, fit.weights = c(10, 1, 10))
#' get.logis.par(q = q, fit.weights = c(1, 100, 1))
#' get.logis.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qlogis(p = c(0.025, 0.975), location = 0, scale = 3)
#' old.par <- graphics::par(mfrow = c(1, 3))
#' get.logis.par(p = c(0.025, 0.975), q = q)
#' get.logis.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1))
#' get.logis.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' graphics::par(old.par)
#'
get.logis.par <- function(p = c(0.025, 0.5, 0.975), q, 
                          show.output = TRUE, plot = TRUE, 
                          tol = 0.001, fit.weights = rep(1, length(p)), 
                          scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::plogis(q = q, 
                                                  location = theta[1], 
                                                  scale = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    lm <- stats::lm(q ~ p)
    suppressWarnings(m <- stats::predict(lm, newdata = list(p = 0.5))[[1]])
    suppressWarnings(s <- stats::sd(q)/pi * sqrt(2))
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(m, s), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(-10000, 0.001), 
                            upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(m, s),
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("location", "scale")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("location", "scale")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("location = ", round(Par["location"], digits = 2))
        main2 <- paste("scale = ", round(Par["scale"], digits = 2))
        main <- paste("Logistic (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qlogis(p = min(p) * scaleX[1], 
                                       location = Par["location"], 
                                       scale = Par["scale"]),
                         stats::qlogis(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                       location = Par["location"], 
                                       scale = Par["scale"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::plogis(Support, Par["location"], Par["scale"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.nbinom.par} returns the parameters of a negative binomial distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least two. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.nbinom.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.nbinom.par
#' @aliases get.nbinom.par
#' @title Fitting parameters of a negative binomial distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.nbinom.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a negative binomial distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pnbinom} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qnbinom(p = c(0.025, 0.5, 0.975), size = 10, prob = 0.5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.nbinom.par(q = q)
#' get.nbinom.par(q = q, fit.weights = c(100, 1, 100))
#' get.nbinom.par(q = q, fit.weights = c(1, 100, 1))
#' get.nbinom.par(q = q, fit.weights = c(10, 1, 10))
#' get.nbinom.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qnbinom(p = c(0.025, 0.5, 0.975), size = 1, prob = 0.5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.nbinom.par(q = q, tol = 0.01)
#' get.nbinom.par(q = q, fit.weights = c(100, 1, 100))
#' get.nbinom.par(q = q, fit.weights = c(1, 100, 1), tol = 0.01)
#' get.nbinom.par(q = q, fit.weights = c(10, 1, 10), tol = 0.01)
#' get.nbinom.par(q = q, fit.weights = c(1, 10, 1), tol = 0.01)
#' graphics::par(old.par)
#'
#' q <- stats::qnbinom(p = c(0.025, 0.5, 0.975), size = 1, prob = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.nbinom.par(q = q)
#' get.nbinom.par(q = q, fit.weights = c(100, 1, 100))
#' get.nbinom.par(q = q, fit.weights = c(1, 100, 1))
#' get.nbinom.par(q = q, fit.weights = c(10, 1, 10))
#' get.nbinom.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qnbinom(p = c(0.025, 0.975), size = 10, prob = 0.5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.nbinom.par(p = c(0.025, 0.975), q = q,)
#' get.nbinom.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1))
#' get.nbinom.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 100))
#' get.nbinom.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' get.nbinom.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 10))
#' graphics::par(old.par)
#'
get.nbinom.par <- function(p = c(0.025, 0.5, 0.975), q,
                           show.output = TRUE, plot = TRUE, 
                           tol = 0.001, fit.weights = rep(1, length(p)), 
                           scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) < 0 | any(abs(q - trunc(q)) != 0)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles must be integer and positive => negative binomial distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::pnbinom(q = q, 
                                                   size = theta[1], 
                                                   prob = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    sizeStart <- mean(q)^2/(stats::sd(q)^2 - mean(q))
    probStart<-(stats::sd(q)^2 - mean(q))/stats::sd(q)^2
    try1 <- try(
        fit <- stats::optim(par = c(sizeStart, probStart), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(0.001, 0.001), 
                            upper = c(10000, 0.999)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(sizeStart, probStart),
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("size", "prob")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("size", "prob")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("size = ", round(Par["size"], digits = 2), sep = "")
        main2 <- paste("prob = ", round(Par["prob"], digits = 2), sep = "")
        main <- paste("negbin (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qnbinom(p = min(p) * scaleX[1], 
                                        size = Par["size"], 
                                        prob = Par["prob"]),
                         stats::qnbinom(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                        size = Par["size"], 
                                        prob = Par["prob"]))
        Support <- seq(min(min(q), Support.lim[1]), max(max(q), Support.lim[2]), length = 200)
        Probability <- stats::pnbinom(Support, Par["size"], Par["prob"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.norm.par} returns the parameters of a normal distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings 
#' must be identical and should be at least two. Using the default \code{p}, the 
#' three corresponding quantiles are the 2.5th percentile, the median and the 
#' 97.5th percentile, respectively. \code{get.norm.par} uses the R function \code{optim} with the
#' method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{BFGS} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item 
#' \code{value} displays the achieved minimal value of the functions that were minimized. 
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component 
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1. 
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation). 
#'
#' @name get.norm.par
#' @aliases get.norm.par
#' @title Fitting parameters of normal distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting),  \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting) 
#' @usage get.norm.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE, 
#'    plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities. 
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a normal distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated 
#' and the theoretical distribution parameters in certain circumstances. This is 
#' because the estimation of the parameters is based on a numerical optimization 
#' method and depends strongly on the initial values. In addition, one must always 
#' keep in mind that a distribution for different combinations of parameters may 
#' look very similar. Therefore, the optimization method cannot always find the 
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or 
#' specified tolerance not achieved", one may try to set the convergence tolerance 
#' to a higher value. It is yet to be noted, that good till very good fits of parameters 
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pnorm} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qnorm(p = c(0.025, 0.5, 0.975), mean = 12, sd = 34)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.norm.par(q = q)
#' get.norm.par(q = q, scaleX = c(0.00001, 0.99999))
#' get.norm.par(q = q, fit.weights = c(10, 1, 10))
#' get.norm.par(q = q, fit.weights = c(1, 10, 1))
#' get.norm.par(q = q, fit.weights = c(100, 1, 100))
#' get.norm.par(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qnorm(p = c(0.025, 0.5, 0.975), mean = 0, sd = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.norm.par(q = q)
#' get.norm.par(q = q, fit.weights = c(10, 1, 10))
#' get.norm.par(q = q, fit.weights = c(1, 10, 1))
#' get.norm.par(q = q, fit.weights = c(100, 1, 100))
#' get.norm.par(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qnorm(p = c(0.025, 0.5, 0.975), mean = 0.1, sd = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.norm.par(q = q)
#' get.norm.par(q = q, fit.weights = c(10, 1, 10))
#' get.norm.par(q = q, fit.weights = c(1, 10, 1))
#' get.norm.par(q = q, fit.weights = c(100, 1, 100))
#' get.norm.par(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qnorm(p = c(0.025, 0.975), mean = 12, sd = 34)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.norm.par(p = c(0.025, 0.975), q = q)
#' get.norm.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' get.norm.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1))
#' get.norm.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 10))
#' get.norm.par(p = c(0.025, 0.975), q = q, fit.weights = c(1, 100))
#' graphics::par(old.par)
#'
get.norm.par <- function(p = c(0.025, 0.5, 0.975), q, 
                         show.output = TRUE, plot = TRUE, 
                         tol = 0.001, fit.weights = rep(1, length(p)),
                         scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    } 
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    lm.fit <- stats::lm(q ~ p)
    suppressWarnings(m <- stats::predict(lm.fit, newdata = list(p = 0.5))[[1]])
    suppressWarnings(s <- (stats::predict(lm.fit, newdata = list(p = 0.975))[[1]] - m)/1.96)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::pnorm(q = q, 
                                                 mean = theta[1], 
                                                 sd = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(m, s), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(-10000, 0.001), 
                            upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(m, s),
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("mean", "sd")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("mean", "sd")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plot graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("mean = ", round(Par["mean"], digits = 2))
        main2 <- paste("sd = ", round(Par["sd"], digits = 2))
        main <- paste("Normal (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qnorm(p = min(p) * scaleX[1], 
                                      mean = Par["mean"], 
                                      sd = Par["sd"]),
                         stats::qnorm(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                      mean = Par["mean"], 
                                      sd = Par["sd"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pnorm(Support, Par["mean"], Par["sd"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), 
                       main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)   
}


################################################################################
#' \code{get.norm.sd} returns the standard deviation of a normal distribution
#' where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities and the number of quantiles must be identical and
#' should be at least two. \code{get.norm.sd} uses the central limit theorem and
#' the linear regression.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{lm} will be shown.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function will be meaningful only if the quantile comes from a normal distribution.
#'
#' @name get.norm.sd
#' @aliases get.norm.sd
#' @title Fitting standard deviation of a normal distribution from one or more quantiles and known mean
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.norm.sd(p = c(0.025, 0.5, 0.975), q, show.output = TRUE, plot = TRUE,
#'                    fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, vector of probabilities.
#' @param q numeric, vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default vaule is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns an estimated standard deviation or missing value
#' @note It should be noted that the data must be normally distributed, or the 
#' central limt theorem must hold for large (enough) samples sizes.  
#' @seealso See \code{pnorm} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qnorm(p = c(0.025, 0.5, 0.975), mean = 0, sd = 2)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.norm.sd(q = q)
#' get.norm.sd(q = q, scaleX = c(0.0001, 0.9999))
#' get.norm.sd(q = q, fit.weights = c(10, 1, 10))
#' get.norm.sd(q = q, fit.weights = c(1, 10, 1))
#' get.norm.sd(q = q, fit.weights = c(100, 1, 100))
#' get.norm.sd(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qnorm(p = c(0.025, 0.5, 0.975), mean = 176, sd = 15)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.norm.sd(q = q)
#' get.norm.sd(q = q, fit.weights = c(10, 1, 10))
#' get.norm.sd(q = q, fit.weights = c(1, 10, 1))
#' get.norm.sd(q = q, fit.weights = c(100, 1, 100))
#' get.norm.sd(q = q, fit.weights = c(1, 100, 1))
#' graphics::par(old.par)
#'
#' ## The estimation model is not suitable for the following quantiles.
#' ## Because the quantile is unsymmetrical, which could not be from a normally distributed data.
#' q <- c(-2, 30, 31)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.norm.sd(q = q)
#' get.norm.sd(q = q, fit.weights = c(10, 1, 10))
#' get.norm.sd(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.0001, 0.9999))
#' get.norm.sd(q = q, fit.weights = c(100, 1, 100))
#' get.norm.sd(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.0001, 0.9999))
#' graphics::par(old.par)
#'
#' ## Estimating from actually exponentially distributed data
#' x.exp <- rexp(n = 10, rate = 5)
#' mean(x.exp)
#' stats::sd(x.exp)
#' q <- quantile(x.exp, c(0.025, 0.5, 0.975))
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.norm.sd(q = q)
#' get.norm.sd(q = q, fit.weights = c(1, 10, 1))
#' get.norm.sd(q = q, fit.weights = c(10, 1, 10))
#' get.norm.sd(q = q, fit.weights = c(1, 100, 1))
#' get.norm.sd(q = q, fit.weights = c(100, 1, 100))
#' graphics::par(old.par)
#'
#' ## other examples
#' q <- stats::qnorm(p = c(0.025, 0.5, 0.975), mean = 1, sd = 1)
#' get.norm.sd(q = q)
#'
#' q <- stats::qnorm(p = c(0.025, 0.5, 0.975), mean = 1, sd = 0.5)
#' get.norm.sd(q = q)
#'
#' q <- stats::qnorm(p = c(0.025, 0.5, 0.975), mean = 0.01, sd = 0.1)
#' get.norm.sd(q = q)
#' 
get.norm.sd <- function(p = c(0.025, 0.5, 0.975), q,
                        show.output = TRUE, plot = TRUE, 
                        fit.weights = rep(1, length(p)), 
                        scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(fit.weights) < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the weighting vector should be positive", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least one quantile must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # estimating procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    q.theor <- stats::qnorm(p)
    lmodel <- stats::lm(q ~ q.theor, weights = fit.weights)
    
    #-----------------------------------------------------------------------------
    # checking output
    #-----------------------------------------------------------------------------
    if (show.output) {
        lmodel$par <- c(lmodel$coeff[1][[1]], lmodel$coeff[2][[1]])
        names(lmodel$par) <- c("mean", "sd")
    }
    
    #-----------------------------------------------------------------------------
    # creating graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(lmodel$par)) & plot) {
        main1 <- paste("mean = ", round(lmodel$par["mean"], digits = 2), 
                       ", sd = ", round(lmodel$par["sd"], digits = 2))
        main <- paste("Normal (", main1, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qnorm(p = min(p) * scaleX[1], 
                                      sd = lmodel$par["sd"], 
                                      mean = lmodel$par["mean"]),
                         stats::qnorm(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                      sd = lmodel$par["sd"], 
                                      mean = lmodel$par["mean"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pnorm(Support, mean = lmodel$par["mean"], sd = lmodel$par["sd"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(lmodel$par)
}



################################################################################
#' \code{get.pert.par} returns the parameters of a pert distribution
#' where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least three. Using the default \code{p},
#' the three corresponding quantiles are the 2.5th percentile, the median and the
#' 97.5th percentile, respectively. \code{get.pert.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{BFGS} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.pert.par
#' @aliases get.pert.par
#' @title Fitting parameters of a pert distribution from four or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.pert.par(p = c(0.025, 0.5, 0.6, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a pert distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{\link[mc2d:ppert]{ppert}} from the package \pkg{mc2d} for distribution 
#' implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- mc2d::qpert(p = c(0.025, 0.5, 0.6, 0.975), min = 0, mode = 3, max = 10, shape = 5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.pert.par(q = q)
#' get.pert.par(q = q, fit.weights = c(100, 1, 1, 100))
#' get.pert.par(q = q, fit.weights = c(10, 1, 1, 10))
#' get.pert.par(q = q, fit.weights = c(1, 100, 1, 1))
#' get.pert.par(q = q, fit.weights = c(1, 10, 1, 1))
#' graphics::par(old.par)
#'
#' q <- mc2d::qpert(p = c(0.025, 0.5, 0.6, 0.975), min = 1, mode = 5, max = 10, shape = 4)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.pert.par(q = q)
#' get.pert.par(q = q, scaleX = c(0.0001, 0.999999))
#' get.pert.par(q = q, fit.weights = c(100, 1, 1, 100))
#' get.pert.par(q = q, fit.weights = c(10, 1, 1, 10))
#' get.pert.par(q = q, fit.weights = c(1, 100, 1, 1))
#' get.pert.par(q = q, fit.weights = c(1, 10, 1, 1))
#' graphics::par(old.par)
#'
#' q <- mc2d::qpert(p = c(0.025, 0.5, 0.6, 0.975), min=-10, mode = 5, max = 10, shape = 4)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.pert.par(q = q)
#' get.pert.par(q = q, fit.weights = c(100, 1, 1, 100))
#' get.pert.par(q = q, fit.weights = c(10, 1, 1, 10))
#' get.pert.par(q = q, fit.weights = c(1, 100, 1, 1))
#' get.pert.par(q = q, fit.weights = c(1, 10, 1, 1))
#' graphics::par(old.par)
#'
#' q <- mc2d::qpert(p = c(0.025, 0.5, 0.6, 0.975), min=-10, mode = 5, max = 10, shape = 0.4)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.pert.par(q = q)
#' get.pert.par(q = q, fit.weights = c(100, 1, 1, 100))
#' get.pert.par(q = q, fit.weights = c(10, 1, 1, 10))
#' get.pert.par(q = q, fit.weights = c(1, 100, 1, 1))
#' get.pert.par(q = q, fit.weights = c(1, 10, 1, 1))
#' graphics::par(old.par)
#'
get.pert.par <- function(p = c(0.025, 0.5, 0.6, 0.975), q,
                         show.output = TRUE, plot = TRUE, 
                         tol = 0.001, fit.weights = rep(1, length(p)), 
                         scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 4) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least four quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(mc2d::ppert(q = q, 
                                                min = theta[1], 
                                                mode = theta[2], 
                                                max = theta[3], 
                                                shape = theta[4]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(1, 5, 10, 4), 
                            minimize, method = "L-BFGS-B",
                            lower = c(-10000, -10000, -10000, 0.001), 
                            upper = c(10000, 10000, 10000, 1000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(1, 5, 10, 4), 
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("min", "mode", "max", "shape")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("min", "mode", "max", "shape")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("min = ", round(Par["min"], digits = 2), sep = "")
        main2 <- paste("mode = ", round(Par["mode"], digits = 2), sep = "")
        main3 <- paste("max = ", round(Par["max"], digits = 2), sep = "")
        main4 <- paste("shape = ", round(Par["shape"], digits = 2), sep = "")
        main <- paste("Pert (", main1, ", ", main2, ", ", main3, ", ", main4, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(mc2d::qpert(p = min(p) * scaleX[1], 
                                     min = Par["min"], 
                                     mode = Par["mode"], 
                                     max = Par["max"], 
                                     shape = Par["shape"]),
                         mc2d::qpert(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                     min = Par["min"], 
                                     mode = Par["mode"], 
                                     max = Par["max"], 
                                     shape = Par["shape"]))
        Support <- seq(min(min(q), Support.lim[1]), max(max(q), Support.lim[2]), length = 200)
        Probability <- mc2d::ppert(Support, Par["min"], Par["mode"], 
                                   Par["max"], shape = Par["shape"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.pois.par} returns the parameters of a Poisson distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least one. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.pois.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.pois.par
#' @aliases get.pois.par
#' @title Fitting parameter of Poisson distribution from one or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.pois.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities
#' @param q numeric, single value or vector of quantiles corresponding to p
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE})
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE})
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001})
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE})
#' @return Returns fitted parameters of a Poisson distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because it is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{ppois} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qpois(p = c(0.025, 0.5, 0.975), lambda = 3)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.pois.par(q = q)
#' get.pois.par(q = q, fit.weights = c(100, 1, 100))
#' get.pois.par(q = q, fit.weights = c(10, 1, 10))
#' get.pois.par(q = q, fit.weights = c(1, 100, 1))
#' get.pois.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qpois(p = c(0.025, 0.5, 0.975), lambda = 4)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.pois.par(q = q)
#' get.pois.par(q = q, fit.weights = c(100, 1, 100))
#' get.pois.par(q = q, fit.weights = c(10, 1, 10))
#' get.pois.par(q = q, fit.weights = c(1, 100, 1))
#' get.pois.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qpois(p = c(0.025, 0.5, 0.975), lambda = 0.5)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.pois.par(q = q, tol = 1)
#' get.pois.par(q = q, fit.weights = c(100, 1, 100), tol = 1)
#' get.pois.par(q = q, fit.weights = c(10, 1, 10), tol = 1)
#' get.pois.par(q = q, fit.weights = c(1, 100, 1))
#' get.pois.par(q = q, fit.weights = c(1, 10, 1), tol = 0.01)
#' graphics::par(old.par)
#'
#' q <- stats::qpois(p = c(0.025, 0.5, 0.975), lambda = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.pois.par(q = q, tol = 0.01)
#' get.pois.par(q = q, fit.weights = c(100, 1, 100), tol = 0.01)
#' get.pois.par(q = q, fit.weights = c(10, 1, 10), tol = 0.01)
#' get.pois.par(q = q, fit.weights = c(1, 100, 1))
#' get.pois.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
get.pois.par <- function(p = c(0.025, 0.5, 0.975), q, 
                         show.output = TRUE, plot = TRUE, 
                         tol = 0.001, fit.weights = rep(1, length(p)), 
                         scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) < 0 | any(abs(q-trunc(q)) != 0)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles must be integer and positive => Poisson distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least one quantile must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(lambda) {
        summand <- stats::ppois(q = q, lambda = lambda) - p
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = 1, 
                            minimize, method = "L-BFGS-B", 
                            lower = 0.001, upper = 10000), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        Par <- NA
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- as.integer(fit$par)
        names(Par) <- c("lambda")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("lambda = ", round(Par["lambda"], digits = 2), sep = "")
        main <- paste("Poisson (", main1, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qpois(p = min(p) * scaleX[1], 
                                      lambda = Par["lambda"]),
                         stats::qpois(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                      lambda = Par["lambda"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::ppois(Support, Par["lambda"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}



################################################################################
#' \code{get.t.par} returns the parameters of a Student's t distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities and the number of quantiles must be identical
#' and should be at least one. Using the default \code{p}, the three corresponding
#' quantiles are the 2.5th percentile, the median and the 97.5th percentile,
#' respectively. \code{get.t.par} uses the R function \code{optim} with the
#' method \code{L-BFGS-B}. If this method fails the optimization method \code{BFGS}
#' will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.t.par
#' @aliases get.t.par
#' @title Fitting parameter of a Student's t distribution from one or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.t.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a Student's t distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{stats::pt} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qt(p = c(0.025, 0.5, 0.975), df = 10)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.t.par(q = q)
#' get.t.par(q = q, fit.weights = c(100, 1, 100))
#' get.t.par(q = q, fit.weights = c(10, 1, 10))
#' get.t.par(q = q, fit.weights = c(1, 100, 1))
#' get.t.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qt(p = c(0.025, 0.5, 0.975), df = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.t.par(q = q, scaleX = c(0.5, 0.5))
#' get.t.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.5, 0.5))
#' get.t.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.5, 0.5))
#' get.t.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.5, 0.5))
#' get.t.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.5, 0.5))
#' graphics::par(old.par)
#'
#' q <- stats::qt(p = c(0.025, 0.5, 0.975), df = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.t.par(q = q, scaleX = c(0.5, 0.5))
#' get.t.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.5, 0.5))
#' get.t.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.5, 0.5))
#' get.t.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.5, 0.5))
#' get.t.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.5, 0.5))
#' graphics::par(old.par)
#'
#' ## example with only one quantile
#' q <- stats::qt(p = c(0.025), df = 3)
#' old.par <- graphics::par(mfrow = c(1, 3))
#' get.t.par(p = c(0.025), q = q)
#' get.t.par(p = c(0.025), q = q, fit.weights = 10)
#' get.t.par(p = c(0.025), q = q, fit.weights = 100)
#' graphics::par(old.par)
#'
get.t.par <- function(p = c(0.025, 0.5, 0.975), q, 
                      show.output = TRUE, plot = TRUE, 
                      tol = 0.001, fit.weights = rep(1, length(p)), 
                      scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least one quantile must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(par) {
        summand <- suppressWarnings(stats::pt(q = q, df = par) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = 1, 
                            minimize, method = "L-BFGS-B", 
                            lower = 0.001, upper = 10000), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = 1, 
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("df")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("df")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("df = ", round(Par["df"], digits = 2))
        main <- paste("Student (", main1, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qt(p = min(p) * scaleX[1], 
                                   df = Par["df"]), 
                         stats::qt(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                   df = Par["df"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pt(Support, Par["df"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}



################################################################################
#' \code{get.tnorm.par} returns the parameters of a truncated normal distribution
#' where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least four. Using the default \code{p}, the
#' four corresponding quantiles are the 2.5th percentile, the median, the 75th
#' percentile and the 97.5th percentile, respectively. \code{get.tnorm.par} uses
#' the R function \code{optim} with the method \code{L-BFGS-B}. If this method
#' fails the optimization method \code{Nelder-Mead} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.tnorm.par
#' @aliases get.tnorm.par
#' @title Fitting parameters of truncated normal distribution from four or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting),  \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.tnorm.par(p = c(0.025, 0.5, 0.75, 0.975), q, show.output = TRUE,
#'    plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a truncated normal distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{\link[msm:ptnorm]{ptnomr}} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- msm::qtnorm(p = c(0.025, 0.5, 0.75, 0.975), mean = 3, sd = 3, lower = 0, upper = 10)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.tnorm.par(q = q)
#' get.tnorm.par(q = q, scaleX = c(0.1, 0.999999))
#' get.tnorm.par(q = q, fit.weights = c(100, 1, 1, 100))
#' get.tnorm.par(q = q, fit.weights = c(10, 1, 1, 10))
#' get.tnorm.par(q = q, fit.weights = c(1, 100, 1, 1))
#' get.tnorm.par(q = q, fit.weights = c(1, 10, 1, 1))
#' graphics::par(old.par)
#'
#' q <- msm::qtnorm(p = c(0.025, 0.5, 0.75, 0.975), mean = 3, sd = 0.1, lower=-1, upper = 4)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.tnorm.par(q = q)
#' get.tnorm.par(q = q, fit.weights = c(100, 1, 1, 100))
#' get.tnorm.par(q = q, fit.weights = c(10, 1, 1, 10))
#' get.tnorm.par(q = q, fit.weights = c(1, 100, 1, 1))
#' get.tnorm.par(q = q, fit.weights = c(1, 10, 1, 1))
#' graphics::par(old.par)
#'
#' q <- msm::qtnorm(p = c(0.025, 0.5, 0.75, 0.975), mean = 0, sd = 1, lower=-2, upper = 2)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.tnorm.par(q = q)
#' get.tnorm.par(q = q, fit.weights = c(100, 1, 1, 100))
#' get.tnorm.par(q = q, fit.weights = c(10, 1, 1, 10))
#' get.tnorm.par(q = q, fit.weights = c(1, 100, 1, 1))
#' get.tnorm.par(q = q, fit.weights = c(1, 10, 1, 1))
#' graphics::par(old.par)
#'
get.tnorm.par <- function(p = c(0.025, 0.5, 0.75, 0.975), q, 
                          show.output = TRUE, plot = TRUE, 
                          tol = 0.001, fit.weights = rep(1, length(p)), 
                          scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 4) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least four quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    lm <- stats::lm(q ~ p)
    suppressWarnings(m <- stats::predict(lm, newdata = list(p = 0.5))[[1]])
    suppressWarnings(s <- (stats::predict(lm, newdata = list(p = 0.975))[[1]] - m)/1.96)
    minimize <- function(theta) {
        summand <- suppressWarnings(msm::ptnorm(q = q, 
                                                mean = theta[1], 
                                                sd = theta[2], 
                                                lower = theta[3], 
                                                upper = theta[4]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(m, s, min(q) - s, max(q) + s), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(-10000, 0.001, -10000, -10000), 
                            upper = c(10000, 10000, 10000, 10000)), 
        silent = TRUE)
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(m, s, min(q) - s, max(q) + s), 
                                minimize, method = "Nelder-Mead"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'Nelder-Mead' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'Nelder-Mead' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("mean", "sd", "lower", "upper")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("mean", "sd", "lower", "upper")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plot graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("mean = ", round(Par["mean"], digits = 2))
        main2 <- paste("sd = ", round(Par["sd"], digits = 2))
        main3 <- paste("lower = ", round(Par["lower"], digits = 2))
        main4 <- paste("upper = ", round(Par["upper"], digits = 2))
        main <- paste("trunc. normal (", main1, ", ", main2, ", ", main3, ", ", main4, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(msm::qtnorm(p = min(p) * scaleX[1], 
                                     mean = Par["mean"], 
                                     sd = Par["sd"], 
                                     lower = Par["lower"], 
                                     upper = Par["upper"]),
                         msm::qtnorm(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                     mean = Par["mean"], 
                                     sd = Par["sd"], 
                                     lower = Par["lower"], 
                                     upper = Par["upper"]))
        Support <- seq(min(min(q), Support.lim[1]), max(max(q), Support.lim[2]), length = 200)
        Probability <- msm::ptnorm(Support, Par["mean"], Par["sd"], Par["lower"], Par["upper"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.triang.par} returns the parameters of a triangular distribution
#' where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least three. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.triang.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}. If this method fails the optimization method
#' \code{BFGS} will be invoked.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.triang.par
#' @aliases get.triang.par
#' @title Fitting parameters of a triangular distribution from three or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.triang.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'     plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a triangular distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{\link[mc2d:ptriang]{ptriang}} from the package \pkg{mc2d} for distribution 
#' implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- mc2d::qtriang(p = c(0.025, 0.5, 0.975), min = 0, mode = 3, max = 10)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.triang.par(q = q)
#' get.triang.par(q = q, fit.weights = c(100, 1, 100))
#' get.triang.par(q = q, fit.weights = c(10, 1, 10))
#' get.triang.par(q = q, fit.weights = c(1, 100, 1))
#' get.triang.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- mc2d::qtriang(p = c(0.025, 0.5, 0.975), min = 1, mode = 5, max = 10)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.triang.par(q = q)
#' get.triang.par(q = q, scaleX = c(0.00001, 0.99999))
#' get.triang.par(q = q, fit.weights = c(100, 1, 100))
#' get.triang.par(q = q, fit.weights = c(10, 1, 10))
#' get.triang.par(q = q, fit.weights = c(1, 100, 1))
#' get.triang.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## bad fit for negative values
#' q <- mc2d::qtriang(p = c(0.025, 0.5, 0.975), min=-20, mode = 5, max = 10)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.triang.par(q = q, tol = 0.1)
#' get.triang.par(q = q)
#' get.triang.par(q = q, fit.weights = c(100, 1, 100))
#' get.triang.par(q = q, fit.weights = c(10, 1, 10))
#' get.triang.par(q = q, fit.weights = c(1, 100, 1), tol = 1)
#' get.triang.par(q = q, fit.weights = c(1, 10, 1), tol = 1)
#' graphics::par(old.par)
#'
#' ## other examples
#' q <- mc2d::qtriang(p = c(0.025, 0.5, 0.975), min=-20, mode = 5, max = 10)
#' get.triang.par(q = q, tol = 0.3)
#'
get.triang.par <- function(p = c(0.025, 0.5, 0.975), q, 
                           show.output = TRUE, plot = TRUE, 
                           tol = 0.001, fit.weights = rep(1, length(p)), 
                           scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 3) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least three quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(mc2d::ptriang(q = q, 
                                                  min = theta[1], 
                                                  mode = theta[2], 
                                                  max = theta[3]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(1, 5, 10), 
                            minimize, method = "L-BFGS-B",
                            lower = c(-10000,-10000,-10000), 
                            upper = c(10000, 10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(1, 5, 10), 
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("min", "mode", "max")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("min", "mode", "max")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("min = ", round(Par["min"], digits = 2), sep = "")
        main2 <- paste("mode = ", round(Par["mode"], digits = 2), sep = "")
        main3 <- paste("max = ", round(Par["max"], digits = 2), sep = "")
        main <- paste("Triangular (", main1, ", ", main2, ", ", main3, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(mc2d::qtriang(p = min(p) * scaleX[1], 
                                       min = Par["min"], 
                                       mode = Par["mode"], 
                                       max = Par["max"]),
                         mc2d::qtriang(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                       min = Par["min"], 
                                       mode = Par["mode"], 
                                       max = Par["max"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- mc2d::ptriang(Support, Par["min"], Par["mode"], Par["max"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
}


################################################################################
#' \code{get.unif.par} returns the parameters of a uniform distribution
#' where the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities and the number of quantiles must be identical
#' and should be at least two. Using the default \code{p}, the three corresponding
#' quantiles are the 2.5th percentile, the median and the 97.5th percentile,
#' respectively.
#' \cr \cr
#' Parameters of the uniform distribution are estimated exactly.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#'
#' @name get.unif.par
#' @aliases get.unif.par
#' @title Fitting parameters of a uniform distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.unif.par(p = c(0.025, 0.975), q, plot = TRUE, scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE})
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE})
#' @return Returns fitted parameters of a uniform distribution.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qunif(p = c(0.025, 0.975), min = 0, max = 5)
#' get.unif.par(q = q)
#' get.unif.par(q = q, scaleX = c(0.001, 0.999))
#'
#' q <- stats::qunif(p = c(0.025, 0.975), min=-6, max = 5)
#' get.unif.par(q = q)
#' 
get.unif.par <- function(p = c(0.025, 0.975), q, 
                         plot = TRUE, 
                         scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p' and/or 'q'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (length(q) != 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, uniform distribution can be fitted only by TWO given percentiles!", call. = FALSE)
    }
    if (length(p) != length(q)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p' and 'q' are not of the same length! The vector of quantile should be of the same length as a corresponding vector of probabilities!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # defining parameters
    #-----------------------------------------------------------------------------
    weight1 <- p[1]/(p[2] - p[1])
    weight2 <- (1 - p[2])/(p[2] - p[1])
    a <- q[1] - weight1 * (q[2] - q[1])
    b <- q[2] + weight2 * (q[2] - q[1])
    Par <- c(min = a, max = b)
    
    #-----------------------------------------------------------------------------
    # plot graphical diagnostics
    #-----------------------------------------------------------------------------
    if (plot) {
        main1 <- paste("min = ", round(Par["min"], digits = 2), sep = "")
        main2 <- paste("max = ", round(Par["max"], digits = 2), sep = "")
        main <- paste("Uniform (", main1, ", ", main2, ")", sep = "")
        Support.lim <- c(stats::qunif(p = min(p) * scaleX[1], 
                                      min = Par["min"], 
                                      max = Par["max"]),
                         stats::qunif(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                      min = Par["min"], 
                                      max = Par["max"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::punif(Support, Par["min"], Par["max"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    #-----------------------------------------------------------------------------
    # plot graphical diagnostics
    #-----------------------------------------------------------------------------
    return(Par)
}



################################################################################
#' \code{get.weibull.par} returns the parameters of a Weibull distribution where
#' the \code{p}th percentiles match with the quantiles \code{q}.
#'
#' The number of probabilities, the number of quantiles and the number of weightings
#' must be identical and should be at least two. Using the default \code{p}, the
#' three corresponding quantiles are the 2.5th percentile, the median and the 97.5th
#' percentile, respectively. \code{get.weibull.par} uses the R function \code{optim}
#' with the method \code{L-BFGS-B}.
#' \cr \cr
#' If \code{show.output = TRUE} the output of the function \code{optim} will be shown.
#' The item \code{convergence} equal to 0 means the successful completion of the
#' optimization procedure, otherwise it indicates a convergence error. The item
#' \code{value} displays the achieved minimal value of the functions that were minimized.
#' \cr \cr
#' The estimated distribution parameters returned by the function \code{optim} are
#' accepted if the achieved value of the minimized function (output component
#' \code{value} of \code{optim}) is smaller than the argument \code{tol}.
#' \cr \cr
#' The items of the probability vector \code{p} should lie between 0 and 1.
#' \cr \cr
#' The items of the weighting vector \code{fit.weights} should be positive values.
#' \cr \cr
#' The function which will be minimized is defined as a sum of squared differences
#' between the given probabilities and the theoretical probabilities of the specified
#' distribution evaluated at the given quantile points (least squares estimation).
#'
#' @name get.weibull.par
#' @aliases get.weibull.par
#' @title Fitting parameters of a Weibull distribution from two or more quantiles
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR),  \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting),  \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage get.weibull.par(p = c(0.025, 0.5, 0.975), q, show.output = TRUE,
#'    plot = TRUE, tol = 0.001, fit.weights = rep(1, length(p)), scaleX = c(0.1, 0.9), ...)
#' @param p numeric, single value or vector of probabilities.
#' @param q numeric, single value or vector of quantiles corresponding to p.
#' @param show.output logical, if \code{TRUE} the \code{optim} result will be printed (default value is \code{TRUE}).
#' @param plot logical, if \code{TRUE} the graphical diagnostics will be plotted (default value is \code{TRUE}).
#' @param tol numeric, single positive value giving the absolute convergence tolerance for reaching zero (default value is \code{0.001}).
#' @param fit.weights numerical vector of the same length as a probabilities vector 
#'    \code{p} containing positive values for weighting quantiles. By default all
#'    quantiles will be weighted by 1.
#' @param scaleX numerical vector of the length 2 containing values (from the open interval (0, 1))
#' for scaling quantile-axis (relevant only if \code{plot = TRUE}). The smaller the left value,
#' the further the graph is extrapolated within the lower percentile, the greater the right
#' value, the further it goes within the upper percentile.
#' @param ...	further arguments passed to the functions \code{plot} and \code{points} (relevant only if \code{plot = TRUE}).
#' @return Returns fitted parameters of a Weibull distribution or missing
#' values (\code{NA}'s) if the distribution cannot fit the specified quantiles.
#' @note It should be noted that there might be deviations between the estimated
#' and the theoretical distribution parameters in certain circumstances. This is
#' because the estimation of the parameters is based on a numerical optimization
#' method and depends strongly on the initial values. In addition, one must always
#' keep in mind that a distribution for different combinations of parameters may
#' look very similar. Therefore, the optimization method cannot always find the
#' "right" distribution, but a "similar" one.
#' \cr \cr
#' If the function terminates with the error massage "convergence error occurred or
#' specified tolerance not achieved", one may try to set the convergence tolerance
#' to a higher value. It is yet to be noted, that good till very good fits of parameters
#' could only be obtained for tolerance values that are smaller than 0.001.
#' @seealso See \code{pweibull} for distribution implementation details.
#' @keywords fitpercentiles
#' @export
#' @examples
#' q <- stats::qweibull(p = c(0.025, 0.5, 0.975), shape = 0.01, scale = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.weibull.par(q = q, scaleX = c(0.1, 0.03))
#' get.weibull.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.1, 0.99))
#' get.weibull.par(q = q, fit.weights = c(10, 1, 10))
#' get.weibull.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.1, 0.03))
#' get.weibull.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.1, 0.03))
#' graphics::par(old.par)
#'
#' q <- stats::qweibull(p = c(0.025, 0.5, 0.975), shape = 0.1, scale = 0.1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.weibull.par(q = q, scaleX = c(0.1, 0.05))
#' get.weibull.par(q = q, fit.weights = c(100, 1, 100), scaleX = c(0.00000001, 0.99999999999))
#' get.weibull.par(q = q, fit.weights = c(10, 1, 10), scaleX = c(0.00000001, 0.99999999999))
#' get.weibull.par(q = q, fit.weights = c(1, 100, 1), scaleX = c(0.00000001, 0.01))
#' get.weibull.par(q = q, fit.weights = c(1, 10, 1), scaleX = c(0.00000001, 0.1))
#' graphics::par(old.par)
#'
#' q <- stats::qweibull(p = c(0.025, 0.5, 0.975), shape = 2, scale = 3)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.weibull.par(q = q)
#' get.weibull.par(q = q, fit.weights = c(100, 1, 100))
#' get.weibull.par(q = q, fit.weights = c(10, 1, 10))
#' get.weibull.par(q = q, fit.weights = c(1, 100, 1))
#' get.weibull.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' q <- stats::qweibull(p = c(0.025, 0.5, 0.975), shape = 1, scale = 1)
#' old.par <- graphics::par(mfrow = c(2, 3))
#' get.weibull.par(q = q)
#' get.weibull.par(q = q, fit.weights = c(100, 1, 100))
#' get.weibull.par(q = q, fit.weights = c(10, 1, 10))
#' get.weibull.par(q = q, fit.weights = c(1, 100, 1))
#' get.weibull.par(q = q, fit.weights = c(1, 10, 1))
#' graphics::par(old.par)
#'
#' ## example with only two quantiles
#' q <- stats::qweibull(p = c(0.025, 0.975), shape = 2, scale = 1)
#' old.par <- graphics::par(mfrow = c(1, 3))
#' get.weibull.par(p = c(0.025, 0.975), q = q)
#' get.weibull.par(p = c(0.025, 0.975), q = q, fit.weights = c(100, 1))
#' get.weibull.par(p = c(0.025, 0.975), q = q, fit.weights = c(10, 1))
#' graphics::par(old.par)
#'
get.weibull.par <- function(p = c(0.025, 0.5, 0.975), q, 
                            show.output = TRUE, plot = TRUE, 
                            tol = 0.001, fit.weights = rep(1, length(p)), 
                            scaleX = c(0.1, 0.9), ...) {
    #-----------------------------------------------------------------------------
    # checking consistency of the input data
    #-----------------------------------------------------------------------------
    if (!is.numeric(p) | !is.numeric(q) | !is.numeric(fit.weights)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, not numerical items in the input vectors 'p', 'q' and/or 'fit.weights'!", call. = FALSE)
    }
    if (prod(order(p) == seq(1:length(p))) == 0 | prod(order(q) == seq(1:length(q))) == 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the vector of probabilities/percentiles is not ordered!", call. = FALSE)
    }
    if (min(p) < 0 | max(p) > 1) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, items of the probability vector should lie between 0 and 1!", call. = FALSE)
    }
    if (min(q) < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, percentiles are out of the domain [0, inf) => Weibull distribution couldn't be fitted!", call. = FALSE)
    }
    if (length(p) != length(q) | length(p) != length(fit.weights) | length(q) != length(fit.weights) ) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, 'p', 'q' and 'fit.weights' are not of the same length! The vectors of quantiles, probabilities and weightings should be of the same length.", call. = FALSE)
    }
    if (length(q) < 2) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, at least two quantiles must be known!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'show.output' should be logical!", call. = FALSE)
    }
    if (!is.logical(plot)) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'plot' should be logical!", call. = FALSE)
    }
    if (!is.numeric(tol) | length(tol) != 1 | tol < 0) {
        # on.exit(return(invisible(NA)))
        stop("INVALID INPUT, the argument 'tol' should be a single positive numerical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # minimizing procedure
    #-----------------------------------------------------------------------------
    fit.weights.original <- fit.weights
    fit.weights <- fit.weights/sum(fit.weights)
    minimize <- function(theta) {
        summand <- suppressWarnings(stats::pweibull(q = q, 
                                                    shape = theta[1], 
                                                    scale = theta[2]) - p)
        summand <- summand * fit.weights
        sum(summand^2)
    }
    fit <- c(); fit$value <- tol + 1
    try1 <- try(
        fit <- stats::optim(par = c(0, 1), 
                            minimize, method = "L-BFGS-B", 
                            lower = c(0.001, 0.001), 
                            upper = c(10000, 10000)), 
        silent = TRUE
    )
    
    #-----------------------------------------------------------------------------
    # checking results
    #-----------------------------------------------------------------------------
    if (is.error(try1) || fit$value >= tol) {
        warning("The fitting procedure 'L-BFGS-B' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE)
        fit <- c(); fit$value <- tol + 1
        try2 <- try(
            fit <- stats::optim(par = c(0, 1),       #1, 
                                minimize, 
                                method = "BFGS"), 
            silent = TRUE)
        if (is.error(try2) || fit$value >= tol) { 
            warning("The fitting procedure 'BFGS' has failed (convergence error occurred or specified tolerance not achieved)!", call. = FALSE) 
            Par <- NA
        } else if (fit$value < tol) {
            message("The fitting procedure 'BFGS' was successful!\n(Used this fallback optimization method because 'L-BFGS-B' has failed...)") 
            Par <- fit$par
            names(Par) <- c("shape", "scale")
            if (show.output) print(fit) 
        }
    } else if (fit$value < tol) {
        message("The fitting procedure 'L-BFGS-B' was successful!") 
        Par <- fit$par
        names(Par) <- c("shape", "scale")
        if (show.output) print(fit) 
    }
    
    #-----------------------------------------------------------------------------
    # plotting graphical diagnostics
    #-----------------------------------------------------------------------------
    if (prod(!is.na(Par)) & plot) {
        main1 <- paste("shape = ", round(Par["shape"], digits = 2))
        main2 <- paste("scale = ", round(Par["scale"], digits = 2))
        main <- paste("Weibull (", main1, ", ", main2, ")", sep = "")
        sub = paste("fit.weights = c(", paste(fit.weights.original, collapse = ", "), ")", sep = "")
        Support.lim <- c(stats::qweibull(p = min(p) * scaleX[1], 
                                         shape = Par["shape"], 
                                         scale = Par["scale"]),
                         stats::qweibull(p = (max(p) + (1 - max(p)) * scaleX[2]), 
                                         shape = Par["shape"], 
                                         scale = Par["scale"]))
        Support <- seq(min(min(q), Support.lim[1]), 
                       max(max(q), Support.lim[2]), 
                       length = 200)
        Probability <- stats::pweibull(Support, Par["shape"], Par["scale"])
        graphics::plot(Support, Probability, type = "l", 
                       xlim = range(Support.lim, q), main = main, 
                       xlab = "Quantiles", sub = sub, ...)
        graphics::points(x = q, y = p, pch = 19, ...)
    }
    
    #-----------------------------------------------------------------------------
    # output
    #-----------------------------------------------------------------------------
    return(Par)
} # end of get.weibull.par()







#*******************************************************************************
#*******************************************************************************
# FUNCTIONS FOT FITTING DISTRIBUTIONS TO DATA
#*******************************************************************************
#*******************************************************************************



################################################################################
################################################################################
#' Fitting amount continuous distributions to given univariate data
#'
#' This function is not intended to be called directly but is internally called
#' in \code{fit.cont}.
#'
#' @name useFitdist
#' @aliases useFitdist
#' @title Fitting amount continuous distributions to given univariate data.
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage useFitdist(data2fit, show.output = TRUE, distributions)
#' @param data2fit numerical vector, data to be fitted.
#' @param show.output logical value, if \code{TRUE} the output will be printed.
#' @param distributions simple character or character vector giving the names of
#' distribution families, that should be fitted to the data. The possible values
#' are: \code{norm}, \code{cauchy}, \code{logis}, \code{beta}, \code{exp},
#' \code{chisq}, \code{unif}, \code{gamma}, \code{lnorm}, \code{weibull},
#' \code{f}, \code{t}, \code{gompertz}, \code{triang}.
#' @return Returns matrix with fitting results. More information...
#' @keywords others
#' @export
#' @examples
#' x1 <- rgamma(374, 4,0.08)
#' res1 <- useFitdist(data2fit = x1)
#' res1
#'
#' x2 <- rbeta(300, shape1 = 1, shape2 = 2)
#' res2 <- useFitdist(data2fit = x2)
#' res2
#'
useFitdist <- function(data2fit, show.output = TRUE,
                       distributions = c("norm", "cauchy", "logis", "beta", "exp", 
                                         "chisq", "unif", "gamma", "lnorm", "weibull", 
                                         "f", "t", "gompertz", "triang")) {
    #-----------------------------------------------------------------------------
    # checking input
    #-----------------------------------------------------------------------------
    if (missing(data2fit)) {
        stop("Argument 'data2fit' ist empty!", call. = FALSE)
    }
    if (!is.null(dim(data2fit))) {
        stop("Argument 'data2fit' should be a simple numerical vector!", call. = FALSE)
    }
    if (!all(is.numeric(data2fit))) {
        stop("One or more items of 'data2fit' is/are not numeric!", call. = FALSE)
    }
    if (any(is.na(data2fit))) {
        stop("'data2fit' contains missing values, fitting procedure is not possible!", call. = FALSE)
    }
    if (!is.logical(show.output)) {
        stop("'show.output' should be a simple logical value!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # define variables
    #-----------------------------------------------------------------------------
    fit.list <- list()
    index <- 1
    
    #-----------------------------------------------------------------------------
    # internal function
    #-----------------------------------------------------------------------------
    trim.whitespace <- function(char) {
        char <- as.character(char)
        n <- length(char)
        out <- NULL
        for (i in 1:n) {
            c <- unlist(strsplit(char[i], split = " "))
            b <- which(c == "")
            if (length(b) > 0) c <- c[-b]
            c <- paste(c, collapse = " ")
            out <- c(out, c)
        }
        return(out)
    }
    
    #-----------------------------------------------------------------------------
    # fitting procedures
    #-----------------------------------------------------------------------------
    #     if (show.output) cat("\n-------------------------------------------------------------------\n") 
    #     if (show.output) cat("Begin fitting distributions... \n") 
    if (show.output) message("\nBegin fitting distributions ---------------------------------------")
    
    # fit normal distributions
    if (is.element("norm", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "norm")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting normal distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Normal"
            index <- index + 1
        } else if (show.output) message("* fitting normal distribution ... failed") 
    }
    
    # fit cauchy distribution
    if (is.element("cauchy", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "cauchy")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting Cauchy  distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Cauchy"
            index <- index + 1
        } else if (show.output) message("* fitting Cauchy distribution ... failed") 
    }
    
    # fit logistic distribution
    if (is.element("logis", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "logis")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting logistic distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Logistic"
            index <- index + 1
        } else if (show.output) message("* fitting logistic distribution ... failed") 
    }
    
    # fit beta distribution
    if (is.element("beta", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "beta")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting beta distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Beta"
            index <- index + 1
        } else if (show.output) message("* fitting beta distribution ... failed") 
    }
    
    # fit exponential distribution
    if (is.element("exp", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "exp")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting exponential distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Exponential"
            index <- index + 1
        } else if (show.output) message("* fitting exponential distribution ... failed") 
    }
    
    # fit exponential distribution
    if (is.element("chisq", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "chisq", start = mean(data2fit))), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting chi-square distribution ... OK")
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Chi-square"
            index <- index + 1
        } else if (show.output) message("* fitting chi-square distribution ... failed")
    }
    
    # fit uniform distribution
    if (is.element("unif", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "unif", method = "mme")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting uniform distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Uniform"
            index <- index + 1
        } else if (show.output) message("* fitting uniform distribution ... failed") 
    }
    
    # fit gamma distribution
    if (is.element("gamma", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "gamma")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting gamma distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Gamma"
            index <- index + 1
        } else if (show.output) message("* fitting gamma distribution ... failed") 
    }
    
    # fit lognormal distribution
    if (is.element("lnorm", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "lnorm")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting lognormal distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Lognormal"
            index <- index + 1
        } else if (show.output) message("* fitting lognormal distribution ... failed") 
    }
    
    # fit weibull distribution
    if (is.element("weibull", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "weibull")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting Weibull distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Weibull"
            index <- index + 1
        } else if (show.output) message("* fitting Weibull distribution ... failed") 
    }
    
    # fit F distribution
    if (is.element("f", distributions)) {if (mean(data2fit) != 1) {
        n <- abs(2 * mean(data2fit)/(mean(data2fit) - 1))
    } else if (mean(data2fit) == 1) {
        n <- 1
    }
        if ((stats::var(data2fit) * (n - 2)^2 * (n - 4) - 2 * n^2) == 0) {
            m <- 1
        } else if ((stats::var(data2fit) * (n - 2)^2 * (n - 4) - 2 * n^2) != 0) {
            m <- abs(2 * n^2 * (n - 2)/(stats::var(data2fit) * (n - 2)^2 * (n - 4) - 2 * n^2))
        }
        if (m == 0) m <- m + 0.001
        if (n == 0) n <- n + 0.001
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "f", start = c(m, n))), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting F-distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "F"
            index <- index + 1
        } else if (show.output) message("* fitting F-distribution ... failed") 
    }
    
    # fit t distribution
    if (is.element("t", distributions)) {
        if (stats::var(data2fit) != 1) {
            start.val <- abs(stats::var(data2fit) * 2/(stats::var(data2fit) - 1))
        } else if (stats::var(data2fit) == 1) {
            start.val <- 1
        }
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "t", start = start.val)), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting Student's t-distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Student"
            index <- index + 1
        } else if (show.output) message("* fitting Student's t-distribution ... failed") 
    }
    
    # fit gompertz distribution
    if (is.element("gompertz", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "gompertz")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting Gompertz distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Gompertz"
            index <- index + 1
        } else if (show.output) message("* fitting Gompertz distribution ... failed") 
    }
    
    # fit triangular distribution
    if (is.element("triang", distributions)) {
        try.result <- try(temp <- suppressWarnings(rriskFitdist.cont(data2fit, "triang")), silent = TRUE)
        if (!inherits(try.result, "try-error")) {
            if (show.output) message("* fitting triangular distribution ... OK") 
            fit.list[[index]] <- temp
            names(fit.list)[index] <- "Triangular"
            index <- index + 1
        } else if (show.output) message("* fitting triangular distribution ... failed") 
    }
    #     if (show.output) cat("End fitting distributions... \n") 
    #     if (show.output) cat("------------------------------------------------------------------- \n") 
    if (show.output) message("End fitting distributions -----------------------------------------\n")
    
    #-----------------------------------------------------------------------------
    # exit, if neither distribution could be fitted
    #-----------------------------------------------------------------------------
    if (length(fit.list) == 0) {
        on.exit(return(invisible(NULL)))
        stop("\n Neither continuous distribution could be fitted to the data \n", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # formatting results
    #-----------------------------------------------------------------------------
    temp <- lapply(fit.list, function(x) {
        result <- data.frame(
            ifelse(!is.null(x$loglik), round(x$loglik, digits = 2), "NULL"),
            ifelse(!is.null(x$aic), round(x$aic, digits = 2), "NULL"),
            ifelse(!is.null(x$bic), round(x$bic, digits = 2), "NULL"),
            ifelse(!is.null(x$chisq), round(x$chisq, digits = 2), "NULL"),
            ifelse(!is.null(x$chisqpvalue), round(x$chisqpvalue, digits = 2), "NULL"),
            ifelse(!is.null(x$ad), round(x$ad, digits = 2), "NULL"),
            ifelse(!is.null(x$adtest), x$adtest, "NULL"),
            ifelse(!is.null(x$ks), round(x$ks, digits = 2), "NULL"),
            ifelse(!is.null(x$kstest), x$kstest, "NULL"))  })
    
    for (i in 1:length(temp)) {
        if (i == 1) {
            res.matrix <- temp[[i]]
        } else {
            res.matrix <- rbind(res.matrix, temp[[i]])
        }
    }
    dimnames(res.matrix) <- list(names(fit.list), c("logL", "AIC", "BIC",
                                                    "Chisq(value)", "Chisq(p)", 
                                                    "AD(value)", "H(AD)", 
                                                    "KS(value)", "H(KS)"))
    res.matrix <- apply(res.matrix, c(1, 2), function(x) trim.whitespace(x))
    
    if (show.output) print(as.data.frame(res.matrix)) 
    
    #-----------------------------------------------------------------------------
    # create output
    #-----------------------------------------------------------------------------
    output <- list(res.matrix, fit.list)
    names(output) <- c("res.matrix", "fit.list")
    
    return(invisible(output))
} # end of useFitdist()



################################################################################
################################################################################
#' @description This function provides a GUI for choosing a most appropriate continuous
#' distribution fitted to given data.
#'
#' @name fit.cont
#' @aliases fit.cont
#' @title GUI for fitting continuous distributions to given data
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting)
#' @usage fit.cont(data2fit)
#' @param data2fit numerical vector, data to be fitted.
#' @return Returns chosen continuous distribution, estimated parameters and data,
#' on which the fitting is based.
#' @note This function is used for defining a Monte-Carlo random variate item
#' (\code{mcrv}) in the \code{rrisk} project.
#' @keywords gui
#' @export
#' @importFrom eha rgompertz
#' @importFrom eha dgompertz
#' @importFrom eha pgompertz
#' @importFrom eha qgompertz
#' @importFrom mc2d rtriang
#' @importFrom mc2d dtriang
#' @importFrom mc2d ptriang
#' @importFrom mc2d qtriang
#' @examples
#' \dontrun{
#'   if ( class(tcltk::tclRequire("Tktable")) == "tclObj" ) {
#'     res1 <- fit.cont(data2fit = rgamma(374, 4, 0.08))
#'     res1
#'
#'     res2 <- fit.cont(data2fit = rbeta(300, shape1 = 1, shape2 = 2))
#'     res2
#'
#'     res3 <- fit.cont(data2fit = mc2d::rtriang(300, min = 1, mode = 3, max = 10))
#'     res3
#'
#'     res4 <- fit.cont(data2fit = stats::rnorm(300))
#'     res4
#'   }
#' }
#' 
fit.cont <- function(data2fit = stats::rnorm(1000)) {
    if (class(tcltk::tclRequire("Tktable")) != "tclObj") {
        stop("Tcl package \"Tktable\" required. Please install it.")
    }
    
    #-----------------------------------------------------------------------------
    # checking input
    #-----------------------------------------------------------------------------
    if (!is.null(dim(data2fit))) {
        stop("Argument 'data2fit' should be a simple numerical vector!", call. = FALSE)
    }
    if (!all(is.numeric(data2fit))) {
        stop("One or more items of 'data2fit' is/are not numeric!", call. = FALSE)
    }
    if (any(is.na(data2fit))) {
        stop("'data2fit' contains missing values, fitting procedure is not possible!", call. = FALSE)
    }
    
    #-----------------------------------------------------------------------------
    # create fitting results matrix
    #-----------------------------------------------------------------------------
    useFitdist.results <- useFitdist(data2fit)
    
    if (is.null(useFitdist.results)) {
        # if (is.element("package:rrisk", search())) # wenn "rrisk" vorhanden, mache weiter. sonst breche ab.
        #{
        #  on.exit(.generate.newitem())
        #  utils::winDialog(type = "ok", "Neither continuous distribution could be fitted to the data")
        #  stop("exit fit.cont() and call generate.newitem()", call. = FALSE)
        #} else
        #{
        on.exit(return(invisible(NULL)))
        utils::winDialog(type = "ok", "Neither continuous distribution could be fitted to the data")
        stop("exit fit.cont() and return NA's", call. = FALSE)
        #}
    }
    
    res.matrix <- useFitdist.results$res.matrix
    fit.list <- useFitdist.results$fit.list
    
    notRejectedIndex <- which(res.matrix[, which(colnames(res.matrix) == "H(KS)")] == "not rejected")
    
    
    newIndizes <- c(notRejectedIndex, setdiff(1:nrow(res.matrix), notRejectedIndex))
    res.matrix <- res.matrix[newIndizes, ]
    #rbind(res.matrix[c(notRejectedIndex, setdiff()), ], res.matrix[-notRejectedIndex, ])
    #res.matrix <- rbind(res.matrix[notRejectedIndex, ], res.matrix[-notRejectedIndex, ])
    
    #-----------------------------------------------------------------------------
    # create tcltk::tclArray for results matrix
    #-----------------------------------------------------------------------------
    res.matrix <- cbind(dimnames(res.matrix)[[1]], res.matrix)
    res.matrix <- rbind(c(dimnames(res.matrix)[[2]]), res.matrix)
    res.matrix[1, 1] <- "Family"
    tclarray <- tcltk::tclArray()
    for (i in 0:(nrow(res.matrix) - 1)) {
        for (j in 0:(ncol(res.matrix) - 1)) {
            tclarray[[i, j]] <- res.matrix[i + 1, j + 1]
        }
    }
    
    #-----------------------------------------------------------------------------
    # create temporer environment
    #-----------------------------------------------------------------------------
    assign("tempEnvir", value = new.env())
    assign("chosenD", value = NA, envir = tempEnvir)
    assign("fittedParams", value = NA, envir = tempEnvir)
    #assign("randomCommand", value = NA, envir = tempEnvir)
    
    #-----------------------------------------------------------------------------
    # function for apdating diagnoctic plot
    #-----------------------------------------------------------------------------
    updatePlot <- function(...) {
        tkrplot::tkrreplot(imgPlot)
        tcltk::tkraise(fitContDistWindow)
    }
    
    #-----------------------------------------------------------------------------
    # function for plotting fitting diagnostics
    #-----------------------------------------------------------------------------
    plotDiagnostics <- function(...) {
        distr <- tcltk::tclvalue(rbValue)
        
        if (distr == "Normal") distname = "norm"
        if (distr == "Beta") distname = "beta"
        if (distr == "Cauchy") distname = "cauchy"
        if (distr == "Logistic") distname = "logis"
        if (distr == "Exponential") distname = "exp"
        if (distr == "Chi-square") distname = "chisq"
        if (distr == "Uniform") distname = "unif"
        if (distr == "Gamma") distname = "gamma"
        if (distr == "Lognormal") distname = "lnorm"
        if (distr == "Weibull") distname = "weibull"
        if (distr == "Student") distname = "t"
        if (distr == "F") distname = "f"
        if (distr == "Gompertz") distname = "gompertz"
        if (distr == "Triangular") distname = "triang"
        
        params <- fit.list[[which(names(fit.list) == distr)]]$estimate
        x <- seq(min(data2fit) - (max(data2fit) - min(data2fit)) * 0.05,
                 max(data2fit) + (max(data2fit) - min(data2fit)) * 0.05, 
                 length = length(data2fit))
        
        rdistname <- paste("r", distname, sep = "")
        ddistname <- paste("d", distname, sep = "")
        pdistname <- paste("p", distname, sep = "")
        qdistname <- paste("q", distname, sep = "")
        
        set.seed(1)
        y <- do.call(rdistname, c(list(n = 300), as.list(params)))
        d <- do.call(ddistname, c(list(x = x), as.list(params)))
        p <- do.call(pdistname, c(list(q = data2fit), as.list(params)))
        pp <- do.call(pdistname, c(list(q = x), as.list(params)))
        
        graphics::par(mfrow = c(1, 4), oma = c(0, 0,2, 0))
        title.text <- paste("Diagnostic plots for", distr, "distribution")
        graphics::hist(data2fit, 
                       main = "Emp. and theor. distributions", 
                       xlab = "data", ylab = "density",
                       probability = TRUE, cex.main = 1, 
                       ylim = c(0, max(d, graphics::hist(data2fit, plot = FALSE)$density)))
        graphics::lines(x, d, lwd = 2, col = "red")
        
        stats::qqplot(y, data2fit, 
                      main = "QQ-plot", 
                      xlab = "theoretical quantiles", ylab = "sample quantiles", 
                      cex.main = 1, pch = 20)
        graphics::abline(0, 1, col = "red", lwd = 2)
        
        graphics::plot(stats::ecdf(data2fit), 
                       main = "Empirical and theoretical CDFs", 
                       xlab = "data", ylab = "CDF", 
                       cex.main = 1, pch = 20)
        graphics::lines(x, pp, lwd = 2, col = "red")
        
        graphics::plot(sort(p), 
                       stats::ecdf(sort(data2fit))(sort(data2fit)), 
                       main = "PP-plot", 
                       xlab = "theoretical probabilities", ylab = "sample probabilities",
                       cex.main = 1, pch = 20)
        graphics::abline(0, 1, col = "red", lwd = 2)
        
        graphics::title(main = title.text, outer = TRUE, cex.main = 1.5)
    } # end of function plot.diagnostics()
    
    #-----------------------------------------------------------------------------
    # what to do on "Ok" button
    #-----------------------------------------------------------------------------
    onOK <- function(...) {
        assign("chosenD", value = tcltk::tclvalue(rbValue), envir = tempEnvir)
        fittedParams <- fit.list[[which(names(fit.list) == tcltk::tclvalue(rbValue))]]$estimate
        assign("fittedParams", value = fittedParams, envir = tempEnvir)
        tcltk::tkdestroy(fitContDistWindow)
    } # end of onOk()
    
    #-----------------------------------------------------------------------------
    # what to do on "Cancel" button
    #-----------------------------------------------------------------------------
    onCancel <- function(...) {
        tcltk::tkdestroy(fitContDistWindow)
    } # end of onCancel()
    
    #-----------------------------------------------------------------------------
    # create dialog window
    #-----------------------------------------------------------------------------
    # fitContDistWindow <- tcltk::tktoplevel(width = 880, height = 200)
    fitContDistWindow <- tcltk::tktoplevel(width = 880, height = 200)
    tcltk::tkwm.title(fitContDistWindow, "Fitting continuous distributions")
    tcltk::tkwm.resizable(fitContDistWindow, 1, 1)  # fixed size, not resizeable
    tcltk::tkwm.maxsize(fitContDistWindow, 1000, 800)
    tcltk::tkwm.minsize(fitContDistWindow, 880, 200)
    allFrame <- tcltk::tkframe(fitContDistWindow)
    fontA <- tcltk::tkfont.create(size = 12, weight = "bold")
    index = 0
    rbValue <- tcltk::tclVar(dimnames(res.matrix)[[1]][2])
    
    # creat image frame
    imageFrame <- tcltk::tkframe(allFrame, relief = "groove", background = "white")
    imgPlot <- tkrplot::tkrplot(imageFrame, fun = plotDiagnostics, 
                                hscale = 2.2, vscale = 0.7)
    tcltk::tkpack(imgPlot)
    tcltk::tkpack(imageFrame)
    
    # create table and buttons frame
    tableButtonsradioFrame <- tcltk::tkframe(allFrame)
    tableButtonsFrame <- tcltk::tkframe(tableButtonsradioFrame)
    
    # crate table
    fitResultTable <- tcltk::tkwidget(
        tableButtonsFrame, "table", variable = tclarray,
        rows = nrow(res.matrix), cols = ncol(res.matrix), 
        background = "white", borderwidth = 2, state = "disabled", 
        titlerows = 1, titlecols = 1, 
        resizeborders = "none", colwidth = 11)
    tcltk::tkpack(fitResultTable, side = "top")
    
    # create buttons frame
    buttonsFrame <- tcltk::tkframe(tableButtonsFrame)
    okButton <- tcltk::ttkbutton(buttonsFrame, width = 10, 
                                 text = "OK", command = onOK)
    CancelButton <- tcltk::ttkbutton(buttonsFrame, width = 10, 
                                     text = "Cancel", command = onCancel)
    tcltk::tkpack(okButton, side = "left", padx = c(20, 20))
    tcltk::tkpack(CancelButton, side = "left", padx = c(20, 20))
    tcltk::tkpack(buttonsFrame, pady = c(12, 0), side = "bottom")
    
    tcltk::tkpack(tableButtonsFrame, side = "right")
    
    # creat radio buttons frame
    radioFrame <- tcltk::tkframe(tableButtonsradioFrame, height = 5)
    rbLabel <- tcltk::tklabel(radioFrame, text = "Distribution", 
                              font = fontA, width = 13)
    tcltk::tkgrid(rbLabel, column = 0, row = index, sticky = "w")
    if (is.element("Normal", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Normal") - 1
        rbNormal <- tcltk::tkradiobutton(
            radioFrame, variable = rbValue, 
            value = "Normal", text = "Normal",
            command = function() updatePlot(distname = "norm"))
        tcltk::tkgrid(rbNormal, column = 0, row = index, 
                      sticky = "w", pady = 0.1)
    }
    if (is.element("Beta", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Beta") - 1
        rbBeta <- tcltk::tkradiobutton(
            radioFrame, value = "Beta", text = "Beta", 
            variable = rbValue,
            command = function() updatePlot(distname = "beta"))
        tcltk::tkgrid(rbBeta, column = 0, row = index, 
                      sticky = "w", pady = 0.)
    }
    if (is.element("Cauchy", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Cauchy") - 1
        rbCauchy <- tcltk::tkradiobutton(
            radioFrame, value = "Cauchy", text = "Cauchy", 
            variable = rbValue,
            command = function() updatePlot(distname = "cauchy"))
        tcltk::tkgrid(rbCauchy, column = 0, row = index, sticky = "w")
    }
    if (is.element("Logistic", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Logistic") - 1
        rbLogistic <- tcltk::tkradiobutton(
            radioFrame, value = "Logistic", 
            text = "Logistic", variable = rbValue,
            command = function() updatePlot(distname = "logis"))
        tcltk::tkgrid(rbLogistic, column = 0, row = index, sticky = "w")
    }
    if (is.element("Exponential", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Exponential") - 1
        rbExp <- tcltk::tkradiobutton(
            radioFrame, value = "Exponential", 
            text = "Exponential", variable = rbValue,
            command = function() updatePlot(distname = "exp"))
        tcltk::tkgrid(rbExp, column = 0, row = index, sticky = "w")
    }
    if (is.element("Chi-square", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Chi-square") - 1
        rbChi2 <- tcltk::tkradiobutton(
            radioFrame, value = "Chi-square", 
            text = "Chi-square", variable = rbValue,
            command = function() updatePlot(distname = "chisq"))
        tcltk::tkgrid(rbChi2, column = 0, row = index, sticky = "w")
    }
    if (is.element("Uniform", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Uniform") - 1
        rbUnif <- tcltk::tkradiobutton(
            radioFrame, value = "Uniform", 
            text = "Uniform", variable = rbValue,
            command = function() updatePlot(distname = "unif"))
        tcltk::tkgrid(rbUnif, column = 0, row = index, sticky = "w")
    }
    if (is.element("Gamma", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Gamma") - 1
        rbGamma <- tcltk::tkradiobutton(
            radioFrame, value = "Gamma", 
            text = "Gamma", variable = rbValue,
            command = function() updatePlot(distname = "gamma"))
        tcltk::tkgrid(rbGamma, column = 0, row = index, sticky = "w")
    }
    if (is.element("Lognormal", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Lognormal") - 1
        rbLognorm <- tcltk::tkradiobutton(
            radioFrame, value = "Lognormal", 
            text = "Lognormal", variable = rbValue,
            command = function() updatePlot(distname = "lnorm"))
        tcltk::tkgrid(rbLognorm, column = 0, row = index, sticky = "w")
    }
    if (is.element("Weibull", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Weibull") - 1
        rbWeibull <- tcltk::tkradiobutton(
            radioFrame, value = "Weibull", 
            text = "Weibull", variable = rbValue,
            command = function() updatePlot(distname = "weibull"))
        tcltk::tkgrid(rbWeibull, column = 0, row = index, sticky = "w")
    }
    if (is.element("F", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "F") - 1
        rbF <- tcltk::tkradiobutton(
            radioFrame, value = "F", 
            text = "F", variable = rbValue,
            command = function() updatePlot(distname = "f"))
        tcltk::tkgrid(rbF, column = 0, row = index, sticky = "w")
    }
    if (is.element("Student", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Student") - 1
        rbStudent <- tcltk::tkradiobutton(
            radioFrame, value = "Student", 
            text = "Student's t", variable = rbValue,
            command = function() updatePlot(distname = "t"))
        tcltk::tkgrid(rbStudent, column = 0, row = index, sticky = "w")
    }
    if (is.element("Gompertz", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Gompertz") - 1
        rbGompertz <- tcltk::tkradiobutton(
            radioFrame, value = "Gompertz", 
            text = "Gompertz", variable = rbValue,
            command = function() updatePlot(distname = "gompertz"))
        tcltk::tkgrid(rbGompertz, column = 0, row = index, sticky = "w")
    }
    if (is.element("Triangular", dimnames(res.matrix)[[1]])) {
        index <- which(rownames(res.matrix) == "Triangular") - 1
        #index <- index + 1
        rbTriangular <- tcltk::tkradiobutton(
            radioFrame, value = "Triangular", 
            text = "Triangular", variable = rbValue,
            command = function() updatePlot(distname = "triang"))
        tcltk::tkgrid(rbTriangular, column = 0, row = index, sticky = "w")
    }
    tcltk::tkpack(radioFrame, side = "left")
    tcltk::tkpack(tableButtonsradioFrame)
    tcltk::tkpack(allFrame, padx = c(8, 8), pady = c(8, 8))
    
    tcltk::tkfocus(fitContDistWindow)
    tcltk::tkwait.window(fitContDistWindow)
    
    #-----------------------------------------------------------------------------
    # generate output
    #-----------------------------------------------------------------------------
    chosenD <- get("chosenD", envir = tempEnvir)
    fittedParams <- get("fittedParams", envir = tempEnvir)
    exitMessage <- "NA"
    if (!is.na(chosenD)) {
        if (chosenD == "Triangular") {
            chosenD <- "triang"
            exitMessage <- "Triangular (triang)"
        } else if (chosenD == "Gompertz") {
            chosenD <- "gompertz"
            exitMessage <- "Gompertz (gompertz)"
        } else if (chosenD == "Student") {
            chosenD <- "t"
            exitMessage <- "Student's t (t)"
        } else if (chosenD == "F") {
            chosenD <- "f"
            exitMessage <- "F (f)"
        } else if (chosenD == "Weibull") {
            chosenD <- "weibull"
            exitMessage <- "Weibull (weibull)"
        } else if (chosenD == "Lognormal") {
            chosenD <- "lnorm"
            exitMessage <- "Log-normal (lnorm)"
        } else if (chosenD == "Gamma") {
            chosenD <- "gamma"
            exitMessage <- "Gamma (gamma)"
        } else if (chosenD == "Uniform") {
            chosenD <- "unif"
            exitMessage <- "Uniform (unif)"
        } else if (chosenD == "Chi-square") {
            chosenD <- "chisq"
            exitMessage <- "Chi-squared (chisq)"
        } else if (chosenD == "Exponential") {
            chosenD <- "exp"
            exitMessage <- "Exponential (exp)"
        } else if (chosenD == "Logistic") {
            chosenD <- "logis"
            exitMessage <- "Logistic (logis)"
        } else if (chosenD == "Cauchy") {
            chosenD <- "cauchy"
            exitMessage <- "Cauchy (cauchy)"
        } else if (chosenD == "Beta") {
            chosenD <- "beta"
            exitMessage <- "Beta (beta)"
        } else if (chosenD == "Normal") {
            chosenD <- "norm"
            exitMessage <- "Normal (norm)"
        }
    }
    output <- list(data2fit, chosenD, fittedParams)
    names(output) <- c("data2fit", "chosenDistr", "fittedParams")
    
    print.on.exit <- function(chosenD) {
        cat(paste("\nChosen continuous distribution is:", exitMessage))
        cat("\nFitted parameters are: \n")
        print(fittedParams)
        cat("\n")
    } # end of fucntion print.on.exit()
    
    on.exit(print.on.exit(chosenD))
    return(invisible(output))
}

################################################################################
################################################################################
#' Fits a univariate distribution by matching moments.
#'
#' This function is an alias of the function \code{mmedist} from the package
#' \pkg{fitdistrplus} (Version 0.1-2). The original function was extended to
#' fitting additional distributions. Parameter of the following distribution
#' families can be estimated in this function: normal, lognormal, exponential,
#' Poisson, gamma, logistic, negative binomial, geometric, Beta and continuous
#' univariate.
#' \cr \cr For more details see the assistance page of the function
#' \code{mmedist} from the package \pkg{fitdistrplus}.
#' \cr \cr This function is not intended to be called directly but is internally
#' called in \code{rriskFitdist.cont}.
#'
#' @name rriskMMEdist
#' @aliases rriskMMEdist
#' @title Fitting univariate distributions by matching moments
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Marie-Laure Delignette-Muller (coauthor of the package \pkg{fitdistrplus}), \cr
#' Christophe Dutang (coauthor of the package \pkg{fitdistrplus})
#' @usage rriskMMEdist(data, distr)
#' @param data a numerical vector.
#' @param distr A character string "name" naming a distribution or directly the
#' density function \code{dname}. The estimated values of the distribution parameters
#' are provided only for the following distributions : "norm", "lnorm", "exp",
#' "pois", "gamma", "logis", "nbinom" , "geom", "beta" and "unif".
#' @return \code{rriskMMEdist} returns the named parameter or a named vector of parameters.
#' @keywords fitdistrplus
#' @export
#' @examples
#' ## Continuous distributions
#' set.seed(1)
#' x1 <- stats::rnorm(500, mean = 2, sd = 0.7)
#' rriskMMEdist(x1, "norm")
#' rriskMMEdist(x1, "exp")
#' rriskMMEdist(x1, "gamma")
#' rriskMMEdist(x1, "logis")
#' rriskMMEdist(x1, "unif")
#' 
#' ## produces an error:
#' # rriskMMEdist(x1, "lnorm")
#' # rriskMMEdist(x1, "beta")
#'
#' ## Discrete distributions
#' set.seed(2)
#' x2 <- rpois(500, lambda = 3)
#' rriskMMEdist(x2, "pois")
#' rriskMMEdist(x2, "nbinom")
#' rriskMMEdist(x2, "geom")
#'
rriskMMEdist <- function(data, distr) {
    if (!is.character(distr)) {
        distname <- substring(as.character(match.call()$distr), 2)
    } else distname <- distr
    if (!is.element(distname, c("norm", "lnorm", "pois", "exp",
                                "gamma", "nbinom", "geom", "beta", 
                                "unif", "logis"))) {
        stop(paste("Method of matching moments is not available for ",
                   distname, " distribution"))
    }
    if (!(is.vector(data) & is.numeric(data) & length(data) > 1)) {
        stop("data must be a numerical vector of length greater than 1")
    }
    if (distname == "norm") {
        n <- length(data)
        sd0 <- sqrt((n - 1)/n) * stats::sd(data)
        mx <- mean(data)
        estimate <- c(mx, sd0)
        names(estimate) <- c("mean", "sd")
        return(estimate)
    }
    if (distname == "lnorm") {
        if (any(data <= 0)) {
            stop("values must be positive to fit a lognormal distribution")
        }
        n <- length(data)
        ldata <- log(data)
        sd0 <- sqrt((n - 1)/n) * stats::sd(ldata)
        ml <- mean(ldata)
        estimate <- c(ml, sd0)
        names(estimate) <- c("meanlog", "sdlog")
        return(estimate)
    }
    if (distname == "pois") {
        estimate <- mean(data)
        names(estimate) <- "lambda"
        return(estimate)
    }
    if (distname == "exp") {
        estimate <- 1/mean(data)
        names(estimate) <- "rate"
        return(estimate)
    }
    if (distname == "gamma") {
        n <- length(data)
        m <- mean(data)
        v <- (n - 1)/n * stats::var(data)
        shape <- m^2/v
        rate <- m/v
        estimate <- c(shape, rate)
        names(estimate) <- c("shape", "rate")
        return(estimate)
    }
    if (distname == "nbinom") {
        n <- length(data)
        m <- mean(data)
        v <- (n - 1)/n * stats::var(data)
        size <- ifelse(v > m, m^2/(v - m), NaN)
        estimate <- c(size, m)
        names(estimate) <- c("size", "mu")
        return(estimate)
    }
    if (distname == "geom") {
        m <- mean(data)
        prob <- ifelse(m > 0, 1/(1 + m), NaN)
        estimate <- prob
        names(estimate) <- "prob"
        return(estimate)
    }
    if (distname == "beta") {
        if (any(data < 0) | any(data > 1)) {
            stop("values must be in [0-1] to fit a beta distribution")
        }
        n <- length(data)
        m <- mean(data)
        v <- (n - 1)/n * stats::var(data)
        aux <- m * (1 - m)/v - 1
        shape1 <- m * aux
        shape2 <- (1 - m) * aux
        estimate <- c(shape1, shape2)
        names(estimate) <- c("shape1", "shape2")
        return(estimate)
    }
    if (distname == "unif") {
        n <- length(data)
        m <- mean(data)
        v <- (n - 1)/n * stats::var(data)
        min1 <- m - sqrt(3 * v)
        max1 <- m + sqrt(3 * v)
        estimate <- c(min1, max1)
        names(estimate) <- c("min", "max")
        return(estimate)
    }
    if (distname == "logis") {
        n <- length(data)
        m <- mean(data)
        v <- (n - 1)/n * stats::var(data)
        scale <- sqrt(3 * v)/pi
        estimate <- c(m, scale)
        names(estimate) <- c("location", "scale")
        return(estimate)
    }
}

################################################################################
################################################################################
#' Fits a univariate distribution by maximum likelihood.
#'
#' This function is an alias of the function \code{mledist} from the package
#' \pkg{fitdistrplus} (Version 0.1-2). The original function was extended to
#' fitting additional distributions. The following continuous distributions can
#' be fitted by this function: normal, exponential, lognormal, logistic, gamma,
#' Weibull, Beta, chi-square, Student's t, F, Cauchy, Gompertz and triangular.
#' And the following discrete distributions can be fitted: (wird ergaenzt).
#' \cr \cr For more details see the assistance page of the function
#' \code{mledist} from the package \pkg{fitdistrplus}.
#' \cr \cr This function is not intended to be called directly but is internally
#' called in \code{rriskFitdist.cont}.
#'
#' @name rriskMLEdist
#' @aliases rriskMLEdist
#' @title Maximum likelihood fitting of univariate distributions
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Marie-Laure Delignette-Muller (coauthor of the package \pkg{fitdistrplus}), \cr
#' Christophe Dutang (coauthor of the package \pkg{fitdistrplus})
#' @usage rriskMLEdist(data, distr, start = NULL, optim.method = "default",
#'    lower = -Inf, upper = Inf, custom.optim = NULL, ...)
#' @param data A numerical vector for non censored data or a dataframe of two
#' columns respectively named left and right, describing each observed value as
#' an interval for censored data. In that case the left column contains either
#' NA for left censored observations, the left bound of the interval for interval
#' censored observations, or the observed value for non-censored observations.
#' The right column contains either NA for right censored observations, the right
#' bound of the interval for interval censored observations, or the observed value
#' for non-censored observations.
#' @param distr A character string "name" naming a distribution (or directly the
#' density function) for which the corresponding density function dname and the
#' corresponding distribution pname must be classically defined. The possible values are:
#' "norm", "exp", "lnorm", "logis", "gamma", "weibull", "beta", "chisq", "t", "f",
#' "cauchy", "gompertz".
#' @param start A named list giving the initial values of parameters of the named
#' distribution. This argument may be omitted for some distributions for which
#' reasonable starting values are computed (see details).
#' @param optim.method	"default" (see details) or optimization method to pass to \code{optim}.
#' @param lower	Left bounds on the parameters for the "L-BFGS-B" method (see \code{optim}).
#' @param upper Right bounds on the parameters for the "L-BFGS-B" method (see \code{optim}).
#' @param custom.optim a function carrying the MLE optimization (see details).
#' @param ...	further arguments passed to the \code{optim} or \code{custom.optim} function.
#' @return \code{rriskMLEdist} returns a list with fitting results containing following informations
#' \item{\code{estimate}}{numeric, a single value or a vector containing estimated parameters.}
#' \item{\code{convergence}}{an integer code for the convergence of \code{optim}.
#' The value \code{0} indicates a successful convergence.}
#' \item{\code{loglik}}{the log-likelihood}
#' \item{\code{hessian}}{a symmetric matrix computed by \code{optim} as an estimate of the Hessian
#' at the solution found or computed in the user-supplied optimization function.
#' It is used in \code{rriskFitdist.cont} to estimate standard errors.}
#' \item{\code{optim.function}}{the name of the optimization function used for maximum likelihood.}
#' @keywords fitdistrplus
#' @export
# @importFrom mc2d dtriang
# @importFrom mc2d ptriang
#' @importFrom eha dgompertz
#' @importFrom eha pgompertz
#' @examples
#' ## a basic fit of some distribution with maximum likelihood estimation
#' set.seed(1)
#' x2 <- rchisq(500, 4)
#' rriskMLEdist(x2, "norm")
#' rriskMLEdist(x2, "exp")
#' rriskMLEdist(x2, "lnorm")
#' rriskMLEdist(x2, "logis")
#' rriskMLEdist(x2, "gamma")
#' rriskMLEdist(x2, "weibull")
#' #rriskMLEdist(x2, "beta")
#' rriskMLEdist(x2, "chisq")
#' rriskMLEdist(x2, "t")
#' rriskMLEdist(x2, "f")
#' rriskMLEdist(x2, "cauchy")
#' rriskMLEdist(x2, "gompertz")
#' 
#' ## produces an error:
#' # rriskMLEdist(x2, "triang")
#' 
rriskMLEdist <- function(data, distr, 
                         start = NULL, optim.method = "default",
                         lower = -Inf, upper = Inf, 
                         custom.optim = NULL, ...) {
    distname <- ifelse(!is.character(distr), 
                       substring(as.character(match.call()$distr), 2),
                       distr)
    ddistname <- paste("d", distname, sep = "")
    #pdistname<- paste("p", distname, sep = "")  #????
    if (!exists(ddistname, mode = "function")) {
        stop(paste("The", ddistname, "function must be defined"))
    }
    if (distname == "unif") {
        stop("Maximum likelihood estimation is not available for the uniform distribution")
    }
    if (is.vector(data)) {
        cens <- FALSE
        if (!(is.numeric(data) & length(data) > 1)) {
            stop("data must be a numerical vector of length greater than 1 for non censored data\n            or a data frame with two columns named left and right and more than one line for censored data")
        }
    } else {
        cens <- TRUE
        censdata <- data
        if (!(is.vector(censdata$left) & 
              is.vector(censdata$right) &
              length(censdata[, 1]) > 1)) {
            stop("data must be a numerical vector of length greater than 1 for non censored data\n        or a dataframe with two columns named left and right and more than one line for censored data")
        }
        pdistname <- paste("p", distname, sep = "")
        if (!exists(pdistname, mode = "function")) {
            stop(paste("The ", pdistname, " function must be defined to apply maximum likelihood to censored data"))
        }
    }
    if (cens) {
        lcens <- censdata[is.na(censdata$left), ]$right
        if (any(is.na(lcens))) {
            stop("An observation cannot be both right and left censored, coded with two NA values")
        }
        rcens <- censdata[is.na(censdata$right), ]$left
        ncens <- censdata[censdata$left == censdata$right & 
                              !is.na(censdata$left) &
                              !is.na(censdata$right), ]$left
        icens <- censdata[censdata$left != censdata$right & 
                              !is.na(censdata$left) &
                              !is.na(censdata$right), ]
        data <- c(rcens, lcens, ncens, (icens$left + icens$right)/2)
    }
    if (is.null(start)) {
        if (distname == "norm") {
            n <- length(data)
            sd0 <- sqrt((n - 1)/n) * stats::sd(data)
            mx <- mean(data)
            start <- list(mean = mx, sd = sd0)
        } else if (distname == "lnorm") {
            if (any(data <= 0)) {
                stop("values must be positive to fit a lognormal distribution")
            }
            n <- length(data)
            ldata <- log(data)
            sd0 <- sqrt((n - 1)/n) * stats::sd(ldata)
            ml <- mean(ldata)
            start <- list(meanlog = ml, sdlog = sd0)
        } else if (distname == "pois") {
            start <- list(lambda = mean(data))
        } else if (distname == "exp") {
            if (any(data < 0)) {
                stop("values must be positive to fit a exponential distribution")
            }
            start <- list(rate = 1/mean(data))
        } else if (distname == "gamma") {
            if (any(data < 0)) {
                stop("values must be positive to fit a Gamma distribution")
            }
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n * stats::var(data)
            start <- list(shape = m^2/v, rate = m/v)
        } else if (distname == "nbinom") {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n * stats::var(data)
            size <- ifelse(v > m, m^2/(v - m), 100)
            start <- list(size = size, mu = m)
        } else if (distname == "geom") {
            m <- mean(data)
            prob <- ifelse(m > 0, 1/(1 + m), 1)
            start <- list(prob = prob)
        } else if (distname == "t") {
            df.start <- 2 * stats::sd(data)^2/(stats::sd(data)^2 - 1)
            start <- list(df = df.start)
        } else if (distname == "beta") {
            if (any(data < 0) | any(data > 1)) {
                stop("values must be in [0-1] to fit a beta distribution")
            }
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n * stats::var(data)
            aux <- m * (1 - m)/v - 1
            start <- list(shape1 = m * aux, shape2 = (1 - m) *
                              aux)
        } else if (distname == "weibull") {
            if (any(data < 0)) {
                stop("values must be positive to fit a Weibull distribution")
            }
            m <- mean(log(data))
            v <- stats::var(log(data))
            shape <- 1.2/sqrt(v)
            scale <- exp(m + 0.572/shape)
            start <- list(shape = shape, scale = scale)
        } else if (distname == "logis") {
            n <- length(data)
            m <- mean(data)
            v <- (n - 1)/n * stats::var(data)
            start <- list(location = m, scale = sqrt(3 * v)/pi)
        } else if (distname == "chisq") {
            if (any(data < 0)) {
                stop("values must be positive to fit a Chi-square distribution")
            }
            start <- list(df = mean(data))
        } else if (distname == "f") {
            if (any(data < 0)) {
                stop("values must be positive to fit a F distribution")
            }
            df2.start <- 2 * mean(data)/(mean(data) - 1)
            start <- list(df1 = 3, df2 = df2.start)
        } else if (distname == "gompertz") {
            if (any(data < 0)) {
                stop("values must be positive to fit a Gompertz distribution")
            }
            scale.start <- stats::sd(data) * sqrt(6)/pi
            gamma.const <- 0.577221566
            shape.start <- (scale.start * exp(mean(data)/scale.start + gamma.const))^(-1)
            start <- list(shape = shape.start, scale = scale.start)
        } else if (distname == "cauchy") {
            start <- list(location = stats::median(data), scale = stats::IQR(data)/2)
        } else stop("Fitting procedure for given distribution is not implemented", 
                    call. = FALSE)
        #if (distname == "triang") {
        #    start <- list(min = min(data), mode = (max(data) - min(data))/2, max = max(data))
        #}
        if (!is.list(start)) {
            stop("'start' must be defined as a named list for this distribution")
        }
    }
    vstart <- unlist(start)
    argddistname <- names(formals(ddistname))
    m <- match(names(start), argddistname)
    if (any(is.na(m))) {
        stop("'start' must specify names which are arguments to 'distr'")
    }
    if (!cens) {
        if ("log" %in% argddistname) {
            fnobj <- function(par, obs, ddistnam) {
                -sum(do.call(ddistnam, c(list(obs), par, log = TRUE)))
            }
        } else {
            fnobj <- function(par, obs, ddistnam) {
                -sum(log(do.call(ddistnam, c(list(obs), par))))
            }
        }
    } else {
        argpdistname <- names(formals(pdistname))
        if (("log" %in% argddistname) & ("log.p" %in% argpdistname)) {
            fnobjcens <- function(par, rcens, lcens, icens, ncens, ddistnam, pdistnam) {
                - sum(do.call(ddistnam, c(list(x = ncens), as.list(par), list(log = TRUE))))
                - sum(do.call(pdistnam, c(list(q = lcens), as.list(par), list(log = TRUE))))
                - sum(do.call(pdistnam, c(list(q = rcens), as.list(par), list(lower.tail = FALSE), list(log = TRUE))))
                - sum(log(do.call(pdistnam, c(list(q = icens$right), as.list(par))) 
                          - do.call(pdistnam, c(list(q = icens$left), as.list(par)))))
            }
        } else {
            fnobjcens <- function(par, rcens, lcens, icens, ncens, ddistnam, pdistnam) {
                - sum(log(do.call(ddistnam, c(list(x = ncens), as.list(par))))) 
                - sum(log(do.call(pdistnam, c(list(q = lcens), as.list(par)))))
                - sum(log(1 - do.call(pdistnam, c(list(q = rcens), as.list(par)))))
                - sum(log(do.call(pdistnam, c(list(q = icens$right), as.list(par))) 
                          - do.call(pdistnam, c(list(q = icens$left), as.list(par)))))
            }
        }
    }
    if (optim.method == "default") {
        if (length(vstart) > 1) {
            meth <- "Nelder-Mead"
            if (ddistname == "dtriang") meth = "SANN"
        } else meth <- "BFGS"
    } else meth = optim.method
    if (is.null(custom.optim)) {
        if (!cens) {
            opttryerror <- try(opt <- stats::optim(par = vstart, fn = fnobj,
                                                   obs = data, ddistnam = ddistname, hessian = TRUE,
                                                   method = meth, lower = lower, upper = upper,
                                                   ...), silent = TRUE)
        } else {
            opttryerror <- try(opt <- stats::optim(par = vstart, fn = fnobjcens,
                                                   rcens = rcens, lcens = lcens, 
                                                   icens = icens, ncens = ncens,
                                                   ddistnam = ddistname, 
                                                   pdistnam = pdistname, hessian = TRUE,
                                                   method = meth, 
                                                   lower = lower, upper = upper, ...),
                               silent = TRUE)
        }
        if (inherits(opttryerror, "try-error")) {
            warning("The function optim encountered an error and stopped")
            return(list(estimate = rep(NA, length(vstart)), 
                        convergence = 100,
                        loglik = NA, 
                        hessian = NA))
        }
        if (opt$convergence > 0) {
            warning("The function optim failed to converge, with the error code ",
                    opt$convergence)
            return(list(estimate = rep(NA, length(vstart)), 
                        convergence = opt$convergence,
                        loglik = NA, 
                        hessian = NA))
        }
        return(list(estimate = opt$par, 
                    convergence = opt$convergence,
                    loglik = -opt$value, 
                    hessian = opt$hessian, 
                    optim.function = "optim"))
    } else {
        if (!cens) {
            opttryerror <- try(
                opt <- custom.optim(fn = fnobj, obs = data, ddistnam = ddistname, 
                                    par = vstart, ...), 
                silent = TRUE)
        } else {
            opttryerror <- try(
                opt <- custom.optim(fn = fnobjcens, rcens = rcens, lcens = lcens, 
                                    icens = icens, ncens = ncens,
                                    ddistnam = ddistname, pdistnam = pdistname, 
                                    par = vstart, ...), 
                silent = TRUE)
        }
        if (inherits(opttryerror, "try-error")) {
            print(opttryerror)
            warning("The customized optimization function encountered an error and stopped")
            return(list(estimate = rep(NA, length(vstart)), 
                        convergence = 100,
                        loglik = NA, 
                        hessian = NA))
        }
        if (opt$convergence > 0) {
            warning("The customized optimization function failed to converge, with the error code ",
                    opt$convergence)
            return(list(estimate = rep(NA, length(vstart)), 
                        convergence = opt$convergence,
                        loglik = NA, 
                        hessian = NA))
        }
        return(list(estimate = opt$par, 
                    convergence = opt$convergence,
                    loglik = -opt$value, 
                    hessian = opt$hessian, 
                    optim.function = custom.optim))
    }
}



################################################################################
################################################################################
#' Fits a univariate distribution by maximum likelihood or by matching moments.
#'
#' This function is an alias of the function \code{fitdist} from the package
#' \pkg{fitdistrplus} (Version 0.1-2). The original function was extended to
#' fitting additional distributions. The following continuous distributions can
#' be fitted by this function: normal, lognormal, logistic, exponential, F,
#' Student's t, Beta, Cauchy, Weibull, Gamma.
#' \cr \cr For more details see the assistance page of the function
#' \code{fitdist} from the package \pkg{fitdistrplus}.
#' \cr \cr This function is not intended to be called directly but is internally
#' called in \code{useFitdist}.
#'
#' @name rriskFitdist.cont
#' @aliases rriskFitdist.cont
#' @title Fitting univariate distributions by maximum likelihood or by matching moments
#' @author Matthias Greiner \email{matthias.greiner@@bfr.bund.de} (BfR), \cr
#' Kristin Tolksdorf \email{kristin.tolksdorf@@bfr.bund.de} (BfR), \cr
#' Katharina Schueller \email{schueller@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting), \cr
#' Natalia Belgorodski \email{belgorodski@@stat-up.de} (\acronym{STAT-UP} Statistical Consulting) \cr
#' Marie-Laure Delignette-Muller (coauthor of the package \pkg{fitdistrplus}) \cr
#' Christophe Dutang (coauthor of the package \pkg{fitdistrplus})
#' @usage rriskFitdist.cont(data, distr, method = c("mle", "mme"), start,
#'    chisqbreaks, meancount, ...)
#' @param data A numerical vector, data to be fitted.
#' @param distr A character string \code{name} naming a distribution for which the corresponding
#' density function \code{dname}, the corresponding distribution function
#' \code{pname} and the corresponding quantile function must be defined, or
#' directly the density function.
#' @param method A character string coding for the fitting method: "mle" for the
#' maximum likelihood estimation and "mme" for the matching moment estimation.
#' @param start A named list giving the initial values of parameters of the named
#' distribution. This argument will not be taken into account if \code{method = "mme"},
#' and may be omitted for some distributions for which reasonable starting values
#' are computed if \code{method = "mle"}.
#' @param chisqbreaks A numerical vector defining the breaks of the cell used to
#' compute the chi-square statistic. If omitted, these breaks are automatically
#' computed from the data in order to reach roughly the same number of observations
#' per cell, roughly equal to the argument \code{meancount}, or slightly more if
#' there are some ties.
#' @param meancount The mean number of observations per cell expected for the
#' definition of the breaks of the cells used to compute the chi-squared statistic.
#' @param ... further arguments to be passed to generic function, or to the
#' function \code{rriskMLEdist}, in order to control the optimization method.
#' @return \code{rriskFitdist.cont} returns a list containing 19 items with following fitting results:
#' \item{\code{estimate}}{numeric, a single value or a vector containing estimated parameters.}
#' \item{\code{method}}{the character string representing the used fitting method ("mle" or "mme").}
#' \item{\code{sd}}{the estimated standard errors or \code{NULL} in case of the "mme" method.}
#' \item{\code{cor}}{the estimated correlation matrix or \code{NULL} in case of the "mme" method.}
#' \item{\code{loglik}}{the log-likelihood or \code{NULL} in case of the "mme" method.}
#' \item{\code{aic}}{the Akaike information criterion or \code{NULL} in case of the "mme" method.}
#' \item{\code{bic}}{the BIC or SBC (Schwarz Bayesian criterion) or \code{NULL} in case of the "mme" method.}
#' \item{\code{n}}{the length of the data set.}
#' \item{\code{data}}{the data set.}
#' \item{\code{distname}}{the name of the estimated distribution.}
#' \item{\code{chisq}}{the Chi-squared statistic or \code{NULL}, if not computed.}
#' \item{\code{chisqbreaks}}{breaks used to define cells in the chi-square statistic.}
#' \item{\code{chisqpvalue}}{p-value of the chi-square statistic or \code{NULL}, if not computed.}
#' \item{\code{chisqdf}}{degree of freedom of the chi-square distribution or \code{NULL}, if not computed.}
#' \item{\code{chisqtable}}{a table with observed and theoretical counts used for the Chi-squared calculations.}
#' \item{\code{ad}}{the Anderson-Darling statistic or \code{NULL}, if not computed.}
#' \item{\code{adtest}}{the decision of the Anderson-Darling test or \code{NULL}, if not computed.}
#' \item{\code{ks}}{the Kolmogorov-Smirnov statistic or \code{NULL}, if not computed.}
#' \item{\code{kstest}}{the decision of the Kolmogorov-Smirnov test or \code{NULL}, if not computed.}
#' @keywords fitdistrplus
#' @export
#' @examples
#' set.seed(1)
#' x <- stats::rnorm(5000, mean = 10, sd = 5)
#' rriskFitdist.cont(x, "norm")
#' rriskFitdist.cont(x, "t")
#' 
rriskFitdist.cont <- function(data, distr, 
                              method = c("mle", "mme"), 
                              start, chisqbreaks, 
                              meancount, ...) {
    if (!is.character(distr)) {
        distname <- substring(as.character(match.call()$distr), 2)
    } else distname <- distr
    ddistname <- paste("d", distname, sep = "")
    if (!exists(ddistname, mode = "function")) {
        stop(paste("the ", ddistname, "function must be defined"))
    }
    pdistname <- paste("p", distname, sep = "")
    if (!exists(pdistname, mode = "function")) {
        stop(paste("the ", pdistname, " function must be defined"))
    }
    if (any(method == "mom")) {
        warning("the method \"mom\" for matching moments is deprecated and is replaced by \"mme\".")
    }
    method <- match.arg(method)
    if (!missing(start) & method == "mme") {
        warning("starting values for parameters will not be used with matching moments")
    }
    if (!(is.vector(data) & is.numeric(data) & length(data) > 1)) {
        stop("data must be a numerical vector of length greater than 1")
    }
    n <- length(data)
    if (method == "mme") {
        estimate <- rriskMMEdist(data, distname)
        sd <- NULL
        loglik <- NULL
        aic <- NULL
        bic <- NULL
        correl <- NULL
    } else {
        if (missing(start)) {
            mle <- rriskMLEdist(data, distname, ...)
        } else mle <- rriskMLEdist(data, distname, start, ...)
        if (mle$convergence > 0) {
            stop("the method \"mle\" failed to estimate the parameters, \n                with the error code ",
                 mle$convergence, "\n")
        }
        estimate <- mle$estimate
        if (!is.null(mle$hessian)) {
            if (all(!is.na(mle$hessian))) {
                varcovar <- solve(mle$hessian)
                sd <- sqrt(diag(varcovar))
                correl <- stats::cov2cor(varcovar)
            } else {
                varcovar <- NA
                sd <- NA
                correl <- NA
            }
        } else {
            varcovar <- NA
            sd <- NA
            correl <- NA
        }
        loglik <- mle$loglik
        npar <- length(estimate)
        aic <- -2 * loglik + 2 * npar
        bic <- -2 * loglik + log(n) * npar
    }
    if (is.element(distname, c("binom", "nbinom", "geom", "hyper", "pois"))) {
        discrete <- TRUE
    } else discrete <- FALSE
    if (missing(chisqbreaks)) {
        if (missing(meancount)) {
            meancount <- round(n/((4 * n)^(2/5)))
        }
        sdata <- sort(data)
        if (length(sdata) > ceiling(1.5 * meancount)) {
            limit <- sdata[meancount]
            sdata <- sdata[sdata > limit]
            chisqbreaks <- limit
        } else {
            warning("The sample is too small to automatically define chisqbreaks")
            chisq <- NULL
            chisqbreaks <- NULL
            chisqpvalue <- NULL
            chisqtable <- NULL
            chisqdf <- NULL
        }
        while (length(sdata) > ceiling(1.5 * meancount)) {
            limit <- sdata[meancount]
            sdata <- sdata[sdata > limit]
            chisqbreaks <- c(chisqbreaks, limit)
        }
    }
    if (!is.null(chisqbreaks)) {
        if (!is.numeric(chisqbreaks)) {
            stop("chisqbreaks must be a numerical vector defining the cell boundaries")
        }
        nbreaks <- length(chisqbreaks)
        pbreaks <- do.call(pdistname, c(list(q = chisqbreaks),
                                        as.list(estimate)))
        Fobsbreaks <- stats::ecdf(data)(chisqbreaks)
        Fobsunder <- c(0, Fobsbreaks[1:nbreaks - 1])
        punder <- c(0, pbreaks[1:nbreaks - 1])
        if (pbreaks[nbreaks] == 1 & Fobsbreaks[nbreaks] == 1) {
            p <- pbreaks - punder
            Fobs <- Fobsbreaks - Fobsunder
        } else {
            p <- c(pbreaks - punder, 1 - pbreaks[nbreaks])
            Fobs <- c(Fobsbreaks - Fobsunder, 1 - Fobsbreaks[nbreaks])
        }
        obscounts <- round(Fobs * n)
        theocounts <- p * n
        chisq <- sum(((obscounts - theocounts)^2)/theocounts)
        chisqdf <- length(obscounts) - 1 - length(estimate)
        if (chisqdf > 0) {
            chisqpvalue <- stats::pchisq(chisq, df = chisqdf, lower.tail = FALSE)
        } else chisqpvalue <- NULL
        chisqtable <- as.table(cbind(obscounts, theocounts))
        for (i in 1:length(obscounts) - 1) {
            rownames(chisqtable)[i] <- paste("<=",
                                             signif(chisqbreaks[i], digits = 4))
        }
        rownames(chisqtable)[length(obscounts)] <- paste(">",
                                                         signif(chisqbreaks[i], digits = 4))
    }
    if (!discrete) {
        s <- sort(data)
        obspu <- seq(1, n)/n
        obspl <- seq(0, n - 1)/n
        theop <- do.call(pdistname, c(list(q = s), as.list(estimate)))
        ks <- max(pmax(abs(theop - obspu), abs(theop - obspl)))
        Dmod <- ks * (sqrt(n) + 0.12 + 0.11/sqrt(n))
        if (n >= 30) {
            kstest <- ifelse(Dmod > 1.358, "rejected", "not rejected")
        } else kstest <- NULL
        ad <- -n - sum((2 * (1:n) - 1) * log(theop) + (2 * n + 1 - 2 * (1:n)) * log(1 - theop))/n
        if ((distname == "norm" | distname == "lnorm") & n >= 5) {
            a2mod <- ad * (1 + 0.75/n + 2.25/n^2)
            adtest <- ifelse(a2mod > 0.752, "rejected", "not rejected")
        } else if (distname == "exp" & n >= 5) {
            a2mod <- ad * (1 + 0.6/n)
            adtest <- ifelse(a2mod > 1.321, "rejected", "not rejected")
        } else if (distname == "gamma" & n >= 5) {
            m <- as.list(estimate)$shape
            interp <- stats::approxfun(c(1, 2, 3, 4, 5, 6, 8, 10, 12, 15, 20), 
                                       c(0.786, 0.768, 0.762, 0.759, 0.758, 0.757, 0.755, 0.754, 0.754, 0.754, 0.753), 
                                       yright = 0.752)
            adtest <- ifelse(ad > interp(m), "rejected", "not rejected")
        } else if (distname == "weibull" & n >= 5) {
            a2mod <- ad * (1 + 0.2/sqrt(n))
            adtest <- ifelse(a2mod > 0.757, "rejected", "not rejected")
        } else if (distname == "logis" & n >= 5) {
            a2mod <- ad * (1 + 0.25/n)
            adtest <- ifelse(a2mod > 0.66, "rejected", "not rejected")
        } else if (distname == "cauchy" & n >= 5) {
            interp <- stats::approxfun(c(5, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 100), 
                                       c(1.77, 3.2, 3.77, 4.14, 4.25,
                                         4.05, 3.57, 3.09, 2.48, 2.14, 1.92, 1.52), 
                                       yright = 1.225)
            adtest <- ifelse(ad > interp(n), "rejected", "not rejected")
        } else adtest <- NULL
        if (length(table(data)) != length(data)) {
            warning("Kolmogorov-Smirnov and Anderson-Darling statistics may not be correct with ties")
        }
    } else {
        ks <- NULL
        kstest <- NULL
        ad <- NULL
        adtest <- NULL
    }
    return(structure(list(estimate = estimate, method = method,
                          sd = sd, cor = correl, loglik = loglik, aic = aic, bic = bic,
                          n = n, data = data, distname = distname, chisq = chisq,
                          chisqbreaks = chisqbreaks, chisqpvalue = chisqpvalue,
                          chisqdf = chisqdf, chisqtable = chisqtable, ad = ad,
                          adtest = adtest, ks = ks, kstest = kstest), 
                     class = "fitdist"))
} # end of rriskFitdist.cont()


