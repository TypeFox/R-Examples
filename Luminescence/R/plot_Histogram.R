#' Plot a histogram with separate error plot
#'
#' Function plots a predefined histogram with an accompanying error plot as
#' suggested by Rex Galbraith at the UK LED in Oxford 2010.
#'
#' If the normal curve is added, the y-axis in the histogram will show the
#' probability density.\cr\cr 
#' A statistic summary, i.e. a collection of statistic measures of 
#' centrality and dispersion (and further measures) can be added by specifying 
#' one or more of the following keywords: \code{"n"} (number of samples),
#' \code{"mean"} (mean De value), \code{"mean.weighted"} (error-weighted mean),
#' \code{"median"} (median of the De values), \code{"sdrel"} (relative standard
#' deviation in percent), \code{"sdrel.weighted"} (error-weighted relative 
#' standard deviation in percent), \code{"sdabs"} (absolute standard deviation),
#' \code{"sdabs.weighted"} (error-weighted absolute standard deviation), 
#' \code{"serel"} (relative standard error), \code{"serel.weighted"} (
#' error-weighted relative standard error), \code{"seabs"} (absolute standard
#' error), \code{"seabs.weighted"} (error-weighted absolute standard error), 
#' \code{"kurtosis"} (kurtosis) and \code{"skewness"} (skewness).
#'
#' @param data \code{\link{data.frame}} or \code{\linkS4class{RLum.Results}}
#' object (required): for \code{data.frame}: two columns: De (\code{data[,1]})
#' and De error (\code{data[,2]}) 
#' @param na.rm \code{\link{logical}} (with default): excludes \code{NA}
#' values from the data set prior to any further operations.
#' @param mtext \code{\link{character}} (optional): further sample information
#' (\link{mtext}).
#' @param cex.global \code{\link{numeric}} (with default): global scaling
#' factor.
#' @param se \code{\link{logical}} (optional): plots standard error points over
#' the histogram, default is \code{FALSE}.
#' @param rug \code{\link{logical}} (optional): adds rugs to the histogram,
#' default is \code{TRUE}.
#' @param normal_curve \code{\link{logical}} (with default): adds a normal
#' curve to the histogram. Mean and sd are calculated from the input data. More
#' see details section.
#' @param summary \code{\link{character}} (optional): add statistic measures of 
#' centrality and dispersion to the plot. Can be one or more of several 
#' keywords. See details for available keywords.
#' @param summary.pos \code{\link{numeric}} or \code{\link{character}} (with
#' default): optional position coordinates or keyword (e.g. \code{"topright"})
#' for the statistical summary. Alternatively, the keyword \code{"sub"} may be
#' specified to place the summary below the plot header. However, this latter
#' option in only possible if \code{mtext} is not used. In case of coordinate
#' specification, y-coordinate refers to the right y-axis.
#' @param colour \code{\link{numeric}} or \link{character} (with default):
#' optional vector of length 4 which specifies the colours of the following
#' plot items in exactly this order: histogram bars, rug lines, normal
#' distribution curve and standard error points\cr (e.g., \code{c("grey",
#' "black", "red", "grey")}).
#' @param \dots further arguments and graphical parameters passed to
#' \code{\link{plot}} or \code{\link{hist}}. If y-axis labels are provided,
#' these must be specified as a vector of length 2 since the plot features two
#' axes (e.g. \code{ylab = c("axis label 1", "axis label 2")}). Y-axes limits
#' (\code{ylim}) must be provided as vector of length four, with the first two
#' elements specifying the left axes limits and the latter two elements giving
#' the right axis limits.
#' @note The input data is not restricted to a special type.
#' @section Function version: 0.4.4
#' @author Michael Dietze, GFZ Potsdam (Germany), \cr Sebastian Kreutzer,
#' IRAMAT-CRP2A, Universite Bordeaux Montaigne (France)
#' @seealso \code{\link{hist}}, \code{\link{plot}}
#' @examples
#'
#' ## load data
#' data(ExampleData.DeValues, envir = environment())
#' ExampleData.DeValues <-
#'   Second2Gray(ExampleData.DeValues$BT998, dose.rate = c(0.0438,0.0019))
#'
#' ## plot histogram the easiest way
#' plot_Histogram(ExampleData.DeValues)
#'
#' ## plot histogram with some more modifications
#' plot_Histogram(ExampleData.DeValues,
#'                rug = TRUE,
#'                normal_curve = TRUE,
#'                cex.global = 0.9,
#'                pch = 2,
#'                colour = c("grey", "black", "blue", "green"),
#'                summary = c("n", "mean", "sdrel"),
#'                summary.pos = "topleft",
#'                main = "Histogram of De-values",
#'                mtext = "Example data set",
#'                ylab = c(expression(paste(D[e], " distribution")),
#'                         "Standard error"),
#'                xlim = c(100, 250),
#'                ylim = c(0, 0.1, 5, 20))
#'
#'
#' @export
plot_Histogram <- function(
  data,
  na.rm = TRUE,
  mtext,
  cex.global,
  se,
  rug,
  normal_curve,
  summary,
  summary.pos,
  colour,
  ...
) {
  
  # Integrity tests ---------------------------------------------------------
  ## check/adjust input data structure
  if(is(data, "RLum.Results") == FALSE &
     is(data, "data.frame") == FALSE) {
    
    stop(paste("[plot_Histogram()] Input data format is neither",
               "'data.frame' nor 'RLum.Results'"))
  } else {
    
    if(is(data, "RLum.Results") == TRUE) {
      data <- get_RLum(data, "data")
    }
  }
  
  ## handle error-free data sets
  if(length(data) < 2) {
    data <- cbind(data, rep(NA, length(data)))
  }
  
  
  ## Set general parameters ---------------------------------------------------
  ## Check/set default parameters
  if(missing(cex.global) == TRUE) {
    cex.global <- 1
  }
  
  if(missing(mtext) == TRUE) {
    mtext <- ""
  }
  
  if(missing(se) == TRUE) {
    se = TRUE
  }
  
  if(missing(rug) == TRUE) {
    rug = TRUE
  }
  
  if(missing(colour) == TRUE) {
    colour = c("white", "black", "red", "black")
  }
  
  if(missing(summary) == TRUE) {
    summary <- ""
  }
  
  if(missing(summary.pos) == TRUE) {
    summary.pos <- "sub"
  }
  
  if(missing(normal_curve) == TRUE) {
    normal_curve = FALSE
  }
  
  ## read out additional arguments list
  extraArgs <- list(...)
  
  ## define fun
  if("fun" %in% names(extraArgs)) {
    fun <- extraArgs$fun
  } else {
    fun <- FALSE
  }
  
  ## optionally, count and exclude NA values and print result
  if(na.rm == TRUE) {
    n.NA <- sum(is.na(data[,1]))
    if(n.NA == 1) {
      print("1 NA value excluded.")
    } else if(n.NA > 1) {
      print(paste(n.NA, "NA values excluded."))
    }
    data <- data[!is.na(data[,1]),]
  }
  
  if("main" %in% names(extraArgs)) {
    main.plot <- extraArgs$main
  } else {
    main.plot <- "Histogram"
  }
  
  if("xlab" %in% names(extraArgs)) {
    xlab.plot <- extraArgs$xlab
  } else {
    xlab.plot <- expression(paste(D[e], " (Gy)"))
  }
  
  if("ylab" %in% names(extraArgs)) {
    ylab.plot <- extraArgs$ylab
  } else {
    ylab.plot <- c("Frequency",
                   "Standard error")
  }
  
  if("breaks" %in% names(extraArgs)) {
    breaks.plot <- extraArgs$breaks
  } else {
    breaks.plot <- hist(x = data[,1],
                        plot = FALSE)$breaks
  }
  
  if("xlim" %in% names(extraArgs)) {
    xlim.plot <- extraArgs$xlim
  } else {
    xlim.plot <- range(breaks.plot)
  }
  
  if("ylim" %in% names(extraArgs)) {
    ylim.plot <- extraArgs$ylim
  } else {
    H.lim <- hist(data[,1],
                  breaks = breaks.plot,
                  plot = FALSE)
    if(normal_curve == TRUE) {
      left.ylim <- c(0, max(H.lim$density))
    } else {
      left.ylim <- c(0, max(H.lim$counts))
    }
    range.error <- try(expr = range(data[,2], na.rm = TRUE),
                       silent = TRUE)
    range.error[1] <- ifelse(is.infinite(range.error[1]), 0, range.error[1])
    range.error[2] <- ifelse(is.infinite(range.error[2]), 0, range.error[2])
    ylim.plot <- c(left.ylim, range.error)
  }
  
  if("pch" %in% names(extraArgs)) {
    pch.plot <- extraArgs$pch
  } else {
    pch.plot <- 1
  }
  ## Set plot area format
  par(mar = c(4.5, 4.5, 4.5, 4.5),
      cex = cex.global)
  
  ## Plot histogram -----------------------------------------------------------
  HIST <- hist(data[,1],
               main = "",
               xlab = xlab.plot,
               ylab = ylab.plot[1],
               xlim = xlim.plot,
               ylim = ylim.plot[1:2],
               breaks = breaks.plot,
               freq = !normal_curve,
               col = colour[1]
  )
  
  ## add title
  title(line = 2,
        main = main.plot)
  
  ## Optionally, add rug ------------------------------------------------------
  if(rug == TRUE) {rug(data[,1], col = colour[2])}
  
  ## Optionally, add a normal curve based on the data -------------------------
  if(normal_curve == TRUE){
    ## cheat the R check routine, tztztz how neat
    x <- NULL
    rm(x)
    
    ## add normal distribution curve
    curve(dnorm(x,
                mean = mean(na.exclude(data[,1])),
                sd = sd(na.exclude(data[,1]))),
          col = colour[3],
          add = TRUE,
          lwd = 1.2 * cex.global)
  }
  
  ## calculate and paste statistical summary
  data.stats <- list(data = data)
  
  ## calculate and paste statistical summary
  De.stats <- matrix(nrow = length(data), ncol = 18)
  colnames(De.stats) <- c("n",
                          "mean",
                          "mean.weighted",
                          "median",
                          "median.weighted",
                          "kde.max",
                          "sd.abs",
                          "sd.rel",
                          "se.abs",
                          "se.rel",
                          "q25",
                          "q75",
                          "skewness",
                          "kurtosis",
                          "sd.abs.weighted",
                          "sd.rel.weighted",
                          "se.abs.weighted",
                          "se.rel.weighted")
  
  for(i in 1:length(data)) {
    statistics <- calc_Statistics(data)
    De.stats[i,1] <- statistics$weighted$n
    De.stats[i,2] <- statistics$unweighted$mean
    De.stats[i,3] <- statistics$weighted$mean
    De.stats[i,4] <- statistics$unweighted$median
    De.stats[i,5] <- statistics$unweighted$median
    De.stats[i,7] <- statistics$unweighted$sd.abs
    De.stats[i,8] <- statistics$unweighted$sd.rel
    De.stats[i,9] <- statistics$unweighted$se.abs
    De.stats[i,10] <- statistics$weighted$se.rel
    De.stats[i,11] <- quantile(data[,1], 0.25)
    De.stats[i,12] <- quantile(data[,1], 0.75)
    De.stats[i,13] <- statistics$unweighted$skewness
    De.stats[i,14] <- statistics$unweighted$kurtosis
    De.stats[i,15] <- statistics$weighted$sd.abs
    De.stats[i,16] <- statistics$weighted$sd.rel
    De.stats[i,17] <- statistics$weighted$se.abs
    De.stats[i,18] <- statistics$weighted$se.rel
    
    ##kdemax - here a little doubled as it appears below again
    De.density <-density(x = data[,1],
                         kernel = "gaussian",
                         from = xlim.plot[1],
                         to = xlim.plot[2])
    
    De.stats[i,6] <- De.density$x[which.max(De.density$y)]
  }
  
  label.text = list(NA)
  
  if(summary.pos[1] != "sub") {
    n.rows <- length(summary)
    
    for(i in 1:length(data)) {
      stops <- paste(rep("\n", (i - 1) * n.rows), collapse = "")
      
      summary.text <- character(0)
      
      for(j in 1:length(summary)) {
        summary.text <- c(summary.text,
                          paste(
                            "",
                            ifelse("n" %in% summary[j] == TRUE,
                                   paste("n = ",
                                         De.stats[i,1],
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("mean" %in% summary[j] == TRUE,
                                   paste("mean = ",
                                         round(De.stats[i,2], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("mean.weighted" %in% summary[j] == TRUE,
                                   paste("weighted mean = ",
                                         round(De.stats[i,3], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("median" %in% summary[j] == TRUE,
                                   paste("median = ",
                                         round(De.stats[i,4], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("median.weighted" %in% summary[j] == TRUE,
                                   paste("weighted median = ",
                                         round(De.stats[i,5], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("kdemax" %in% summary[j] == TRUE,
                                   paste("kdemax = ",
                                         round(De.stats[i,6], 2),
                                         " \n ",
                                         sep = ""),
                                   ""),
                            ifelse("sdabs" %in% summary[j] == TRUE,
                                   paste("sd = ",
                                         round(De.stats[i,7], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("sdrel" %in% summary[j] == TRUE,
                                   paste("rel. sd = ",
                                         round(De.stats[i,8], 2), " %",
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("seabs" %in% summary[j] == TRUE,
                                   paste("se = ",
                                         round(De.stats[i,9], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("serel" %in% summary[j] == TRUE,
                                   paste("rel. se = ",
                                         round(De.stats[i,10], 2), " %",
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("skewness" %in% summary[j] == TRUE,
                                   paste("skewness = ",
                                         round(De.stats[i,13], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("kurtosis" %in% summary[j] == TRUE,
                                   paste("kurtosis = ",
                                         round(De.stats[i,14], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("sdabs.weighted" %in% summary[j] == TRUE,
                                   paste("abs. weighted sd = ",
                                         round(De.stats[i,15], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("sdrel.weighted" %in% summary[j] == TRUE,
                                   paste("rel. weighted sd = ",
                                         round(De.stats[i,16], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("seabs.weighted" %in% summary[j] == TRUE,
                                   paste("abs. weighted se = ",
                                         round(De.stats[i,17], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            ifelse("serel.weighted" %in% summary[j] == TRUE,
                                   paste("rel. weighted se = ",
                                         round(De.stats[i,18], 2),
                                         "\n",
                                         sep = ""),
                                   ""),
                            sep = ""))
      }
      
      summary.text <- paste(summary.text, collapse = "")
      
      label.text[[length(label.text) + 1]] <- paste(stops,
                                                    summary.text,
                                                    stops,
                                                    sep = "")
    }
  } else {
    for(i in 1:length(data)) {
      
      summary.text <- character(0)
      
      for(j in 1:length(summary)) {
        summary.text <- c(summary.text,
                          ifelse("n" %in% summary[j] == TRUE,
                                 paste("n = ",
                                       De.stats[i,1],
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("mean" %in% summary[j] == TRUE,
                                 paste("mean = ",
                                       round(De.stats[i,2], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("mean.weighted" %in% summary[j] == TRUE,
                                 paste("weighted mean = ",
                                       round(De.stats[i,3], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("median" %in% summary[j] == TRUE,
                                 paste("median = ",
                                       round(De.stats[i,4], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("median.weighted" %in% summary[j] == TRUE,
                                 paste("weighted median = ",
                                       round(De.stats[i,5], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("kdemax" %in% summary[j] == TRUE,
                                 paste("kdemax = ",
                                       round(De.stats[i,6], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdrel" %in% summary[j] == TRUE,
                                 paste("rel. sd = ",
                                       round(De.stats[i,8], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdabs" %in% summary[j] == TRUE,
                                 paste("abs. sd = ",
                                       round(De.stats[i,7], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("serel" %in% summary[j] == TRUE,
                                 paste("rel. se = ",
                                       round(De.stats[i,10], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("seabs" %in% summary[j] == TRUE,
                                 paste("abs. se = ",
                                       round(De.stats[i,9], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("skewness" %in% summary[j] == TRUE,
                                 paste("skewness = ",
                                       round(De.stats[i,13], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("kurtosis" %in% summary[j] == TRUE,
                                 paste("kurtosis = ",
                                       round(De.stats[i,14], 2),
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdabs.weighted" %in% summary[j] == TRUE,
                                 paste("abs. weighted sd = ",
                                       round(De.stats[i,15], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("sdrel.weighted" %in% summary[j] == TRUE,
                                 paste("rel. weighted sd = ",
                                       round(De.stats[i,16], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("seabs.weighted" %in% summary[j] == TRUE,
                                 paste("abs. weighted se = ",
                                       round(De.stats[i,17], 2), " %",
                                       " | ",
                                       sep = ""),
                                 ""),
                          ifelse("serel.weighted" %in% summary[j] == TRUE,
                                 paste("rel. weighted se = ",
                                       round(De.stats[i,18], 2), " %",
                                       " | ",
                                       sep = ""),
                                 "")
        )
      }
      
      summary.text <- paste(summary.text, collapse = "")
      
      label.text[[length(label.text) + 1]]  <- paste(
        "  ",
        summary.text,
        sep = "")
    }
    
    ## remove outer vertical lines from string
    for(i in 2:length(label.text)) {
      label.text[[i]] <- substr(x = label.text[[i]],
                                start = 3,
                                stop = nchar(label.text[[i]]) - 3)
    }
  }
  
  ## remove dummy list element
  label.text[[1]] <- NULL
  
  ## convert keywords into summary placement coordinates
  if(missing(summary.pos) == TRUE) {
    summary.pos <- c(xlim.plot[1], ylim.plot[2])
    summary.adj <- c(0, 1)
  } else if(length(summary.pos) == 2) {
    summary.pos <- summary.pos
    summary.adj <- c(0, 1)
  } else if(summary.pos[1] == "topleft") {
    summary.pos <- c(xlim.plot[1], ylim.plot[2])
    summary.adj <- c(0, 1)
  } else if(summary.pos[1] == "top") {
    summary.pos <- c(mean(xlim.plot), ylim.plot[2])
    summary.adj <- c(0.5, 1)
  } else if(summary.pos[1] == "topright") {
    summary.pos <- c(xlim.plot[2], ylim.plot[2])
    summary.adj <- c(1, 1)
  }  else if(summary.pos[1] == "left") {
    summary.pos <- c(xlim.plot[1], mean(ylim.plot[1:2]))
    summary.adj <- c(0, 0.5)
  } else if(summary.pos[1] == "center") {
    summary.pos <- c(mean(xlim.plot), mean(ylim.plot[1:2]))
    summary.adj <- c(0.5, 0.5)
  } else if(summary.pos[1] == "right") {
    summary.pos <- c(xlim.plot[2], mean(ylim.plot[1:2]))
    summary.adj <- c(1, 0.5)
  }else if(summary.pos[1] == "bottomleft") {
    summary.pos <- c(xlim.plot[1], ylim.plot[1])
    summary.adj <- c(0, 0)
  } else if(summary.pos[1] == "bottom") {
    summary.pos <- c(mean(xlim.plot), ylim.plot[1])
    summary.adj <- c(0.5, 0)
  } else if(summary.pos[1] == "bottomright") {
    summary.pos <- c(xlim.plot[2], ylim.plot[1])
    summary.adj <- c(1, 0)
  }
  
  ## add summary content
  for(i in 1:length(data.stats)) {
    if(summary.pos[1] != "sub") {
      text(x = summary.pos[1],
           y = summary.pos[2],
           adj = summary.adj,
           labels = label.text[[i]],
           col = colour[2],
           cex = cex.global * 0.8)
    } else {
      if(mtext == "") {
        mtext(side = 3,
              line = 1 - i,
              text = label.text[[i]],
              col = colour[2],
              cex = cex.global * 0.8)
      }
    }
  }
  
  ## Optionally, add standard error plot --------------------------------------
  if(sum(is.na(data[,2])) == length(data[,2])) {
    se <- FALSE
  }
  
  if(se == TRUE) {
    par(new = TRUE)
    plot.data <- data[!is.na(data[,2]),]
    
    plot(x = plot.data[,1],
         y = plot.data[,2],
         xlim = xlim.plot,
         ylim = ylim.plot[3:4],
         pch = pch.plot,
         col = colour[4],
         main = "",
         xlab = "",
         ylab = "",
         axes = FALSE,
         frame.plot = FALSE
    )
    axis(side = 4,
         labels = TRUE,
         cex = cex.global
    )
    mtext(ylab.plot[2],
          side = 4,
          line = 3,
          cex = cex.global)
    
    #    par(new = FALSE)
  }
  
  ## Optionally add user-defined mtext
  mtext(side = 3,
        line = 0.5,
        text = mtext,
        cex = 0.8 * cex.global)
  
  ## FUN by R Luminescence Team
  if(fun==TRUE){sTeve()}
  
}
