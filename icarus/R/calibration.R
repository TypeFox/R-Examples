# copyright (C) 2014-2016 A.Rebecq
# This function executes easy calibration with just data and matrix of margins

#########
#' Calibration on margins
#' @description
#' Performs calibration on margins with several methods and customizable parameters
#' @param data The dataframe containing the survey data
#' @param marginMatrix The matrix giving the margins for each column variable included
#' in the calibration problem
#' @param colWeights The name of the column containing the initial weights in the survey
#' dataframe
#' @param method The method used to calibrate. Can be "linear", "raking", "logit", "truncated"
#' @param bounds Two-element vector containing the lower and upper bounds for bounded methods
#' ("truncated" and "logit")
#' @param costs The penalized calibration method will be used, using costs defined by this
#' vector. Must match the number of rows of marginMatrix. Negative of non-finite costs are given
#' an infinite cost (coefficient of C^-1 matrix is 0)
#' @param gap Only useful for penalized calibration. Sets the maximum gap between max and min
#' calibrated weights / initial weights ratio (and thus is similar to the "bounds"
#' parameter used in regular calibration)
#' @param popTotal Precise the total population if margins are defined by relative value in
#' marginMatrix (percentages)
#' @param pct If TRUE, margins for categorical variables are considered to
#' be entered as percentages. popTotal must then be set. (FALSE by default)
#' @param scale If TRUE, stats (including bounds) on ratio calibrated weights / initial weights are
#' done on a vector multiplied by the weighted non-response ratio (ratio population total /
#' total of initial weights). Has same behavior as "ECHELLE=0" in Calmar.
#' @param description If TRUE, output stats about the calibration process as well as the
#' graph of the density of the ratio calibrated weights / initial weights
#' @param maxIter The maximum number of iterations before stopping
#' @param check performs a few check about the dataframe. TRUE by default
#' @param uCostPenalized Unary cost by which every cost is "costs" column is multiplied
#' @param lambda The initial ridge lambda used in penalized calibration. By default, the initial
#' lambda is automatically chosen by the algorithm, but you can speed up the search for the optimum
#' if you already know a lambda close to the lambda_opt corresponding to the gap you set. Be careful,
#' the search zone is reduced when a lambda is set by the user, so the program may not converge
#' if the lambda set is too far from the lambda_opt.
#' @param precisionBounds Only used for calibration on minimum bounds. Desired precision
#' for lower and upper reweighting factor, both bounds being as close to 1 as possible
#' @param forceSimplex Only used for calibration on tight bounds.Bisection algorithm is used
#' for matrices whose size exceed 1e8. forceSimplex = TRUE forces the use of the simplex algorithm
#' whatever the size of the problem (you might want to set this parameter to TRUE if you
#' have a large memory size)
#' @param forceBisection Only used for calibration on tight bounds. Forces the use of the bisection
#' algorithm to solve calibration on tight bounds
#' @param colCalibratedWeights The name of the column of final calibrated weights
#' @param exportDistributionImage File name to which the density plot shown when
#' description is TRUE is exported. Requires package "ggplot2"
#' @param exportDistributionTable File name to which the distribution table of before/after
#' weights shown when description is TRUE is exported. Requires package "xtable"
#' 
#' @examples
#' N <- 300 ## population total
#' ## Horvitz Thompson estimator of the mean: 1.666667
#' weightedMean(data_ex2$cinema, data_ex2$poids, N) 
#' ## Enter calibration margins:
#' mar1 <- c("categ",3,80,90,60)
#' mar2 <- c("sexe",2,140,90,0)
#' mar3 <- c("service",2,100,130,0)
#' mar4 <- c("salaire", 0, 470000,0,0)
#' margins <- rbind(mar1, mar2, mar3, mar4)
#' ## Compute calibrated weights with raking ratio method
#' wCal <- calibration(data=data_ex2, marginMatrix=margins, colWeights="poids"
#'                             , method="raking", description=FALSE)
#' ## Calibrated estimate: 2.471917
#' weightedMean(data_ex2$cinema, wCal, N)
#'
#' @references Deville, Jean-Claude, and Carl-Erik Sarndal. "Calibration estimators in survey sampling." 
#' Journal of the American statistical Association 87.418 (1992): 376-382.
#' @references Bocci, J., and C. Beaumont. "Another look at ridge calibration." 
#' Metron 66.1 (2008): 5-20.
#' @references Vanderhoeft, Camille. Generalised calibration at statistics Belgium: SPSS Module G-CALIB-S and current practices. 
#' Inst. National de Statistique, 2001.
#' @references Le Guennec, Josiane, and Olivier Sautory. "Calmar 2: Une nouvelle version 
#' de la macro calmar de redressement d'echantillon par calage." Journees de Methodologie Statistique, 
#' Paris. INSEE (2002).
#' @return column containing the final calibrated weights
#'
#' @export
calibration = function(data, marginMatrix, colWeights, method="linear", bounds=NULL
                       , costs=NULL, gap=NULL, popTotal=NULL, pct=FALSE, scale=NULL, description=TRUE
                       , maxIter=2500, check=TRUE, uCostPenalized=1, lambda=NULL, precisionBounds=1e-4, forceSimplex=FALSE, forceBisection=FALSE
                       , colCalibratedWeights="calWeights", exportDistributionImage=NULL, exportDistributionTable=NULL) {
  
  # By default, scale is TRUE when popTotal is not NULL, false otherwise
  if(is.null(popTotal)) {
    scale <- FALSE
  } else {
    scale <- TRUE
  }
  
  ## Clean costs for penalized calibration
  #   if(!is.null(costs)) {
  #     costs <- cleanCosts(costs, uCostPenalized)
  #   }
  
  if(check) {
    
    # Check NAs on calibration variables
    matrixTestNA = missingValuesMargins(data, marginMatrix)
    testNA = as.numeric(matrixTestNA[,2])
    if(sum(testNA) > 0) {
      print(matrixTestNA)
      stop("NAs found in calibration variables")
    }
    
    # check if number of modalities in calibration variables matches marginMatrix
    if(!checkNumberMargins(data, marginMatrix)) stop("Error in number of modalities.")
    
  }
  
  marginCreation <- createFormattedMargins(data, marginMatrix, popTotal, pct)
  matrixCal = marginCreation[[2]]
  formattedMargins = marginCreation[[1]]
  
  # Same rule as in "Calmar" for SAS : if scale is TRUE,
  # calibration is done on weights adjusted for nonresponse
  # (uniform adjustment)
  weights <- as.numeric(data.matrix(data[colWeights]))
  
  ## For tests
  #   Xs_glob <<- matrixCal
  #   total_glob <<- formattedMargins
  #   d_glob <<- weights
  
  if(scale) {
    if(is.null(popTotal)) {
      stop("When scale is TRUE, popTotal cannot be NULL")
    }
    weights <- weights*(popTotal / sum(data.matrix(data[colWeights])) )
    
  }
  
  # TODO : add catching of "no convergence" exception
  if(is.null(costs)) {
    
    g <- NULL
    
    if( (is.numeric(bounds)) || (method != "min") ) {
      g <- calib(Xs=matrixCal, d=weights, total=formattedMargins, method=method, bounds=bounds, maxIter=maxIter)
    } else {
      if( (bounds == "min") || (method == "min")) {
        g <- minBoundsCalib(Xs=matrixCal, d=weights, total=formattedMargins
                            , q=rep(1,nrow(matrixCal)), maxIter=maxIter, description=description, precisionBounds=precisionBounds, forceSimplex=forceSimplex, forceBisection=forceBisection)
      }
    }
    
    data[colCalibratedWeights] = g*weights
    
  } else {
    
    # Forbid popTotal null when gap is selected
    if(!is.null(gap) && is.null(popTotal)) {
      warning("popTotal NULL when gap is selected is a risky setting !")
    }
    
    ## Format costs
    costsFormatted <- formatCosts(costs, marginMatrix, popTotal)
    
    wCal = penalizedCalib(Xs=matrixCal, d=weights, total=formattedMargins, method=method
                          , bounds=bounds, costs=costsFormatted, uCostPenalized=uCostPenalized
                          , maxIter=maxIter, lambda=lambda, gap=gap)
    data[colCalibratedWeights] = data.matrix(wCal)
    g = wCal / weights
  }
  
  if(description) {
    writeLines("")
    writeLines("################### Summary of before/after weight ratios ###################")
  }
  # popTotalComp is popTotal computed from sum of calibrated weights
  popTotalComp <- sum(data[colCalibratedWeights])
  
  weightsRatio = g
  
  if(description) {
    
    writeLines(paste("Calibration method : ",method, sep=""))
    
    if(! (method %in% c("linear","raking")) && ! is.null(bounds) ) {
      
      if(is.numeric(bounds)) {
        writeLines(paste("\t L bound : ",bounds[1], sep=""))
        writeLines(paste("\t U bound : ",bounds[2], sep=""))
      }
      
      if( (bounds == "min") || (method == "min") ) {
        writeLines(paste("\t L bound : ",round(min(g),4), sep=""))
        writeLines(paste("\t U bound : ",round(max(g),4), sep=""))
      }
      
    }
    
    writeLines(paste("Mean : ",round(mean(weightsRatio),4), sep=""))
    quantileRatio <- round(stats::quantile(weightsRatio, probs=c(0,0.01,0.1,0.25,0.5,0.75,0.9,0.99,1)),4)
    print(quantileRatio)
    
  }
  
  ## Export in TeX
  if(!is.null(exportDistributionTable)) {
    
    if (!requireNamespace("xtable", quietly = TRUE)) {
      stop("Package xtable needed for exportDistributionTable to work. Please install it.",
           call. = FALSE)
    }
    
    # Linear or raking ratio
    if(is.null(bounds)) {
      
      newNames <- names(quantileRatio)
      newNames <- c(newNames,"Mean")
      
      statsRatio <- c(quantileRatio,mean(quantileRatio))
      names(statsRatio) <- newNames
      
    } else {
      newNames <- names(quantileRatio)
      newNames <- c("L",newNames,"U","Mean")
      
      statsRatio <- c(bounds[1],quantileRatio,bounds[2],mean(quantileRatio))
      names(statsRatio) <- newNames
    }
    
    latexQuantiles <- xtable::xtable(as.data.frame(t(statsRatio)))
    
    # Notice that there is one extra column in align(latexQuantiles)
    # since we haven't specified yet to exclide rownames
    if(is.null(bounds)) {
      xtable::align(latexQuantiles) <- "|c|ccccccccc||c|"
    } else {
      xtable::align(latexQuantiles) <- "|c|c|ccccccccc|c||c|"
    }
    
    
    print(latexQuantiles,  include.rownames = FALSE, include.colnames = TRUE,
           floating = FALSE, file=exportDistributionTable)
  }
  
  
  if(description) {
    writeLines("")
    writeLines("################### Comparison Margins Before/After calibration ###################")
    print(calibrationMarginStats(data=data, marginMatrix=marginMatrix, popTotal=popTotal, pct=pct, colWeights=colWeights, colCalibratedWeights=colCalibratedWeights))
  }
  
  # Plot density of weights ratio
  if(description) {
    
    if(requireNamespace("ggplot2")) {
      
      densityPlot = ggplot2::ggplot(data.frame(weightsRatio), ggplot2::aes(x=weightsRatio)) + ggplot2::geom_density(alpha=0.5, fill="#FF6666", size=1.25, adjust=2) + ggplot2::theme_bw()
      print(densityPlot)
      
      if(!is.null(exportDistributionImage)) {
        ggplot2::ggsave(densityPlot, file=exportDistributionImage)
      }
      
    } else {
      warning("Require package ggplot2 to plot weights ratio")
    }
  }
  
  return(g*weights)
}
