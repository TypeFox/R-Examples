#' create \code{\link{sdcProblem-class}}-objects
#'
#' Function \code{\link{makeProblem}} is used to create
#' \code{\link{sdcProblem-class}}-objects.
#'
#' @param data a data frame featuring at least one column for each desired
#' dimensional variable. Optionally the input data can feature variables
#' that contain information on cell counts, weights that should be used during
#' the cut and branch algorithm, additional numeric variables or variables that
#' hold information on sampling weights.
#' @param dimList a named list with each list element being either a data-frame or a link to a .csv-file containing the complete level-hierarchy of a dimensional variable using a top-to-bottom approach. The list names correspond to variable names that must exist in argument \code{data}. The level-hierarchy must be specified as follows:
#' \itemize{
#' \item list-element is a data-frame that must contain exactly 2 columns with the first column specifying levels and the second column holding variable-codes.
#' \itemize{
#' \item first column: a character vector specifying levels with each vector element being a string only containing of '@@'s from length 1 to n. If a vector element consists of \code{i}-chars, the corresponding code is of level \code{i}. The code '@@' (one character) equals the grand total (level=1).
#' \item second column: a character vector specifying level codes
#' }
#' \item list-element is full path to a .csv-file with two columns seperated by semicolons (;) having the same structure as the data.frame described above
#' }
#' @param dimVarInd numeric vector (or NULL) defining the column-indices of dimensional variables (defining the table) within argument \code{data}
#' @param freqVarInd numeric vector (or NULL) defining the column-indices of a variable holding counts within argument \code{data}
#' @param numVarInd numeric vector (or NULL) defining the column-indices of additional numeric variables available in argument \code{data}
#' @param weightInd numeric vector of length 1 (or NULL) defining the column-index of a variable holding weights that should be used during as objective coefficients during the cut and branch algorithm to protect primary sensitive cells within argument \code{data}
#' @param sampWeightInd numeric vector of length 1 (or NULL) defining the column-index of a variable holding sampling weights within argument \code{data}
#'
#' @return a \code{\link{sdcProblem-class}}-object
#'
#' @examples
#' # loading micro data
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/microData1.RData", sep="")
#' microData <- get(load(fn))
#'
#' # having a look at the data structure
#' str(microData)
#'
#' # we can observe that we have a micro data set consisting of two spanning
#' # variables ('region' and 'gender') and one numeric variable ('val')
#'
#' # specify structure of hierarchical variable 'region'
#' # levels 'A' to 'D' sum up to a Total
#' dim.region <- data.frame(
#'  levels=c('@@','@@@@','@@@@','@@@@','@@@@'),
#'  codes=c('Total', 'A','B','C','D'),
#'  stringsAsFactors=FALSE)
#'
#' # specify structure of hierarchical variable 'gender'
#' # levels 'male' and 'female' sum up to a Total
#' dim.gender <- data.frame(
#'  levels=c('@@','@@@@','@@@@'),
#'  codes=c('Total', 'male','female'),
#'  stringsAsFactors=FALSE)
#'
#' # create a list with each element being a data-frame containing information
#' # on a dimensional variables
#' dimList <- list(dim.region, dim.gender)
#'
#' # name the list:
#' # - first list-element: corresponds to variable 'region'
#' # - second list-element: corresponds to variable 'gender'
#' names(dimList) <- c('region', 'gender')
#'
#' # specify the indices where dimensional variables are located
#' # within the input data
#'
#' # - variable 'region': first column
#' # - variable 'gender': second column
#' dimVarInd <- c(1,2)
#'
#' # third column containts a numeric variable
#' numVarInd <- 3
#'
#' # no variables holding counts, numeric values, weights or sampling
#' # weights are available in the input data
#' freqVarInd <- weightInd <- sampWeightInd <- NULL
#'
#' # creating an object of class \code{\link{sdcProblem-class}}
#' problem <- makeProblem(
#'  data=microData,
#'  dimList=dimList,
#'  dimVarInd=dimVarInd,
#'  freqVarInd=freqVarInd,
#'  numVarInd=numVarInd,
#'  weightInd=weightInd,
#'  sampWeightInd=sampWeightInd)
#'
#' # what do we have?
#' print(class(problem))
#' @rdname makeProblem
#' @export makeProblem
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
makeProblem <- function(data, dimList, dimVarInd, freqVarInd=NULL, numVarInd=NULL, weightInd=NULL,sampWeightInd=NULL) {
  # returns an object of class 'sdcProblem'
  # 'doPrep()' is the old function 'newDimInfo()'
  # since it also recodes inputData eventually, it was renamed
  doPrep <- function(inputData, inputDims) {
    if ( any(sapply(inputDims, class) != "dimVar") ) {
      stop("Error: all elements of 'inputDims' must be of class 'dimVar'!\n")
    }
    if ( class(inputData) != "dataObj") {
      stop("Error: 'inputData' be of class 'dataObj'!\n")
    }

    varNames <- g_var_name(inputData)
    varNamesInDims <- sapply(1:length(dimList), function(x) {
      g_varname(dimList[[x]])
    })

    if ( !all(varNamesInDims %in% varNames) ) {
      stop("makeProblem::doPrep() mismatch in variable names in 'inputData' and 'inputDims'!\n")
    }

    rawData <- g_raw_data(inputData)

    # variable names in dataObj
    vNamesInData <- g_var_name(inputData)

    # vNames in inputDims
    vNamesInDimList <- sapply(1:length(inputDims), function(x) {
      g_varname(inputDims[[x]])
    })

    # variables not used
    vNotUsed <- setdiff(vNamesInDimList, varNames)
    if ( length(vNotUsed) > 0 ) {
      removeIndex <- match(vNotUsed, vNamesInDimList)
      inputDims <- inputDims[-c(removeIndex)]
      vNamesInDimList <- sapply(1:length(inputDims), function(x) {
        g_varname(inputDims[[x]])
      })

      if ( any(vNamesInDimList != varNames) ) {
        stop("Error: Matching failed!\n")
      }
    }

    posIndex <- match(vNamesInData, vNamesInDimList)
    dimVarInd <- g_dimvar_ind(inputData)
    if ( length(posIndex) < 1 ) {
      stop("Error: matching of variable names failed. Please check 'inputData' and/or 'inputDims'!\n")
    } else {
      if ( any(is.na(posIndex)) ) {
        dimVarInd <- setdiff(dimVarInd, which(is.na(posIndex)))
        vNamesInData <- vNamesInData[dimVarInd]
        inputDims <- inputDims[na.omit(posIndex)]
      } else {
        # correct order
        inputDims <- inputDims[posIndex]
      }
    }

    ss <- list()
    for ( i in seq_along(dimVarInd) ) {
      remove.vals <- FALSE
      remove_ind <- NULL
      if ( !c_has_default_codes(inputDims[[i]], input=rawData[[dimVarInd[i]]]) ) {
        dups <- g_dups(inputDims[[i]])
        if ( length(dups) > 0 ) {
          dupsUp <- g_dups_up(inputDims[[i]])
          for ( k in length(dups):1 ) {
            ind <- which(rawData[[dimVarInd[i]]]==dups[k])
            if ( length(ind) > 0 ) {
              if ( length(which(rawData[[dimVarInd[i]]]==dupsUp[k])) > 0 ) {
                remove.vals <- TRUE
                remove_ind <- c(remove_ind, ind)
              } else {
                rawData[[dimVarInd[i]]][ind] <- dupsUp[k]
              }
            }
          }
          if ( remove.vals ) {
            rawData <- rawData[-unique(remove_ind)]
          }
          s_raw_data(inputData) <- list(rawData)
        }
        ss[[i]] <- c_standardize(inputDims[[i]], input=rawData[[dimVarInd[i]]])
      } else {
        ss[[i]] <- rawData[[dimVarInd[i]]]
      }
      # remove entries in ss[[1]...ss[[i-1]]
      if ( remove.vals ) {
        if ( i > 1 ) {
          for ( z in 1:(i-1)) {
            ss[[z]] <- ss[[z]][-remove_ind]
          }
        }
      }
    }
    strID <- pasteStrVec(as.vector(unlist(ss)), length(posIndex))

    info <- lapply(inputDims, function(x) {
      sum(g_structure(x))
    })
    strInfo <- list()
    for ( i in 1:length(inputDims) ) {
      sumCur <- info[[i]]
      if ( i == 1 ) {
        strInfo[[i]] <- c(1, sumCur)
      } else {
        strInfo[[i]] <- c(1+max(strInfo[[c(i-1)]]), max(strInfo[[c(i-1)]])+sumCur)
      }
    }

    dimInfoObj <- new("dimInfo",
      dimInfo=inputDims,
      strID=strID,
      strInfo=strInfo,
      vNames=vNamesInData,# because of ordering
      posIndex=dimVarInd # because dimVars are re-ordered according to input data!
    )
    return(list(inputData=inputData, dimInfoObj=dimInfoObj))
  }

  for ( i in seq_along(dimList) ) {
    dimList[[i]] <- init.dimVar(input=list(input=dimList[[i]], vName=names(dimList)[i]))
  }

  ## generate inputData from data
  inputData <- init.dataObj(input=list(inputData=data, dimVarInd=dimVarInd, freqVarInd=freqVarInd, numVarInd=numVarInd, weightInd=weightInd,sampWeightInd=sampWeightInd))

  ## check if all variable names listed in inputDims exist in the
  ## specified dimensions of the input data
  varNames <- g_var_name(inputData)
  varNamesInDims <- sapply(1:length(dimList), function(x) {
    g_varname(dimList[[x]])
  })

  if ( !all(varNamesInDims %in% varNames) ) {
    stop("makeProblem:: mismatch in variable names in 'inputData' and 'inputDims'!\n")
  }

  ## calculate the dimInfoObj and eventually recode inputData
  ## (eventually recode rawData slot of inputData if "rawData" contains "wrong" dups)
  out <- doPrep(inputData, dimList)

  ## use output of doPrep() to calculate an object of class "sdcProblem"
  prob <- c_calc_full_prob(input=list(objectA=out$inputData, objectB=out$dimInfoObj))
  prob
}

#' perform primary suppression in \code{\link{sdcProblem-class}}-objects
#'
#' Function \code{\link{primarySuppression}} is used to identify and suppress primary
#' sensitive table cells in \code{\link{sdcProblem-class}} objects.
#' Argument \code{type} allows to select a rule that should be used to identify
#' primary sensitive cells. At the moment it is possible to identify and
#' suppress sensitive table cells using the frequency-rule, the nk-dominance
#' rule and the p-percent rule.
#'
#' @param object a \code{\link{sdcProblem-class}} object
#' @param type character vector of length 1 defining the primary suppression rule. Allowed types are:
#' \itemize{
#' \item \code{freq}: apply frequency rule with parameters \code{maxN} and \code{allowZeros}
#' \item \code{nk}: apply nk-dominance rule with parameters \code{n}, \code{k} and \code{numVarInd}
#' \item \code{p}: apply p-percent rule with parameters \code{p} and \code{numVarInd}
#' \item \code{pq}: apply pq-rule with parameters \code{p} and \code{q}
#' }
#' @param ... parameters used in the identification of primary sensitive cells. Parameters that can be modified|changed are:
#' \itemize{
#' \item \code{maxN}: numeric vector of length 1 used when applying the frequency rule. All cells having counts <= \code{maxN} are set as primary suppressed. The default value of \code{maxN} is 3.
#' \item \code{allowZeros}: logical vector of length 1 specifying if empty cells (count==0) should be considered sensitive when using the frequency rule. The default value of \code{allowZeros} is 'FALSE' so that empty cells are not considered primary sensitive by default.
#' \item \code{p}: numeric vector of length 1 specifying parameter \code{p} that is used when applying the p-percent rule with default value of 80.
#' \item \code{pq}: numeric vector of length 2 specifying parameters \code{p} and \code{q} that are used when applying the pq-rule with the default being c(25, 50).
#' \item \code{n}: numeric vector of length 1 specifying parameter \code{n} that is used when applying the nk-dominance rule. Parameter \code{n} is set to 2 by default.
#' \item \code{k}: numeric vector of length 1 specifying parameter \code{k} that is used when applying the nk-dominance rule. Parameter \code{n} is set to 85 by default.
#' \item \code{numVarInd}: numeric vector of length 1 specifying the index of the numerical variable that should be used to identify cells that are dominated by 2 (p-percent rule) or n (nk-dominance)-rule. If \code{type} is either 'nk', 'p' or 'pq', it is mandatory to specify \code{numVarInd}.
#' }
#' @return a \code{\link{sdcProblem-class}} object
#'
#' @examples
#' # load micro data
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/microData1.RData", sep="")
#' microData <- get(load(fn))
#'
#' # load problem (as it was created in the example in \code{\link{makeProblem}})
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problem.RData", sep="")
#' problem <- get(load(fn))
#'
#' # we have a look at the frequency table by gender and region
#' xtabs(rep(1, nrow(microData)) ~ gender + region, data=microData)
#'
#' # cell with region=='A' and gender=='female' has 2 units contributing to it
#' # this cell should be considered sensitive according the the freq-rule with 'maxN' equal to 2!
#' p1 <- primarySuppression(problem, type='freq', maxN=2)
#'
#' # we can also apply a p-percent rule with parameter 'p' being 30 as below.
#' # This is only possible if we are dealing with micro data and we also have to specify the index of
#' # a numeric variable.
#' p2 <- primarySuppression(problem, type='p', p=30, numVarInd=1)
#'
#' # looking at anonymization states we see, that one cell is primary suppressed (sdcStatus=='u')
#' # and the remaining cells are possible candidates for secondary suppression (sdcStatus=='s') given
#' # the frequency rule with parameter 'maxN=2'.
#' # Applying the p-percent rule with parameter 'p=30' resulted in two primary suppressions.
#' data.frame(p1.sdc=getInfo(p1, type='sdcStatus'), p2.sdc=getInfo(p2, type="sdcStatus"))
#'
#' @rdname primarySuppression
#' @export primarySuppression
#' @note the nk-dominance rule, the p-percent rule and the pq-rule can only be applied if micro data have been used as input data to function \code{\link{makeProblem}}.
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
primarySuppression <- function(object, type, ...) {
  start.time <- proc.time()
  if ( !type %in% c('nk', 'freq', 'p', 'pq') ) {
    stop("valid types are 'nk', 'freq', 'p' or 'pq'!\n")
  }

  numVarsIndices <- g_numvar_ind(g_dataObj(object))
  paraList <- genParaObj(selection='control.primary', numVarIndices=numVarsIndices, ...)

  if ( type == "freq") {
    object <- c_rule_freq(object, input=paraList)
  }

  if ( type == "nk" ) {
    if ( is.na(paraList$numVarInd) ) {
      stop("argument 'numVarInd' must be specified!\n")
    }
    object <- c_rule_nk(object, input=paraList)
  }

  if ( type == "p") {
    if ( is.na(paraList$numVarInd) ) {
      stop("argument 'numVarInd' must be specified!\n")
    }
    object <- c_rule_p(object, input=paraList)
  }

  if ( type == "pq") {
    if ( is.na(paraList$numVarInd) ) {
      stop("argument 'numVarInd' must be specified!\n")
    }
    object <- c_rule_pq(object, input=paraList)
  }

  elapsed.time <- g_elapsedTime(object) + (proc.time() - start.time)[3]
  s_elapsedTime(object) <- elapsed.time
  return(object)
}

#' protecting \code{\link{sdcProblem-class}} objects
#'
#' Function \code{\link{protectTable}} is used to protect primary sensitive table cells
#' (that usually have been identified and set using
#' \code{\link{primarySuppression}}). The function protects primary
#' sensitive table cells according to the method that has been chosen and the
#' parameters that have been set. Additional parameters that are used to control
#' the protection algorithm are set using parameter \code{...}.
#'
#' @param object a \code{\link{sdcProblem-class}} object that has created using \code{\link{makeProblem}} and has been modified by \code{\link{primarySuppression}}
#' @param method a character vector of length 1 specifying the algorithm that should be used to protect the primary sensitive table cells. Allowed values are:
#' \itemize{
#' \item \code{OPT}: protect the complete problem at once using a cut and branch algorithm. The optimal algorithm should be used for small problem-instances only.
#' \item \code{HITAS}: split the overall problem in smaller problems. These problems are protected using a top-down approach.
#' \item \code{HYPERCUBE}: protect the complete problem by protecting sub-tables with a fast heuristic that is based on finding and suppressing geometric structures (n-dimensional cubes) that are required to protect primary sensitive table cells.
#' \item \code{SIMPLEHEURISTIC}: heuristic, quick procedure which might be applied to very large problem instances
#' }
#' @param ... parameters used in the protection algorithm that has been selected. Parameters that can be changed are:
#' \itemize{
#' \item general parameters include:
#' \itemize{
#' \item \code{verbose}: logical vector of length 1 defining if verbose output should be produced. Parameter \code{verbose} defaults to 'FALSE'
#' \item \code{save}: logical vector of length 1 defining if temporary results should be saved in the current working directory (TRUE) or not (FALSE). Parameter \code{save} defaults to 'FALSE' }
#' \item parameters used for HITAS|OPT procedures:
#' \itemize{
#' \item \code{solver}: character vector of length 1 defining the solver to be used. Currently available choices are limited to 'glpk'.
#' \item \code{timeLimit}: numeric vector of length 1 (or NULL) defining a time limit in minutes after which the cut and branch algorithm should stop and return a possible non-optimal solution. Parameter \code{safe} has a default value of 'NULL'
#' \item \code{maxVars}: a numeric vector of length 1 (or NULL) defining the maximum problem size in terms of decision variables for which an optimization should be tried. If the number of decision variables in the current problem are larger than parameter \code{maxVars}, only a possible non-optimal, heuristic solution is calculated. Parameter \code{safe} has a default value of 'NULL'
#' \item \code{fastSolution}: logical vector of length 1 defining if or if not the cut and branch algorithm will be started or if the possibly non-optimal heuristic solution is returned independent of parameter \code{maxVars}. Parameter \code{fastSolution} has a default value of 'FALSE'
#' \item \code{fixVariables}: logical vector of length 1 defining whether or not it should be tried to fix some variables to zero or one based on reduced costs early in the cut and branch algorithm. Parameter \code{fixVariables} has a default value of 'TRUE'
#' \item \code{approxPerc}: numeric vector of length 1 that defines a percentage for which a integer solution of the cut and branch algorithm is accepted as optimal with respect to the upper bound given by the (relaxed) solution of the master problem. Its default value is set to '10'
#' \item \code{useC}: boolean vector of length 1 defining if c++ implementation of the secondary cell suppression problem should be used, defaults to FALSE}
#' \item parameters used for HYPERCUBE procedure:
#' \itemize{
#' \item \code{protectionLevel}: numeric vector of length 1 specifying the required protection level for the HYPERCUBE-procedure. Its default value is 80
#' \item \code{suppMethod}: character vector of length 1 defining the rule on how to select the 'optimal' cube to protect a single sensitive cells. Possible choices are:
#' \itemize{
#' \item \code{minSupps}: minimize the number of additional secondary suppressions (this is also the default setting).
#' \item \code{minSum}: minimize the sum of counts of additional suppressed cells
#' \item \code{minSumLogs}: minimize the log of the sum of additional suppressed cells}
#' \item suppAdditionalQuader: logical vector of length 1 specfifying if additional cubes should be suppressed if any secondary suppressions in the 'optimal' cube are 'singletons'. Parameter \code{suppAdditionalQuader} has a default value of 'FALSE'}
#' \item parameter used for protectLinkedTables():
#' \itemize{
#' \item \code{maxIter}: numeric vector of length 1 specifying the maximal number of interations that should be make while trying to protect common cells of two different tables. The default value of parameter \code{maxIter} is 10}
#' }
#'
#' @return an \code{\link{safeObj-class}} object
#' @examples
#' # load problem (as it was created after performing primary suppression
#' # in the example of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # protect the table using the 'HITAS' algorithm with verbose output
#' protectedData <- protectTable(problem, method='HITAS', verbose=TRUE, useC=TRUE)
#'
#' # showing a summary
#' summary(protectedData)
#'
#' # looking at the final table with result suppression pattern
#' print(getInfo(protectedData, type='finalData'))
#' @rdname protectTable
#' @export protectTable
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
protectTable <- function(object, method, ...) {
  if ( !method %in% c('HITAS', 'OPT', 'HYPERCUBE', 'SIMPLEHEURISTIC') ) {
    stop("valid methods are 'SIMPLEHEURISTIC', 'HITAS', 'HYPERCUBE' or 'OPT'!\n")
  }

  paraList <- genParaObj(selection='control.secondary', method=method, ...)
  if ( length(g_primSupps(object@problemInstance)) == 0 ) {
    return(c_finalize(object=object, input=paraList))
  }

  if ( method == 'SIMPLEHEURISTIC' ) {
    out <- c_quick_suppression(object, input=paraList)
    out <- out$object
  } else {
    if ( paraList$useC ) {
      if ( method == "OPT" ) {
        out <- c_opt_cpp(object=object, input=paraList)
      }
      if ( method == "HITAS" ) {
        out <- c_hitas_cpp(object=object, input=paraList)
      }
    } else {
      out <- c_anon_worker(object, input=paraList)
    }
  }
  invisible(c_finalize(object=out, input=paraList))
}

#' attacking primary suppressed cells and calculating current lower and upper bounds
#'
#' Function \code{\link{attack}} is used to calculate lower and upper bounds for a given
#' sdcProblem (stored as object of class \code{\link{sdcProblem-class}}).
#' For all calculations the current suppression pattern is used when calculating solutions of the
#' attacker's problem.
#'
#' @param object an object of class \code{\link{sdcProblem-class}}
#' @param verbose a logical vector specifying if output should be verbose (TRUE) or not (FALSE)
#' @return a data.frame with column 'index' holding indices of primary suppressed cells and columns
#' 'bounds_min' and 'bounds_max' featuring calculated lower and upper bounds for each cell.
#' Column 'protected' shows if a given cell is accordingly protected (TRUE) or not (FALSE).
#'
#' @examples
#' # load problem (as it was created after performing primary suppression
#' # in the example of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # calculate current lower|upper bounds given current suppression pattern
#' # (in this case consisting of primary suppressions only)
#' attack(problem, verbose=FALSE)
#' @rdname attack
#' @export attack
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
attack <- function(object, verbose=FALSE) {
  out <- csp_cpp(sdcProblem=object, attackonly=TRUE, verbose=verbose)
  return(out)
}

#' query information from objects
#'
#' Function \code{\link{getInfo}} is used to query information from objects of class
#' \code{\link{sdcProblem-class}}, \code{\link{problemInstance-class}} or \code{\link{safeObj-class}}
#'
#' @param object a \code{\link{sdcProblem-class}} object, \code{\link{problemInstance-class}} object or \code{\link{safeObj-class}} object.
#' @param type a character vector of length 1 specifying the information which should be returned.
#' \itemize{
#' \item if argument \code{object} is of class \code{sdcProblem-class} or \code{\link{problemInstance-class}}, valid choices are:
#' \itemize{
#' \item \code{lb}: slot 'lb' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{ub}: slot 'ub' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{LPL}: slot 'LPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{SPL}: slot 'SPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{UPL}: slot 'UPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{sdcStatus}:  slot 'sdcStatus' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{freq}: slot 'freq' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{strID}: slot 'strID' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{numVars}: slot 'numVars' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{w}: slot 'w' of input \code{object} if it is of class \code{\link{problemInstance-class}} or this slot within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}} }
#' \item if argument \code{object} is of class \code{\link{safeObj-class}}, valid choices are:
#' \itemize{
#' \item \code{finalData}: slot 'finalData' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{nrNonDuplicatedCells}: slot 'nrNonDuplicatedCells' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{nrPrimSupps}: slot 'nrPrimSupps' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{nrSecondSupps}: slot 'nrSecondSupps' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{nrPublishableCells}: slot 'nrPublishableCells' of input \code{object} of class \code{\link{safeObj-class}}
#' \item \code{suppMethod}: slot 'suppMethod' of input \code{object} of class \code{\link{safeObj-class}}}
#' }
#'
#' @return manipulated data dependend on arguments \code{object} and \code{type}
#'
#' @examples
#' # load problem (as it was created in the example
#' # of \code{\link{makeProblem}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problem.RData", sep="")
#' problem <- get(load(fn))
#'
#' # problem is an object of class \code{\link{sdcProblem-class}}
#' print(class(problem))
#'
#' for ( slot in c('lb','ub','LPL','SPL','UPL','sdcStatus',
#'   'freq', 'strID', 'numVars', 'w') ) {
#'   cat('slot', slot,':\n')
#'   print(getInfo(problem, type=slot))
#' }
#'
#' # extracting information for objects of class \code{\link{safeObj-class}}
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/protectedData.RData", sep="")
#' protectedData <- get(load(fn))
#' for ( slot in c('finalData', 'nrNonDuplicatedCells', 'nrPrimSupps',
#'   'nrSecondSupps', 'nrPublishableCells', 'suppMethod') ) {
#'   cat('slot', slot,':\n')
#'   print(getInfo(protectedData, type=slot))
#' }
#' @rdname getInfo
#' @export getInfo
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
getInfo <- function(object, type) {
  if ( !class(object) %in% c('sdcProblem', 'problemInstance', 'safeObj') ) {
    stop("getInfo:: argument 'object' must be of class 'sdcProblem', 'problemInstance' or 'safeObj'!\n")
  }

  if ( class(object) == 'safeObj' ) {
    if ( !type %in% c('finalData', 'nrNonDuplicatedCells', 'nrPrimSupps', 'nrSecondSupps', 'nrPublishableCells', 'suppMethod') ) {
      stop("getInfo:: type must be one of 'finalData', 'nrNonDuplicatedCells', 'nrPrimSupps', 'nrSecondSupps', 'nrPublishableCells' or 'suppMethod'!\n")
    }
    return(get.safeObj(object, type=type, input=list()))
  }
  else {
    if ( !type %in% c('lb','ub','LPL','SPL','UPL','sdcStatus', 'freq', 'strID', 'numVars', 'w') ) {
      stop("getInfo:: check argument 'type'!\n")
    }
    if ( class(object) == 'sdcProblem' ) {
      pI <- g_problemInstance(object)
    } else {
      pI <- object
    }
    return(get.problemInstance(pI, type=type))
  }
}

#' set information of \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} objects
#'
#' Function \code{\link{getInfo}} is used to query information from
#' \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} objects
#'
#' @param object an object of class \code{\link{sdcProblem-class}} or \code{\link{problemInstance-class}}
#' @param type a character vector of length 1 specifying the the information that should be changed or modified, valid choices are:
#' \itemize{
#' \item \code{lb}: slot 'lb' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'lb' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{ub}: slot 'ub' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'ub' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{LPL}: slot 'LPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'LPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{SPL}: slot 'SPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'SPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{UPL}: slot 'UPL' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'UPL' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}}
#' \item \code{sdcStatus}:  slot 'sdcStatus' of input \code{object} if it is of class \code{\link{problemInstance-class}} or slot 'sdcStatus' within slot 'problemInstance' if \code{object} is of class \code{\link{sdcProblem-class}} }
#' @param index numeric vector defining cell-indices for which which values in a specified slot should be changed|modified
#' @param input numeric or character vector depending on argument \code{type} with its length matching the length of argument \code{index}
#' \itemize{
#' \item character vector if type matches 'sdcStatus'
#' \item a numeric vector if type matches 'lb', 'ub', 'LPL', 'SPL' or 'UPL'
#' }
#'
#' @return a \code{\link{sdcProblem-class}}- or \code{\link{problemInstance-class}} object
#'
#' @examples
#' # load primary suppressed data (created in the example of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # which is the overall total?
#' index.tot <- which.max(getInfo(problem, 'freq'))
#' index.tot
#'
#' # we see that the cell with index.tot==1 is the overall total and its
#' # anonymization state of the total can be extracted as follows:
#' print(getInfo(problem, type='sdcStatus')[index.tot])
#'
#' # we want this cell to never be suppressed
#' problem <- setInfo(problem, type='sdcStatus', index=index.tot, input='z')
#'
#' # we can verify this:
#' print(getInfo(problem, type='sdcStatus')[index.tot])
#'
#' # changing slot 'UPL' for all cells
#' inp <- data.frame(strID=getInfo(problem,'strID'), UPL_old=getInfo(problem,'UPL'))
#' inp$UPL_new <- inp$UPL_old+1
#' problem <- setInfo(problem, type='UPL', index=1:nrow(inp), input=inp$UPL_new)
#'
#' @rdname setInfo
#' @export setInfo
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setInfo <- function(object, type, index, input) {
  if ( !class(object) %in% c('sdcProblem', 'problemInstance') ) {
    stop("setInfo:: argument 'object' must be of class 'sdcProblem' or 'problemInstance'!\n")
  }

  if ( !type %in% c('lb','ub','LPL','SPL','UPL','sdcStatus') ) {
    stop("setInfo:: check argument 'type'!\n")
  }

  if ( class(object) == "sdcProblem" ) {
    pI <- g_problemInstance(object)
  } else {
    pI <- object
  }

  pI <- set.problemInstance(pI, type=type, input=list(index=index, values=input))

  if ( class(object) == "sdcProblem" ) {
    s_problemInstance(object) <- pI
  } else {
    object <- pI
  }
  object
}

#' change anonymization status of a specific cell
#'
#' Function \code{\link{changeCellStatus}} allows to change|modify the anonymization state
#' of single table cells for objects ofs class \code{\link{sdcProblem-class}}.
#'
#' @param object an object of class \code{\link{sdcProblem-class}}
#' @param characteristics a character vector specifying characteristics of the table cell that should be identified for each dimensional variable defining the table
#' @param varNames a character vector specifying variable names of dimensional variables defining the tables
#' @param rule character vector of length 1 specifying a valid anonymization code ('u', 'z', 'x', 's') to which the the cell under consideration should be set.
#' @param verbose logical vector of length 1 defining verbosity, defaults to 'FALSE'
#'
#' @return a \code{\link{sdcProblem-class}} object
#'
#' @examples
#' # load primary suppressed data (as created in the example
#' # of \code{\link{primarySuppression}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/problemWithSupps.RData", sep="")
#' problem <- get(load(fn))
#'
#' # we want to mark the cell region='D' and gender='male' primary sensitive
#' characteristics <- c('D', 'male')
#' varNames <- c('region', 'gender')
#' verbose <- TRUE
#' rule <- 'u'
#'
#' # looking at the distribution of anonymization states before...
#' print(table(getInfo(problem, 'sdcStatus')))
#'
#' # setting the specific cell as primary sensitive
#' problem <- changeCellStatus(problem, characteristics, varNames, rule, verbose)
#'
#' # having a second look at the anonymization states
#' print(table(getInfo(problem, 'sdcStatus')))
#'
#' @rdname changeCellStatus
#' @export changeCellStatus
#' @note Important: the \code{i}-th element of argument \code{characteristics} is uses as the desired characteristic for the dimensional variable specified at the \code{i}-th position of argument \code{varNames}!
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
changeCellStatus <- function(object, characteristics, varNames, rule, verbose=FALSE) {
  if ( class(object) != 'sdcProblem' ) {
    stop("changeCellStatus:: argument 'object' must be of class 'sdcProblem'!\n")
  }

  paraList <- list()
  paraList$names <- varNames
  paraList$codes <- characteristics
  paraList$verbose <- verbose

  cellID <- c_cellID(object, input=paraList)

  pI <- g_problemInstance(object)
  s_sdcStatus(pI) <- list(index=cellID, vals=rule)
  s_problemInstance(object) <- pI

  if ( paraList$verbose ) {
    cat('--> The cell with ID=', cellID,'and Frequency',g_freq(pI)[cellID], 'has been set to', rule,'.\n')
  }
  object
}

#' query information for a specific cell in \code{\link{safeObj-class}} objects
#'
#' Function \code{\link{cellInfo}} is used to query information for a single table cell
#' for objects of class \code{\link{safeObj-class}}.
#'
#' @param object an object of class \code{\link{safeObj-class}}
#' @param characteristics a character vector specifying characteristics of the table cell that should be identified for each dimensional variable defining the table
#' @param varNames a character vector specifying variable names of dimensional variables defining the tables
#' @param verbose logical vector of length 1 defining verbosity, defaults to 'FALSE'
#'
#' @return a list containing the following calculated information
#' \itemize{
#' \item \code{cellID}: numeric vector of length 1 specifying the index of the cell within the final result dataset
#' \item \code{data}: a data.frame containing a single row with the index of the table cell of interest
#' \item \code{primSupp}: logical vector of length 1 that is 'TRUE' if the cell is a primary sensitive cell and 'FALSE' otherwise
#' \item \code{secondSupp}: logical vector of length 1 that is 'TRUE' if the cell is a secondary suppressed cell and 'FALSE' otherwise
#' }
#'
#' @examples
#' # load protected data (as created in the example
#' # of \code{\link{protectTable}})
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/protectedData.RData", sep="")
#' protectedData <- get(load(fn))
#' characteristics <- c('male', 'D')
#' varNames <- c('gender', 'region')
#' info <- cellInfo(protectedData, characteristics, varNames, verbose=FALSE)
#'
#' # show the info about this cell
#' str(info)
#'
#' @rdname cellInfo
#' @export cellInfo
#' @note Important: the \code{i}-th element of argument \code{characteristics} is uses as the desired characteristic for the dimensional variable specified at the \code{i}-th position of argument \code{varNames}!
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
cellInfo <- function(object, characteristics, varNames, verbose=FALSE) {
  paraList <- list()
  paraList[[1]] <- varNames
  paraList[[2]] <- characteristics
  paraList[[3]] <- verbose
  g_getCellInfo(object, input=paraList)
}

#' protect two \code{\link{sdcProblem-class}} objects that have common cells
#'
#' \code{\link{protectLinkedTables}} can be used to protect tables, that have
#' common cells. It is of course required that after the anonymization process
#' has finished, all common cells have the same anonymization state in both
#' tables.
#'
#' @param objectA a \code{\link{sdcProblem-class}} object
#' @param objectB a \code{\link{sdcProblem-class}} object
#' @param commonCells a list object defining common cells in \code{objectA} and \code{objectB}. For each variable that has one or more common codes in both tables, a list element needs to be specified.
#' \itemize{
#' \item List-elements of length 3: Variable has exact same levels and structure in both tables
#' \itemize{
#' \item \code{first element}: character vector of length 1 specifying the variable name in argument \code{objectA}
#' \item \code{second element}: character vector of length 1 specifying the variable name in argument \code{objectB}
#' \item \code{third element}: character vector of length 1 being with keyword \code{ALL} }
#' \item List-elements of length 4: Variable has different codes and levels in tables \code{objectA} and \code{objectB}
#' \itemize{
#' \item \code{first element}: character vector of length 1 specifying the variable name in argument \code{objectA}
#' \item \code{second element}: character vector of length 1 specifying the variable name in argument \code{objectB}
#' \item \code{third element}: character vector defining codes within \code{objectA}
#' \item \code{fourth element}: character vector with length that equals the length of the third list-element. The vector defines codes of the variable in \code{objectB} that match the codes given in the third list-element for \code{objectA}.
#' }
#' }
#' @param method a character vector of length 1 specifying the algorithm that should be used to protect the primary sensitive table cells. Allowed values are:
#' \itemize{
#' \item \code{HITAS}:
#' \item \code{SIMPLEHEURISTIC}:
#' \item \code{OPT}: }
#' @param ... additional arguments to control the secondary cell suppression algorithm. For details, see \code{\link{protectTable}}.
#'
#' @return a list of length 2 with each list-element being an \code{\link{safeObj-class}} object
#'
#' @examples
#' \dontrun{
#' # load micro data for further processing
#' sp <- searchpaths()
#' fn <- paste(sp[grep("sdcTable", sp)], "/data/microData2.RData", sep="")
#' microData <- get(load(fn))
#'
#' # table1: defined by variables 'gender' and 'ecoOld'
#' microData1 <- microData[,c(2,3,5)]
#'
#' # table2: defined by variables 'region', 'gender' and 'ecoNew'
#' microData2 <- microData[,c(1,2,4,5)]
#'
#' # we need to create information on the hierarchies
#' # variable 'region': exists only in microDat2
#' dim.region <- data.frame(h=c('@@','@@@@','@@@@'), l=c('Tot', 'R1','R2'))
#'
#' # variable 'gender': exists in both datasets
#' dim.gender <- data.frame(h=c('@@','@@@@','@@@@'), l=c('Tot', 'm','f'))
#'
#' # variable 'ecoOld': exists only in microDat1
#' dim.ecoOld <- data.frame(
#'  h=c('@@','@@@@','@@@@@@','@@@@@@','@@@@','@@@@@@','@@@@@@'),
#'  l=c('Tot','A','Aa','Ab','B','Ba','Bb'))
#'
#' # variable 'ecoNew': exists only in microDat2
#' dim.ecoNew <- data.frame(
#'  h=c('@@','@@@@','@@@@@@','@@@@@@','@@@@@@','@@@@','@@@@@@','@@@@@@','@@@@@@'),
#'  l=c('Tot','C','Ca','Cb','Cc','D','Da','Db','Dc'))
#'
#' # creating objects holding information on dimensions
#' dimList1 <- list(gender=dim.gender, ecoOld=dim.ecoOld)
#' dimList2 <- list(region=dim.region, gender=dim.gender, ecoNew=dim.ecoNew)
#'
#' # creating input objects for further processing. For details have a look at
#' # \code{\link{makeProblem}}.
#' problem1 <- makeProblem(data=microData1, dimList=dimList1, dimVarInd=c(1,2),
#'          numVarInd=3)
#' problem2 <- makeProblem(data=microData2, dimList=dimList2, dimVarInd=c(1,2,3),
#'          numVarInd=4)
#'
#' # the cell specified by gender=='Tot' and ecoOld=='A'
#' # is one of the common cells! -> we mark it as primary suppression
#' problem1 <- changeCellStatus(problem1, characteristics=c('Tot', 'A'),
#'      varNames=c('gender','ecoOld'), rule='u', verbose=FALSE)
#'
#' # the cell specified by region=='Tot' and gender=='f' and ecoNew=='C'
#' # is one of the common cells! -> we mark it as primary suppression
#' problem2 <- changeCellStatus(problem2, characteristics=c('Tot', 'f', 'C'),
#'  varNames=c('region','gender', 'ecoNew'), rule='u', verbose=FALSE)
#'
#' # specifying input to define common cells
#' commonCells <- list()
#'
#' # variable "gender"
#' commonCells$v.gender <- list()
#' commonCells$v.gender[[1]] <- 'gender' # variable name in 'problem1'
#' commonCells$v.gender[[2]] <- 'gender' # variable name in 'problem2'
#' # 'gender' has equal characteristics on both datasets -> keyword 'ALL'
#' commonCells$v.gender[[3]] <- 'ALL'
#'
#' # variable: ecoOld and ecoNew
#' commonCells$v.eco <- list()
#' commonCells$v.eco[[1]] <- 'ecoOld'   # variable name in 'problem1'
#' commonCells$v.eco[[2]] <- 'ecoNew'   # variable name in 'problem2'
#'
#' # vector of common characteristics: A and B in variable 'ecoOld' in 'problem1'
#' commonCells$v.eco[[3]] <- c("A","B")
#' # correspond to characteristics 'C' and 'D' in variable 'ecoNew' in 'problem2'
#' commonCells$v.eco[[4]] <- c("C","D")
#'
#' # protect the linked data
#' result <- protectLinkedTables(problem1, problem2, commonCells, method='HITAS', verbose=TRUE)
#'
#' # having a look at the results
#' result.tab1 <- result[[1]]
#' result.tab2 <- result[[2]]
#' summary(result.tab1)
#' summary(result.tab2)
#' }
#'
#' @rdname protectLinkedTables
#' @export protectLinkedTables
#' @seealso \code{\link{protectTable}}
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
protectLinkedTables <- function(objectA, objectB, commonCells, method, ...) {
  f.calcCommonCellIndices <- function(input1, input2, commonCells) {
    id1 <- id2 <- NULL
    dat1 <- g_df(input1); dat1$id1 <- 1:nrow(dat1)
    dat2 <- g_df(input2); dat2$id2 <- 1:nrow(dat2)
    dat1$strID <- dat1$sdcStatus <- NULL
    dat2$strID <- dat2$sdcStatus <- NULL
    dI1 <- g_dimInfo(input1)
    dI2 <- g_dimInfo(input2)

    # restrict to totals in non-overlapping variables in dataset1
    totvars <- setdiff(dI1@vNames, sapply(commonCells, function(x) x[[1]] ))
    if ( length(totvars) > 0 ) {
      for ( i in 1:length(totvars)) {
        cmd <- paste0("dat1 <- dat1[",totvars[i],"=='",dI1@dimInfo[[totvars[i]]]@codesDefault[1],"']")
        eval(parse(text=cmd))
      }
    }
    # restrict to totals in non-overlapping variables in dataset2
    totvars <- setdiff(dI2@vNames, sapply(commonCells, function(x) x[[2]] ))
    if ( length(totvars) > 0 ) {
      for ( i in 1:length(totvars)) {
        cmd <- paste0("dat2 <- dat2[",totvars[i],"=='",dI2@dimInfo[[totvars[i]]]@codesDefault[1],"']")
        eval(parse(text=cmd))
      }
    }

    for ( i in 1:length(commonCells)) {
      if ( length(commonCells[[i]]) == 4 ) {
        cmd1 <- paste0("dat1 <- dat1[,tmpxxx_V",i,":=",commonCells[[i]][[1]],"_o]")
        cmd2 <- paste0("dat2 <- dat2[,tmpxxx_V",i,":=",commonCells[[i]][[2]],"_o]")
        eval(parse(text=cmd1))
        eval(parse(text=cmd2))

        # recode different codes to those of dataset1
        codes <- commonCells[[i]][[4]]
        codesX <- commonCells[[i]][[3]]
        for  ( z in 1:length(codes)) {
          if ( codes[z] != codesX[z] ) {
            cmd <- paste0("dat2 <- dat2[tmpxxx_V",i,"=='",codes[z],"',tmpxxx_V",i,":='",codesX[z],"']")
            eval(parse(text=cmd))
          }
        }
      } else {
        # nothing to do, codes are the same!
        cmd1 <- paste0("dat1[,tmpxxx_V",i,":=",commonCells[[i]][[1]],"_o]")
        cmd2 <- paste0("dat2[,tmpxxx_V",i,":=",commonCells[[i]][[2]],"_o]")
        eval(parse(text=cmd1))
        eval(parse(text=cmd2))
      }
    }

    kV1 <- names(dat1)[grep("tmpxxx", names(dat1))]; setkeyv(dat1, kV1)
    kV2 <- names(dat2)[grep("tmpxxx", names(dat2))]; setkeyv(dat2, kV2)
    mm <- merge(dat1, dat2); setkey(mm, id1)
    if ( any(mm$freq.x != mm$freq.y) ) {
      stop("Error: common cells must have same values!\n")
    }
    return(list(commonInd1=mm$id1, commonInd2=mm$id2))
  }

  f.checkCommonCells <- function(suppPattern1, suppPattern2, commonCellIndices) {
    indOK <- TRUE
    if ( any(suppPattern1[commonCellIndices[[1]]] != suppPattern2[commonCellIndices[[2]]]) )
      indOK <- FALSE
    return(indOK)
  }

  ### arguments ###
  if ( !method %in% c('HITAS', 'OPT', 'SIMPLEHEURISTIC') ) {
    stop("valid methods are 'HITAS', 'OPT' or 'SIMPLEHEURISTIC'!\n")
  }

  paraList <- genParaObj(selection='control.secondary', method=method, ...)

  ### first run
  if ( method == "SIMPLEHEURISTIC" ) {
    outA <- c_quick_suppression(objectA, input=paraList)$object
    outB <- c_quick_suppression(objectB, input=paraList)$object
  } else {
    if ( paraList$useC ) {
      if ( method == "OPT" ) {
        outA <- c_opt_cpp(object=objectA, input=paraList)
        outB <- c_opt_cpp(object=objectB, input=paraList)
      }
      if ( method == "HITAS" ) {
        outA <- c_hitas_cpp(object=objectA, input=paraList)
        outB <- c_hitas_cpp(object=objectB, input=paraList)
      }
    } else {
      outA <- c_anon_worker(object=objectA, input=paraList)
      outB <- c_anon_worker(object=objectB, input=paraList)
    }
  }

  pI.A <- g_problemInstance(outA)
  pI.B <- g_problemInstance(outB)
  # calc original primary suppressions

  origPrimSupp1Index <- g_primSupps(pI.A)
  origPrimSupp2Index <- g_primSupps(pI.B)

  # no primary suppressions
  if ( length(origPrimSupp1Index) + length(origPrimSupp2Index) == 0 ) {
    if ( paraList$verbose ) {
      cat("\n===> no primary suppressions. All common cells have the same anonymity-status! [Finished]\n")
    }
    outA <- c_finalize(object=outA, input=paraList)
    outB <- c_finalize(object=outB, input=paraList)
    return(list(outObj1=outA, outObj2=outB))
  }

  # calculate commonCells:
  commonCellIndices <- f.calcCommonCellIndices(outA, outB, commonCells)

  # suppression patterns after the first run
  suppPatternA <- g_suppPattern(pI.A)
  suppPatternB <- g_suppPattern(pI.B)

  indOK <- f.checkCommonCells(suppPatternA, suppPatternB, commonCellIndices)
  counter <- 1
  if ( !indOK ) {
    if ( paraList$verbose ) {
      cat("we have to start the iterative procedure!\n")
    }
    runInd <- TRUE
    while ( runInd ) {
      x <- cbind(suppPatternA[commonCellIndices[[1]]], suppPatternB[commonCellIndices[[2]]])
      index <- list()
      i1 <- which(x[,1] == 0 & x[,2]==1)
      i2 <- which(x[,1] == 1 & x[,2]==0)
      index[[1]] <- commonCellIndices[[1]][i1]
      index[[2]] <- commonCellIndices[[2]][i2]

      for ( j in 1:2 ) {
        if ( length(index[[j]]) > 0 ) {
          if ( j == 1 ) {
            pI.A <- g_problemInstance(outA)
            s_sdcStatus(pI.A) <- list(index=index[[j]], vals=rep("u", length(index[[j]])))
            s_problemInstance(outA) <- pI.A
            s_indicesDealtWith(outA) <- NULL
            s_startJ(outA) <- 1
            s_startI(outA) <- 1
            outA <- c_anon_worker(outA, input=paraList)
          } else {
            pI.B <- g_problemInstance(outB)
            s_sdcStatus(pI.B) <- list(index=index[[j]], vals=rep("u", length(index[[j]])))
            s_problemInstance(outB) <- pI.B
            s_indicesDealtWith(outB) <- NULL
            s_startJ(outB) <- 1
            s_startI(outB) <- 1
            outB <- c_anon_worker(outB, input=paraList)
          }
        }
      }

      suppPatternA <- g_suppPattern(g_problemInstance(outA))
      suppPatternB <- g_suppPattern(g_problemInstance(outB))

      cbind(suppPatternA[commonCellIndices[[1]]], suppPatternB[commonCellIndices[[2]]])
      indOK <- f.checkCommonCells(suppPatternA, suppPatternB, commonCellIndices)
      if ( indOK )
        runInd <- FALSE
      if ( counter > paraList$maxIter ) {
        runInd <- FALSE
        warning("iterative procedure did not converge! --> returning NULL")
        return(NULL)
      }
      counter <- counter + 1
    }
  }
  if ( paraList$verbose ) {
    cat("\n===> all common cells have the same anonymity-state in both tables after",counter,"iterations! [Finished]\n")
  }
  outA <- c_finalize(object=outA, input=paraList)
  outB <- c_finalize(object=outB, input=paraList)
  return(list(outObj1=outA, outObj2=outB))
}
