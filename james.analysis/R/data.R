
############################# BINARY DATA DOCS ##############################

#' Results of example algorithm comparison
#' 
#' Contains results of an example analysis performed with the 'JAMES' extensions 
#' module. The performance of two algorithms is compared (random descent and 
#' parallel tempering) for a core selection problem in which the mean 
#' entry-to-nearest-entry distance is maximized. Four different data sets have 
#' been analyzed. Details about the performed analysis are provided at the 
#' website (see below).
#' 
#' @format S3 object of class "james", as if produced by 
#'   \code{\link{readJAMES}}.
#' @source \url{http://www.jamesframework.org/examples/#analysis}
#' @seealso \code{\link{readJAMES}}
#'   
#' @examples
#' # load data
#' data(james)
#' summary(james)
#' 
#' # plot convergence curves for coconut data set
#' plotConvergence(james, problem = "coconut", min.time = 1000, max.time = 100000)
#' 
#' # create box plots of solution values (quality) and convergence times
#' boxplot(james, problem = "coconut")
#' boxplot(james, problem = "coconut", type = "time")
#' 
#' # extract solution values and convergence times for parallel tempering and random descent
#' values.pt <- getBestSolutionValues(james, problem = "coconut", search = "Parallel Tempering")
#' times.pt <- getConvergenceTimes(james, problem = "coconut", search = "Parallel Tempering")
#' values.rd <- getBestSolutionValues(james, problem = "coconut", search = "Random Descent")
#' times.rd <- getConvergenceTimes(james, problem = "coconut", search = "Random Descent")
#' 
#' # perform wilcoxon test to compare distributions across algorithms
#' values.test <- wilcox.test(values.pt, values.rd)
#' values.test
#' times.test <- wilcox.test(times.pt, times.rd)
#' times.test
#' 
#' # adjust p-values for multiple testing
#' p.adjust(c(values.test$p.value, times.test$p.value))
#' 
"james"

#################################### I/O ####################################

#' Read analysis results from JSON file
#' 
#' Read results from a JSON file produced by the analysis tools in the 'JAMES' 
#' extensions module.
#' 
#' @param file string: path to a JSON file containing results produced by the 
#'   analysis tools from the 'JAMES' extensions module.
#' @return S3 object of class "james" containing the results of running a number
#'   of searches on a set of problems, where each search has been repeatedly 
#'   applied for a number of runs. Data can be manipulated and extracted using 
#'   the provided functions (see below).
#'   
#' @seealso Example data: \code{\link{james}}.
#'   
#'   Data access and manipulations methods: \code{\link{reduceJAMES}}, 
#'   \code{\link{mergeJAMES}}, \code{\link{getProblems}}, 
#'   \code{\link{getSearches}}, \code{\link{getSearchRuns}}, 
#'   \code{\link{getNumSearchRuns}}, \code{\link{getBestSolutionValues}}, 
#'   \code{\link{getBestSolutions}}, \code{\link{getConvergenceTimes}}.
#'   
#'   Plot functions: \code{\link{plotConvergence}}, \code{\link{boxplot.james}}.
#'   
#' @examples
#' # get path to raw JSON file included in package distribution
#' json.file <- system.file("extdata", "james.json", package = "james.analysis")
#' 
#' # read results from file
#' james <- readJAMES(json.file)
#' summary(james)
#' 
#' # plot convergence curves for coconut data set
#' plotConvergence(james, problem = "coconut", min.time = 1000, max.time = 100000)
#' 
#' # create box plots of solution values (quality) and convergence times
#' boxplot(james, problem = "coconut")
#' boxplot(james, problem = "coconut", type = "time")
#' 
#' # extract solution values and convergence times for parallel tempering and random descent
#' values.pt <- getBestSolutionValues(james, problem = "coconut", search = "Parallel Tempering")
#' times.pt <- getConvergenceTimes(james, problem = "coconut", search = "Parallel Tempering")
#' values.rd <- getBestSolutionValues(james, problem = "coconut", search = "Random Descent")
#' times.rd <- getConvergenceTimes(james, problem = "coconut", search = "Random Descent")
#' 
#' # perform wilcoxon test to compare distributions across algorithms
#' values.test <- wilcox.test(values.pt, values.rd)
#' values.test
#' times.test <- wilcox.test(times.pt, times.rd)
#' times.test
#' 
#' # adjust p-values for multiple testing
#' p.adjust(c(values.test$p.value, times.test$p.value))
#' 
#' @export
readJAMES <- function(file) {
  
  # check input
  if (is.null(file) || is.na(file) || !is.character(file)
      || length(file) != 1 || nchar(file) == 0){
    stop("'file' should be a non-empty string")
  }
  # check if file exists
  if(!file.exists(file)){
    msg <- sprintf("'file' %s does not exist", file)
    stop(msg)
  }
  # read JSON file
  results <- rjson::fromJSON(file = file)
  
  # check structure
  
  unexpectedStructure <- function(){
    msg <- sprintf("'unexpected JSON structure in file' %s", file)
    stop(msg)
  }
  
  # check: non empty list of problems
  if(!is.list(results) || length(results) == 0){
    unexpectedStructure()
  }
  # check problem names
  problem.names <- names(results)
  if(is.null(problem.names) || anyDuplicated(problem.names)){
    unexpectedStructure()
  }
  for(n in problem.names){
    if(is.na(n) || nchar(n) == 0){
      unexpectedStructure()
    }
  }
  # check each problem
  for(problem in results){
    # check: non empty list of searches
    if(!is.list(problem) || length(problem) == 0){
      unexpectedStructure()
    }
    # check search names
    search.names <- names(problem)
    if(is.null(search.names) || anyDuplicated(search.names)){
          unexpectedStructure()
        }
    for(n in search.names){
      if(is.na(n) || nchar(n) == 0){
        unexpectedStructure()
      }
    }
    # check each search
    for(search in problem){
      # check: non empty list of search runs
      if(!is.list(search) || length(search) == 0){
        unexpectedStructure()
      }
      # check each run
      for(run in search){
        # check: list with 2 or 3 elements
        if(!is.list(run) || !any(length(run) == c(2,3))){
          unexpectedStructure()
        }
        # check names
        run.components <- names(run)
        if(is.null(run.components) || anyDuplicated(run.components)){
          unexpectedStructure()
        }
        for(n in run.components){
          if(is.na(n) || nchar(n) == 0){
            unexpectedStructure()
          }
        }
        if(!("times" %in% run.components)
           || !("values" %in% run.components)
           || (length(run) == 3) && !("best.solution" %in% run.components)){
          unexpectedStructure()
        }
        # check that number of times and values agree
        if(length(run$times) != length(run$values)){
          unexpectedStructure()
        }
        # check that times are increasing
        if(is.unsorted(run$times)){
          unexpectedStructure()
        }
        # check that values are increasing or decreasing
        if(is.unsorted(run$values) && is.unsorted(rev(run$values))){
          unexpectedStructure()
        }
      }
    }
  }
  
  # add class name
  class(results) <- c("james", class(results))
  
  return(results)
}

############################### MERGE - REDUCE ############################## 

#' Merge analysis results
#' 
#' Merge results from different analyses. If runs of the same search applied to
#' the same problem are found in both data sets, these runs are merged into a
#' single list. This is a generic S3 method.
#' 
#' @param data1 results from the first analysis (of same class as \code{data2}).
#' @param data2 results from the second analysis (of same class as \code{data1}).
#'   
#' @return merged data (assigned classes are retained)
#'   
#' @export
mergeJAMES <- function(data1, data2){
  UseMethod("mergeJAMES")
}
#' @export
mergeJAMES.james <- function(data1, data2){
  # check input
  if(!inherits(data2, "james")){
    stop("'data2' should also be of class \"james\"")
  }
  # copy first data set
  merged <- data1
  # merge results from second data set into first
  for(problem in getProblems(data2)){
    for(search in getSearches(data2, problem)){
      # retrieve existing search runs, empty list if none
      existing.runs <- tryCatch(
                          getSearchRuns(merged, problem, search),
                          error = function(e){
                            list()
                          }
                       )
      # retrieve new search runs
      new.runs <- getSearchRuns(data2, problem, search)
      # merge runs
      merged.runs <- c(existing.runs, new.runs)
      # update merged data
      if(length(getProblems(merged, filter = glob2rx(problem))) == 0){
        # first search for this problem: add problem
        merged[[problem]] <- list()
      }
      merged[[problem]][[search]] <- merged.runs
    }
  }
  return(merged)
}

#' Reduce analysis results to selected problems and searches
#' 
#' Reduce the given \code{data} by filtering the analyzed problems and applied 
#' searches based on the given list of names or \link{regular expression} 
#' (pattern matching is done with \code{\link{grep}}). This is a generic S3 
#' method.
#' 
#' @param data data object containing the analysis results
#' @param problems \link{regular expression} or list of strings. Only those 
#'   problems that match the regular expression or occur in the list are 
#'   retained.
#' @param searches \link{regular expression} or list of strings. Only those 
#'   searches that match the regular expression or occur in the list are 
#'   retained.
#' @param ... any additional arguments are passed to \code{\link{grep}}.
#'   
#' @return Reduced data set containing only those problems and searches whose
#'   names match the respective regular expression or occur in the respective
#'   list of strings. Assigned classes are retained.
#'   
#' @export
reduceJAMES <- function(data, problems = ".*", searches = ".*", ...){
  UseMethod("reduceJAMES")
}
#' @export
reduceJAMES.james <- function(data, problems = ".*", searches = ".*", ...){
  # determine which problems to drop
  problem.names <- getProblems(data)
  if(is.list(problems)){
    drop.problems <- setdiff(problem.names, problems)
  } else {
    drop.problems <- grep(pattern = problems, problem.names, value = TRUE, invert = TRUE, ...)
  }
  # drop those problems
  data[drop.problems] <- NULL
  # iterate over retained problems to filter searches
  for(problem in getProblems(data)){
    # determine which searches to drop
    search.names <- getSearches(data, problem)
    if(is.list(searches)){
      drop.searches <- setdiff(search.names, searches)
    } else {
      drop.searches <- grep(pattern = searches, search.names, value = TRUE, invert = TRUE, ...)
    }
    # drop those searches
    data[[problem]][drop.searches] <- NULL
    # drop problem if all searches have been dropped
    if(length(data[[problem]]) == 0){
      data[[problem]] <- NULL
    }
  }
  # check if data has been retained
  if(length(data) == 0){
    stop("no data is retained")
  }
  return(data)
}

################################## GETTERS ################################## 

#' Get names of analyzed problems
#' 
#' Extracts the names of the problems for which analysis results are contained 
#' in the given data. This is a generic S3 method.
#' 
#' Problem names are sorted using \link[naturalsort]{naturalsort}. If a 
#' \code{filter} is set, only those problem names matching the given 
#' \link{regular expression} are returned (pattern matching is done with 
#' \code{\link{grep}}).
#' 
#' @param data data object containing the analysis results
#' @param filter \link{regular expression} (optional). Only problem names that 
#'   match the given regex are returned, if any.
#' @param ... any additional arguments are passed to \code{\link{grep}}.
#'   
#' @return Sorted vector of strings containing the names of all analyzed 
#'   problems that occur in the given data and match the applied filter (if
#'   any).
#'   
#' @export
getProblems <- function(data, filter, ...){
  UseMethod("getProblems")
}
#' @export
getProblems.james <- function(data, filter, ...){
  problem.names <- names(data)
  # filter, if set
  if(!missing(filter)){
    problem.names <- grep(pattern = filter, problem.names, value = TRUE, ...)
  }
  # sort
  problem.names <- naturalsort::naturalsort(problem.names)
  return(problem.names)
}

#' Get names of applied searches
#' 
#' Extracts the names of all searches that have been applied to the given 
#' \code{problem}.This is a generic S3 method.
#' 
#' Search names are sorted using \link[naturalsort]{naturalsort}. If the
#' \code{data} contains results for a single problem only, the argument 
#' \code{problem} can be omitted. If a \code{filter} is set, only those search 
#' names matching the given \link{regular expression} are returned (pattern 
#' matching is done with \code{\link{grep}}).
#' 
#' @param data data object containing the analysis results
#' @param problem name of the analyzed problem. Can be omitted if the 
#'   \code{data} contains results for a single problem only.
#' @param filter \link{regular expression} (optional). Only search names that 
#'   match the given regex are returned, if any.
#' @param ... any additional arguments are passed to \code{\link{grep}}.
#'   
#' @return Sorted vector of strings containing the names of all searches that
#'   have been applied to the given problem and match the applied filter (if
#'   any).
#'   
#' @export
getSearches <- function(data, problem, filter, ...){
  UseMethod("getSearches")
}
#' @export
getSearches.james <- function(data, problem, filter, ...){
  # fall back to single problem if no problem specified
  if(missing(problem)){
    problem <- getSingleProblem(data)
  }
  # extract problem data
  problem.data <- data[[problem]]
  # check if data is available
  if(is.null(problem.data)){
    msg <- sprintf("'data' does not contain results for problem \"%s\"", problem)
    stop(msg)
  }
  # get all applied searches
  search.names <- names(problem.data)
  # filter, if set
  if(!missing(filter)){
    search.names <- grep(pattern = filter, search.names, value = TRUE, ...);
  }
  # sort
  search.names <- naturalsort::naturalsort(search.names)
  return(search.names)
}

#' Get search run data
#' 
#' Extract the data corresponding to the subsequent runs of a specific 
#' \code{search} being applied to a specific \code{problem}. This is a generic 
#' S3 method.
#' 
#' If the \code{data} contains results for a single problem only, the argument 
#' \code{problem} can be omitted. Likewise, if for the considered \code{problem}
#' results are available for a single search only, the argument \code{search} 
#' can be omitted.
#' 
#' @param data data object containing the analysis results
#' @param problem name of the analyzed problem. Can be omitted if the 
#'   \code{data} contains results for a single problem only.
#' @param search name of the applied search. Can be omitted if the \code{data} 
#'   contains results for a single search only (for the considered 
#'   \code{problem}).
#'   
#' @return A list containing one element for each search run.
#'   
#'   Each run has at least two elements \code{time} and \code{values}, which are
#'   both numeric vectors. The \code{time} vector indicates when the best 
#'   solution was updated during search and the new best solution's value is 
#'   found at the respective index in \code{values}. Times are expressed in 
#'   milliseconds since starting the search. A time of -1 indicates that the 
#'   search was not yet running, which e.g. occurs when a local search adopts a 
#'   random current solution during initialization. Times are always positive 
#'   (or -1) and increasing. Values are either increasing (in case of
#'   maximization) or decreasing (in case of minimization).
#'   
#'   If contained in the given \code{data}, a run also has an element 
#'   \code{best.solution} representing the final best solution found during that
#'   search run. The last element of \code{values} then indicates the value of 
#'   this best solution. When writing results obtained from the analysis tools 
#'   in the 'JAMES' extensions module to a JSON file, one should provide a JSON 
#'   converter for the solution type of the analyzed problems if it is desired 
#'   that the actual best found solutions are contained in the output file.
#'   
#' @export
getSearchRuns <- function(data, problem, search){
  UseMethod("getSearchRuns")
}
#' @export
getSearchRuns.james <- function(data, problem, search){
  # fall back to single problem if not specified
  if(missing(problem)){
    problem <- getSingleProblem(data)   
  }
  # fall back to single search for considered problem, if not specified
  if(missing(search)){
    search <- getSingleSearch(data, problem)
  }
  # extract search runs
  runs <- data[[problem]][[search]]
  # check if data is available
  if(is.null(runs)){
    msg <- sprintf("'data' does not contain results for search \"%s\" being applied to problem \"%s\"",
                   search,
                   problem)
    stop(msg)
  }
  return(runs)
}

#' Get number of applied search runs
#' 
#' Get the number of applied runs of the given \code{search} when solving the 
#' given \code{problem}. This is a generic S3 method.
#' 
#' If the \code{data} contains results for a single problem only, the argument 
#' \code{problem} can be omitted. Likewise, if -- for the considered 
#' \code{problem} -- results are available for a single search only, the 
#' argument \code{search} can be omitted.
#' 
#' @param data data object containing the analysis results
#' @param problem name of the analyzed problem. Can be omitted if the 
#'   \code{data} contains results for a single problem only.
#' @param search name of the applied search. Can be omitted if the \code{data} 
#'   contains results for a single search only (for the considered 
#'   \code{problem}).
#'   
#' @return numeric: number of applied search runs
#'   
#' @export
getNumSearchRuns <- function(data, problem, search){
  UseMethod("getNumSearchRuns")
}
#' @export
getNumSearchRuns.james <- function(data, problem, search){
  # extract search runs
  runs <- getSearchRuns(data, problem, search)
  # count runs
  num.runs <- length(runs)
  return(num.runs)
}

#' Get values of best found solutions
#' 
#' Get the values of the best found solutions during all runs of the given 
#' \code{search} applied to the given \code{problem}. This is a generic S3 
#' method.
#' 
#' If the \code{data} contains results for a single problem only, the argument 
#' \code{problem} can be omitted. Likewise, if for the considered \code{problem}
#' results are available for a single search only, the argument \code{search} 
#' can be omitted.
#' 
#' @param data data object containing the analysis results
#' @param problem name of the analyzed problem. Can be omitted if the 
#'   \code{data} contains results for a single problem only.
#' @param search name of the applied search. Can be omitted if the \code{data} 
#'   contains results for a single search only (for the considered 
#'   \code{problem}).
#'   
#' @return Numeric vector containing the values of the best found solutions 
#'   during each run.
#'   
#' @export
getBestSolutionValues <- function(data, problem, search){
  UseMethod("getBestSolutionValues")
}
#' @export
getBestSolutionValues.james <- function(data, problem, search){
  # extract search runs
  runs <- getSearchRuns(data, problem, search)
  # get best solution values
  best.solution.values <- sapply(runs, function(run){ tail(run$values, n=1) })
  return(best.solution.values)
}

#' Get best found solutions
#' 
#' Get the best found solutions during the different runs of the given 
#' \code{search} applied to the given \code{problem}. This is a generic S3 
#' method.
#' 
#' If the \code{data} contains results for a single problem only, the argument 
#' \code{problem} can be omitted. Likewise, if for the considered \code{problem}
#' results are available for a single search only, the argument \code{search}
#' can be omitted.
#' 
#' When writing results obtained from the analysis tools in the 'JAMES' extensions
#' module to a JSON file, one should provide a JSON converter for the solution 
#' type of the analyzed problems if it is desired that the actual best found 
#' solutions are contained in the output file. Therefore, these solutions might 
#' not be available for all problems, searches or search runs. In case a best 
#' solution is missing for a search run, the corresponding entry in the returned
#' list will be set to \code{NA}. It is possible that a list of only \code{NA} 
#' is returned.
#' 
#' @param data data object containing the analysis results
#' @param problem name of the analyzed problem. Can be omitted if the 
#'   \code{data} contains results for a single problem only.
#' @param search name of the applied search. Can be omitted if the \code{data} 
#'   contains results for a single search only (for the considered 
#'   \code{problem}).
#'   
#' @return List containing the best found solutions during each run. May contain
#'   \code{NA} values.
#'   
#' @export
getBestSolutions <- function(data, problem, search){
  UseMethod("getBestSolutions")
}
#' @export
getBestSolutions.james <- function(data, problem, search){
  # extract search runs
  runs <- getSearchRuns(data, problem, search)
  # get best solutions
  best.solutions <- lapply(runs, function(run){ run$best.solution })
  # replace NULL with NA
  best.solutions[sapply(best.solutions, is.null)] <- NA
  return(best.solutions)
}

#' Get convergence times
#' 
#' Get the convergence times of the different runs of the given \code{search} 
#' applied to the given \code{problem} (in milliseconds). This is a generic S3 
#' method.
#' 
#' If the \code{data} contains results for a single problem only, the argument 
#' \code{problem} can be omitted. Likewise, if for the considered \code{problem}
#' results are available for a single search only, the argument \code{search} 
#' can be omitted.
#' 
#' The convergence time of a search run is defined as the time at which a 
#' certain value threshold is crossed. This threshold is computed from the given
#' convergence ratio \code{r} as follows: if values are being maximized, 
#' \code{thr = (1-r)*min + r*max}; else, \code{thr = (1-r)*max + r*min}, where 
#' \code{min} and \code{max} are the minimum and maximum observed value in the 
#' considered search run, respectively. In case of maximization, a search run is
#' said to have converged as soon as it reaches a value which is larger than or 
#' equal to the threshold \code{thr}. In case of minimization, convergence 
#' occurs when the values drop below the threshold. The convergence ratio 
#' \code{r} defaults to 0.99. If set to 1, a search run is said to have 
#' converged when the final best solution is found.
#' 
#' @param data data object containing the analysis results
#' @param problem name of the analyzed problem. Can be omitted if the 
#'   \code{data} contains results for a single problem only.
#' @param search name of the applied search. Can be omitted if the \code{data} 
#'   contains results for a single search only (for the considered 
#'   \code{problem}).
#' @param r convergence ratio. Defaults to 0.99. Numeric value in [0,1].
#'   
#' @return Numeric vector containing the convergence times of each run (in
#'   milliseconds). All convergence times are greater than or equal to -1.
#'   
#' @export
getConvergenceTimes <- function(data, problem, search, r = 0.99){
  UseMethod("getConvergenceTimes")
}
#' @export
getConvergenceTimes.james <- function(data, problem, search, r = 0.99){
  
  # check input
  if(!is.numeric(r) || r < 0 || r > 1){
    stop("'r' should be a \"numeric\" value in [0,1]")
  }
  
  # extract search runs
  runs <- getSearchRuns(data, problem, search)
  
  # check whether values are maximized or minimized
  minimizing <- isMinimizing(data, problem)
  
  # compute convergence times
  conv.times <- sapply(runs, function(run){
    max.value <- max(run$values)
    min.value <- min(run$values)
    if(minimizing){
      thr <- (1-r)*max.value + r*min.value
      conv.time <- min(run$times[run$values <= thr])
    } else {
      thr <- (1-r)*min.value + r*max.value
      conv.time <- min(run$times[run$values >= thr])
    }
    return(conv.time)
  })
  
  return(conv.times)
}

################################## SUMMARY ################################## 

#' @export
summary.james <- function(object, ...){
  
  results <- object
  
  # determine column widths
  col.widths <- c(8, 7, 5, 11, 8, 8, 8)
  problem.names <- getProblems(results)
  search.names <- unique(unlist(sapply(problem.names, function(p){getSearches(results, p)})))
  col.widths[1] <- max(col.widths[1], max(sapply(problem.names, nchar)))
  col.widths[2] <- max(col.widths[2], max(sapply(search.names, nchar)))
  format <- sprintf("%%-%ds  %%-%ds  %%%ds  %%%ds  %%%ds  %%%ds  %%%ds\n",
                    col.widths[1], col.widths[2], col.widths[3],
                    col.widths[4], col.widths[5], col.widths[6],
                    col.widths[7])
  
  # construct header lines
  lines <- lapply(col.widths, function(w){ paste(rep("-", w), collapse = "") })
  
  # print header
  printf(format, "Problem:", "Search:", "Runs:", "Mean value:", "St. dev:", "Median:", "IQR:")
  printf(format, lines[[1]], lines[[2]], lines[[3]], lines[[4]], lines[[5]], lines[[6]], lines[[7]])
  
  # update format for real value printing (last two columns)
  format <- sprintf("%%-%ds  %%-%ds  %%%ds  %%11.3g  %%8.3g  %%8.3g  %%8.3g\n",
                    col.widths[1], col.widths[2], col.widths[3])
  
  # print info of applied searches per problem
  for (problem in problem.names){
    search.names <- getSearches(results, problem)
    for (search in search.names){
      nruns <- getNumSearchRuns(results, problem, search)
      values <- getBestSolutionValues(results, problem, search)
      mean.value <- mean(values)
      sd.value <- sd(values)
      median.value <- median(values)
      iqr.value <- IQR(values)
      printf(format, problem, search, nruns, mean.value, sd.value, median.value, iqr.value)
    }
  }
  
}







