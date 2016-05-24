
############################## CONVERGENCE CURVES ##############################

#'Plot convergence curves
#'
#'Creates a plot showing the convergence curve of each search that has been 
#'applied to the given \code{problem}, aggregated over all search runs (mean or 
#'median). This is a generic S3 method.
#'
#'If the \code{data} contains results for a single problem only, the argument 
#'\code{problem} can be omitted. If desired to plot convergence curves for a 
#'selection of the applied searches, use \code{\link{reduceJAMES}} to extract 
#'the respective data.
#'
#'The curves are plotted using \code{\link{matplot}}. More information about the
#'graphical parameters are provided in the documentation of this function. By 
#'default, a legend is added to the plot. This can be omitted by setting 
#'\code{legend = FALSE}. If desired, a custom legend may then be added. It is 
#'possible to zoom in on a specific region of the plot using the parameters 
#'\code{min.time} and \code{max.time}.
#'
#'Any additional parameters are passed to \code{\link{matplot}}.
#'
#'@seealso \code{\link{matplot}}
#'  
#'@param data data object containing the analysis results
#'@param problem name of the problem for which the plot is made. Can be omitted 
#'  if the \code{data} contains results for a single problem only.
#'@param type one of \code{"mean"} (default) or \code{"median"}. Determines how 
#'  the values from the different search runs are aggregated.
#'@param col color(s) of the plotted lines and/or symbols, used cyclically when 
#'  providing a vector. Defaults to \code{"black"}.
#'@param plot.type defaults to \code{"s"} (staircase). See \code{\link{matplot}}
#'  and \code{\link{plot}} for more information about the possible plot types.
#'@param lty line type(s), used cyclically when providing a vector. Line types 
#'  default to 1:n where n is the number of plotted curves.
#'@param title plot title. Defaults to \code{"Convergence curve(s)"}.
#'@param subtitle plot subtitle. By default, a subtitle will be added that 
#'  states the name of the problem for which the plot was made. If no substitle 
#'  is desired, set \code{subtitle = ""}.
#'@param xlab x-axis label. Defaults to \code{"Time"} followed by the time unit 
#'  within brackets (abbreviated).
#'@param ylab y-axis label. Defaults to \code{"Value"}.
#'@param time.unit one of \code{"milliseconds"} (default), \code{"seconds"}, 
#'  \code{"minutes"} or \code{"hours"}. Determines the time unit of the values 
#'  on the \code{x-}axis.
#'@param min.time zoom in on the part of the curve(s) after this point in time 
#'  on the \code{x-}axis, according to the specified \code{time.unit}.
#'@param max.time zoom in on the part of the curve(s) before this point in time
#'  on the \code{x-}axis, according to the specified \code{time.unit}.
#'@param legend logical: indicates whether a legend should be added to the plot.
#'  Defaults to \code{TRUE}.
#'@param legend.pos position of the legend, specified as a keyword  keyword from
#'  the list \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}, 
#'  \code{"left"}, \code{"topleft"}, \code{"top"}, \code{"topright"}, 
#'  \code{"right"} and \code{"center"}. Defaults to \code{"bottomright"} in case
#'  values are being maximized or \code{"topright"} in case of minimization.
#'@param legend.inset inset distance(s) from the margins as a fraction of the 
#'  plot region. If a single value is given, it is used for both margins; if two
#'  values are given, the first is used for \code{x-}distance, the second for 
#'  \code{y-}distance. Defaults to \code{c(0.02, 0.05)}.
#'@param legend.names names to be shown in the legend. Defaults to the search 
#'  names obtained from calling \code{\link{getSearches}} for the given 
#'  \code{data} and \code{problem}.
#'@param ... optional other arguments passed to \code{\link{matplot}}.
#'  
#'@export
plotConvergence <- function(data, problem, type = c("mean", "median"),
                            col = "black", plot.type = "s", lty,
                            title = "Convergence curve(s)", subtitle,
                            xlab, ylab = "Value",
                            time.unit = c("milliseconds", "seconds", "minutes", "hours"),
                            min.time, max.time, legend = TRUE, legend.pos,
                            legend.inset = c(0.02, 0.05), legend.names,
                            ...){
  UseMethod("plotConvergence")
}
#' @export
plotConvergence.james <- function(data, problem, type = c("mean", "median"),
                                  col = "black", plot.type = "s", lty,
                                  title = "Convergence curve(s)", subtitle,
                                  xlab, ylab = "Value",
                                  time.unit = c("milliseconds", "seconds", "minutes", "hours"),
                                  min.time, max.time, legend = TRUE, legend.pos,
                                  legend.inset = c(0.02, 0.05), legend.names,
                                  ...){
  
  # fall back to single problem if not specified
  if(missing(problem)){
    problem <- getSingleProblem(data)
  }
  
  # check input
  agg.function <- get(match.arg(type))
  
  time.unit <- match.arg(time.unit)
  
  if(!missing(min.time) && !is.numeric(min.time)){
    stop("'min.time' should be \"numeric\"")
  }
  if(!missing(max.time) && !is.numeric(max.time)){
    stop("'max.time' should be \"numeric\"")
  }
  if(!missing(min.time) && !missing(max.time) && min.time >= max.time){
    stop("'min.time' should be smaller than 'max.time'")
  }
  
  if(!is.logical(legend)){
    stop("'legend' should be \"logical\"")
  }
  
  # get names of applied searches (sorted in natural order)
  searches <- getSearches(data, problem)
  num.searches <- length(searches)
  # aggregate runs for each search
  agg.curves <- list()
  for(search in searches){
    # extract search runs
    runs <- getSearchRuns(data, problem, search)
    # reduce runs to last observed value at same time
    runs <- lapply(runs, function(run){
      reduced.times <- unique(run$times)
      reduced.values <- sapply(reduced.times, function(t){
        tail(run$values[run$times == t], 1)
      })
      run$times <- reduced.times
      run$values <- reduced.values
      return(run)
    })
    # aggregate runs
    n.runs <- length(runs)
    runs.cur.indices <- rep(1, n.runs)
    runs.last.indices <- sapply(runs, function(run){length(run$times)})
    runs.cur.values <- sapply(runs, function(run){run$values[1]})
    agg.times <- c()
    agg.values <- c()
    while(max(runs.last.indices - runs.cur.indices) >= 0){
      # calculate first point in time when a new value is reached in any run
      next.times <- sapply(1:n.runs, function(r){
        run <- runs[[r]]
        cur.index <- runs.cur.indices[r]
        return(run$times[cur.index])
      })
      min.t <- min(next.times, na.rm = TRUE)
      # gather indices of runs with an update
      updated.runs.indices <- sapply(1:n.runs, function(r){
        run <- runs[[r]]
        cur.index <- runs.cur.indices[r]
        t <- run$times[cur.index]
        return(ifelse(!is.na(t) && t == min.t, r, NA))
      })
      updated.runs.indices <- updated.runs.indices[!is.na(updated.runs.indices)]
      # extract new values for updated runs
      updated.runs.values <- sapply(updated.runs.indices, function(r){
        run <- runs[[r]]
        cur.index <- runs.cur.indices[r]
        return(run$values[cur.index])
      })
      # update values
      runs.cur.values[updated.runs.indices] <- updated.runs.values
      # increase current indices of updated runs
      runs.cur.indices[updated.runs.indices] <- runs.cur.indices[updated.runs.indices] + 1
      # register aggregated value and corresponding time
      agg.times <- c(agg.times, min.t)
      agg.values <- c(agg.values, agg.function(runs.cur.values))
    }
    
    # rescale times to requested time unit
    agg.times <- agg.times / getTimeUnitMillis(time.unit)
    
    # zoom in on specific time interval, if desired
    if(!missing(min.time)){
      # save values before minimum time
      pre <- agg.values[agg.times < min.time]
      # truncate
      agg.values <- agg.values[agg.times >= min.time]
      agg.times <- agg.times[agg.times >= min.time]
      # prepend start value, if available
      if(length(pre) > 0){
        agg.values <- c(tail(pre, 1), agg.values)
        agg.times <- c(min.time, agg.times)
      }
    }
    if(!missing(max.time)){
      # save values after maximum time
      post <- agg.values[agg.times > max.time]
      # truncate
      agg.values <- agg.values[agg.times <= max.time]
      agg.times <- agg.times[agg.times <= max.time]
      # append final value, if available
      if(length(post) > 0){
        agg.values <- c(agg.values, head(post, 1))
        agg.times <- c(agg.times, max.time)
      }
    }
    
    # register aggregated curve
    avg.curve <- list(times = agg.times, values = agg.values)
    agg.curves <- append(agg.curves, list(avg.curve))
    
  }
  
  # extract maximum number of data points and final time point (across all curves)
  max.num.points <- max(sapply(agg.curves, function(curve){length(curve$times)}))
  final.time.point <- max(sapply(agg.curves, function(curve){tail(curve$times, 1)}))
  # initialize matrix of time vectors (columns)
  times.matrix <- matrix(NA, nrow = max.num.points, ncol = num.searches)
  # initialize matrix of value vectors (columns)
  values.matrix <- matrix(NA, nrow = max.num.points, ncol = num.searches)
  
  # fill both matrices
  for(i in 1:num.searches){
    curve <- agg.curves[[i]]
    num.points <- length(curve$times)
    times.matrix[1:num.points, i] <- curve$times
    values.matrix[1:num.points, i] <- curve$values
    # repeate final value until final time point, if necessary
    if(num.points < max.num.points){
      times.matrix[num.points+1, i] <- final.time.point
      values.matrix[num.points+1, i] <- values.matrix[num.points, i]
    }
  }

  # set default line types
  if(missing(lty)){
    lty <- 1:num.searches
  }
  # set default subtitle
  if(missing(subtitle)){
    subtitle <- sprintf("Problem: %s", problem)
  }
  # set default x-axis label
  if(missing(xlab)){
    xlab <- sprintf("Time (%s)", getTimeUnitAbbreviation(time.unit))
  }
  # plot curves
  matplot(x = times.matrix, y = values.matrix,
          col = col, type = plot.type, lty = lty,
          main = title, xlab = xlab, ylab = ylab,
          ...)
  # add subtitle
  mtext(subtitle, line = 0.25)
  # add legend, if desired
  if(legend){
    # set default position if not specified
    if(missing(legend.pos)){
      legend.pos <- ifelse(isMinimizing(data, problem), "topright", "bottomright")
    }
    # set default legend names if not specified
    if(missing(legend.names)){
      legend.names = searches
    }
    # add legend
    legend(legend.pos, legend = legend.names,
           inset = legend.inset, col = col, lty = lty)
  }
  
}

################################## BOX PLOTS ##################################

#' Solution quality and convergence time box plots
#' 
#' Produce box-and-whisker plots for the searches that have been applied to the 
#' given problem, visualizing the distribution of the best found solution's 
#' value (solution quality) or time until convergence in subsequent search runs.
#' 
#' If the data \code{x} contains results for a single problem only, the argument
#' \code{problem} can be omitted. If desired to produce box plots for a 
#' selection of the applied searches, use \code{\link{reduceJAMES}} to extract 
#' the respective data.
#' 
#' Convergence times are computed with \code{\link{getConvergenceTimes}}.
#' 
#' The plots are made using the generic \code{\link{boxplot}} method called on a
#' list of vectors containing the distribution samples for each search.
#' 
#' Any additional parameters are passed to \code{\link{boxplot}}.
#' 
#' @seealso \code{\link{getConvergenceTimes}}, \code{\link{boxplot}}
#'   
#' @param x data object containing the analysis results
#' @param problem name of the problem for which the plot is made. Can be omitted
#'   if the data \code{x} contains results for a single problem only.
#' @param type one of \code{"quality"} (default) or \code{"time"}. If set to 
#'   \code{"quality"}, the final solution's value is reported; if set to 
#'   \code{"time"}, the time until convergence is reported. In both cases, the 
#'   respective distribution of values found during the different search runs is
#'   visualized. In the latter case, the argument \code{r} is used to decide 
#'   when a search run has converged.
#' @param r convergence ratio, only used if \code{type} is \code{"time"}. 
#'   Defaults to 0.99. Should be a numeric value in [0,1]. This parameter is 
#'   passed to \code{\link{getConvergenceTimes}}.
#' @param time.unit one of \code{"milliseconds"} (default), \code{"seconds"}, 
#'   \code{"minutes"} or \code{"hours"}. Only used if \code{type} is 
#'   \code{"time"}. Determines the time unit of the convergence times on the
#'   \code{y-}axis.
#' @param title plot title. Defaults to \code{"Solution quality"} or 
#'   \code{"Convergence time"} when \code{type} is set to \code{"quality"} or 
#'   \code{"time"}, respectively.
#' @param subtitle plot subtitle. By default, a subtitle is added that states 
#'   the name of the problem for which the plot is made. If \code{type} is 
#'   \code{"time"} the subtitle also mentions the applied convergence ratio 
#'   \code{r}. If no subtitle is desired set \code{subtitle = ""}.
#' @param ylab y-axis label. Defaults to \code{"Value"} or \code{"Time"} (with 
#'   the time unit indicated between brackets) when \code{type} is set to 
#'   \code{"quality"} or \code{"time"}, respectively.
#' @param names names to be shown on the x-axis under the box plots. Defaults to
#'   the search names obtained from calling \code{\link{getSearches}} for the 
#'   given data \code{x} and \code{problem}.
#' @param ... any additional parameters are passed to \code{\link{boxplot}}.
#'   
#' @importFrom graphics boxplot
#' @export
boxplot.james <- function(x, problem, type = c("quality", "time"), r = 0.99,
                          time.unit = c("milliseconds", "seconds", "minutes", "hours"),
                          title, subtitle, ylab, names, ...){
  
  # rename data object
  data <- x
  
  # check input
  type <- match.arg(type)
  time.unit <- match.arg(time.unit)
  
  # fall back to single problem if missing
  if(missing(problem)){
    problem <- getSingleProblem(data)
  }
  
  # get searches
  searches <- getSearches(data, problem)
  
  # set function to extract values from search runs
  extract <- switch(type,
                    
    # extract final values
    "quality" = getBestSolutionValues,
    
    # extract convergence times (using specified convergence ratio)
    "time" = function(...){getConvergenceTimes(..., r = r) / getTimeUnitMillis(time.unit)}
  )
  
  # extract values from runs of each search
  search.dists <- lapply(searches, function(search){
    run.values <- extract(data, problem, search)
    return(run.values)
  })
  
  # set default title, subtitle and ylab if missing
  if(type == "quality"){
    if(missing(title)){
      title <- "Solution quality"
    }
    if(missing(subtitle)){
      subtitle <- sprintf("Problem: %s", problem)
    }
    if(missing(ylab)){
      ylab <- "Value"
    }
  } else {
    if(missing(title)){
      title <- "Convergence time"
    }
    if(missing(subtitle)){
      subtitle <- sprintf("Problem: %s (convergence ratio: %.16g)", problem, r)
    }
    if(missing(ylab)){
      ylab <- sprintf("Time (%s)", getTimeUnitAbbreviation(time.unit))
    }
  }
  # set default names if missing
  if(missing(names)){
    names <- searches
  }
  
  # create box pot
  boxplot(search.dists, main = title, ylab = ylab, names = names, ...)
  # add subtitle
  mtext(subtitle, line = 0.25)
  
}






