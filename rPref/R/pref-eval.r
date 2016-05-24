

#' Preference Selection
#' 
#' Evaluates a preference on a given data set, i.e., 
#' returns the maximal elements of a data set for a given preference order.
#' 
#' @param df A data frame or, for a grouped preference selection, a grouped data frame. See below for details.
#' @param pref The preference order constructed via \code{\link{complex_pref}} and \code{\link{base_pref}}. 
#'             All variables occurring in the definition of \code{pref} must be either columns of the data frame \code{df} 
#'             or variables/functions of the environment where \code{pref} was defined.
#' @param ... Additional (optional) parameters for top(-level)-k selections:
#'  \describe{
#'    \item{\code{top}}{A \code{top} value of k means that the k-best tuples of the data set are returned. 
#'     This may be non-deterministic, see below for details.}
#'     \item{\code{at_least}}{An \code{at_least} value of k returns the top-k tuples and additionally all tuples which are 
#'     not dominated by the worst tuple (i.e. the minima) of the Top-k set. 
#'     The number of tuples returned is greater or equal than
#'      \code{at_least}. In contrast to top-k, this is deterministic.}
#'     \item{\code{top_level}}{A \code{top_level} value of k returns all tuples from the k-best levels. See below for the definition of a level.}
#'     \item{\code{and_connected}}{Logical value, which is only relevant if more than one of the above \{\code{top}, \code{at_least}, \code{top_level}\} 
#'     values is given, otherwise it will be ignored. 
#'     Then \code{and_connected = TRUE} (which is the default) means that all top-conditions 
#'     must hold for the returned tuples: 
#'     Let \code{cond1} and \code{cond2} be top-conditions like \code{top=2} or \code{top_level=3}, then
#'     \code{psel([...], cond1, cond2)} is equivalent to the intersection of \code{psel([...], cond1)} and \code{psel([...], cond2)}. If we have
#'     \code{and_connected = FALSE}, these conditions are or-connected. 
#'     This corresponds to the union of \code{psel([...], cond1)} and \code{psel([...], cond2)}.}
#'     \item{\code{show_level}}{Logical value. If \code{TRUE}, a column \code{.level} 
#'     is added to the returned data frame, containing all level values. 
#'     If at least one of the \{\code{top}, \code{at_least}, \code{top_level}\} values are given,
#'     then \code{show_level} is \code{TRUE} by default for the \code{psel} function. 
#'     Otherwise, and for \code{psel.indices} in all cases, this option is \code{FALSE} by default.}
#' }
#' 
#' @details
#' 
#' The difference between the two variants of the preference selection is:
#' 
#' \itemize{
#' \item The \code{psel} function returns a subset of the data set which are the maxima according to the given preference. 
#' \item The function \code{psel.indices} returns just the row indices of the maxima 
#' (except top-k queries with \code{show_level = TRUE}, see top-k preference selection).
#' Hence \code{psel(df, pref)} is equivalent to \code{df[psel.indices(df, pref),]} for non-grouped data frames. 
#' }
#' 
#' @section Top-k Preference Selection:
#' 
#' For a given \code{top} value of k the k best elements and their level values are returned. The level values are determined as follows:
#' 
#' \itemize{
#' \item{All the maxima of a data set w.r.t. a preference have level 1.}
#' \item{The maxima of the remainder, i.e., the data set without the level 1 maxima, have level 2.}
#' \item{The n-th iteration of "Take the maxima from the remainder" leads to tuples of level n.}
#' }
#' 
#' By default, \code{psel.indices} does not return the level values. By setting \code{show_level = TRUE} this function
#' returns a data frame with the columns '.indices' and '.level'. 
#' Note that, if none of the top-k values \{\code{top}, \code{at_least}, \code{top_level}\} is set,
#' then all level values are equal to 1. 
#' 
#' By definition, a top-k preference selection is non-deterministic. 
#' A top-1 query of two equivalent tuples (equivalence according to \code{pref})
#' can return on both of these tuples. 
#' For example, a \code{top=1} preference selection on the tuples {(a=1, b=1), (a=1, b=2)}
#' w.r.t. \code{low(a)} preference can return either the 'b=1' or the 'b=2' tuple. 
#' 
#' On the contrary, a preference selection using \code{at_least} is deterministic by adding all tuples having the same level as the worst level 
#' of the corresponding top-k query. This means, the result is filled with all tuples being not worse than the top-k result. 
#' A preference selection with top-level-k returns all tuples having level k or better. 
#'  
#' If the \code{top} or \code{at_least} value is greater than the number of elements in \code{df} 
#' (i.e., \code{nrow(df)}), or \code{top_level} is greater than the highest level in \code{df},
#' then all elements of \code{df} will be returned without further warning.
#' 
#' @section Grouped Preference Selection:
#' 
#' Using \code{psel} it is also possible to perform a preference selection where the maxima are calculated for every group separately. 
#' The groups have to be created with \code{\link{group_by}} from the dplyr package. The preference selection preserves the grouping, i.e.,
#' the groups are restored after the preference selection.
#' 
#' For example, if the \code{summarize} function from dplyr is applied to
#' \code{psel(group_by(...), pref)}, the summarizing is done for the set of maxima of each group. 
#' This can be used to e.g., calculate the number of maxima in each group, see the examples below.
#' 
#' A \{\code{top}, \code{at_least}, \code{top_level}\} preference selection
#' is applied to each group separately.
#' A \code{top=k} selection returns the k best tuples for each group. 
#' Hence if there are 3 groups in \code{df}, each containing at least 2 elements, 
#' and we have \code{top = 2}, then 6 tuples will be returned.
#' 
#' @section Parallel Computation:
#' 
#' On multi-core machines the preference selection can be run in parallel using a divide-and-conquer approach. 
#' Depending on the data set, this may be faster than a single-threaded computation.
#' To activate parallel computation within rPref the following option has to be set:
#' 
#' \code{options(rPref.parallel = TRUE)}
#' 
#' If this option is not set, rPref will use single-threaded computation by default.
#' 
#' @seealso See \code{\link{complex_pref}} on how to construct a Skyline preference. 
#' See \code{\link{plot_front}} on how to plot the Pareto front of a Skyline.
#' 
#' 
#' @name psel
#' @importFrom dplyr is.grouped_df group_by_
#' @export
#' 
#' @examples
#' 
#' # Skyline and top-k/at-least skyline
#' psel(mtcars, low(mpg) * low(hp))
#' psel(mtcars, low(mpg) * low(hp), top = 5)
#' psel(mtcars, low(mpg) * low(hp), at_least = 5)
#' 
#' # visualize the skyline in a plot
#' sky1 <- psel(mtcars, high(mpg) * high(hp))
#' plot(mtcars$mpg, mtcars$hp)
#' points(sky1$mpg, sky1$hp, lwd=3)
#' 
#' # grouped preference with dplyr
#' library(dplyr)
#' psel(group_by(mtcars, cyl), low(mpg))
#' 
#' # return size of each maxima group
#' summarise(psel(group_by(mtcars, cyl), low(mpg)), n())
#' 
psel <- function(df, pref, ...) {  
  
  vars = list(...)
  
  # Check if top-k query
  is_top <- any(names(vars) %in% c('top', 'at_least', 'top_level'))
  
  # Store actual show_level value
  if (is_top) {
    # top value is set => show_level is TRUE by Default
    show_level <- is.null(vars$show_level) || isTRUE(vars$show_level)
  } else {
    # top value is not set => show_level is FALSE by Default
    show_level <- isTRUE(vars$show_level)
  }
  vars$show_level <- show_level
  
  # Call psel.indices with all parameters from ...
  tres <- psel.indices(df, pref, .dots = vars)
  
  # Extract indices
  if (!show_level) indices <- tres # psel.indices.top returned just indices when show_level = FALSE
  else indices <- tres[['.indices']]
  
  # Preference subset
  if (!dplyr::is.grouped_df(df)) res <- df[indices,,drop=FALSE] # Usual preference selection
  else res <- dplyr::group_by_(df[indices,,drop=FALSE], .dots = attr(df, 'vars')) # Grouped preference selection - regroup!
  
  # Add levels if show_level is true
  if (!show_level) return(res)
  else {
    res[['.level']] <- tres[['.level']]
    return(res)
  }
  
  # Do the preference selection
  if (!dplyr::is.grouped_df(df)) return(df[res,]) # Usual preference selection
  else return(dplyr::group_by_(df[res,], .dots = attr(df, 'vars'))) # Grouped preference selection - regroup!
}



#' @export
#' @rdname psel
#' @importFrom dplyr is.grouped_df group_by_
#' @importFrom RcppParallel defaultNumThreads
psel.indices <- function(df, pref, ...) {
  
  # Note that also data.table fulfills is.data.frame(.) = TRUE
  if (!is.data.frame(df))          stop("Preference selection needs a data.frame as first argument.")
  if (!is.actual.preference(pref)) stop("Preference selection needs a preference as second argument.")
    
  # ** Get all special arguments for top-k selections
  
  vars <- list(...)
  
  # Get all variables from '...' in psel
  if ('.dots' %in% names(vars)) vars <- vars[['.dots']] 
    
  is_top <- any(c('top', 'at_least', 'top_level') %in% names(vars))
  
  # Show-level operates indepndant of is_top!
  show_level    <- isTRUE(vars[['show_level']]) # FALSE by Default (for psel.indices)
  
  if (is_top) { 
    # Logical options
    and_connected <- isTRUE(vars[['and_connected']]) || is.null(vars[['and_connected']])  # TRUE by default
    
    # No defaults, but value -1 if null
    top       <- check.singleint.null2minus(vars[['top']],       'top')
    top_level <- check.singleint.null2minus(vars[['top_level']], 'top_level')
    at_least  <- check.singleint.null2minus(vars[['at_least']],  'at_least')
    
    # Check if all "-1" (were NULL before) - no valid top-k selection!
    if (all(c(top, show_level, at_least) == -1)) stop('The parameters top, top_level, at_least must not be all NULL!')
  }
  
  # ** Check for additional (wrong) arguments
  
  unused_names = setdiff(names(vars), c('top', 'at_least', 'top_level', 'and_connected', 'show_level', '.dots', ''))
  
  if (!is.null(unused_names) && length(unused_names) > 0)
    warning(paste0('The following arguments passed to psel are no preference selection parameters and will be ignored: ', paste(unused_names, collapse = ', ')))

  if ((length(vars) > 0 && is.null(names(vars))) || ('' %in% names(vars)))
    warning("Unnamed arguments were passed to '...' in psel. They will be ignored.")
  
  # ** get grouping 
  
  # (grouped data.tables have been removed from dplyr 0.4.3.9)
  is_grouped <- dplyr::is.grouped_df(df)
  if (is_grouped) {
    group_indices <- attr(df, 'indices')
    if (is.null(group_indices)) 
      stop("Could not find grouping indices in grouped data.frame!")
  }
  
  # ** Calculate score/serial pref
 
  # Precalculate score values for given preference, get_corevals must be called before serialize!
  scores <- pref$get_scorevals(1, df)$scores
  pref_serial <- pref$serialize()
  
  # ** Get options
  
  # Get alpha value, default is 1
  alpha <- getOption("rPref.scalagon.alpha", default = 1)
  
  # Use parallel computation? Default is FALSE!
  if (isTRUE(getOption("rPref.parallel", default = FALSE)))
    # Get number of threads vom RcppParallel for parallelization
    Npar <- RcppParallel::defaultNumThreads()
  else
    Npar <- 1
  
  # ** Finally do the (top-k) preference selection
  
  if (!is_top) {
    # Do the preference selection - not-top-k
    if (!is_grouped) { # Usual preference selection (not grouped)
      res <- pref_select_impl(scores, pref_serial, Npar, alpha) # non parallel for Npar=1
    } else { # Grouped preference selection
      res <- grouped_pref_sel_impl(group_indices, scores, pref_serial, Npar, alpha) 
    }
    
    if (!show_level) # just return indices
      # All C indices start at 0, and all R indices start at 1
      return(res + 1)
    else
      # Add level values for is_top = FALSE, i.e., all level values are 1
      return(data.frame(.indices = res + 1, .level = 1))
    
  } else {
    
    # Do the top-k preference selection
    if (!is_grouped) { # Usual preference selection (not grouped)
      res <- pref_select_top_impl(scores, pref_serial, Npar, alpha, 
                                  top, at_least, top_level, and_connected, show_level) 
    } else { # Grouped preference selection
      res <- grouped_pref_sel_top_impl(group_indices, scores, pref_serial, Npar, alpha,
                                       top, at_least, top_level, and_connected, show_level) 
    }
    
    # All C indices start at 0, and all R indices start at 1
    res[['.indices']] <- res[['.indices']] + 1
    
    if (!show_level) return(res[['.indices']]) # Return just indices
    return(res) # Return Data.Frame(".indices", ".level")
  }
}


# Helper
check.singleint.null2minus <- function(val, name) {
  if (is.null(val)) 
    return(-1)
  else if (!(is.numeric(val) && length(val) == 1 && round(val) == val && val >= 0))
    stop(paste0('Error in preference selection: Parameter ', name, ' must be a positive single integer value'), call. = FALSE)
  else
    return(val)
}

