globalVariables(".")

#' Append an object to an environment
#' 
#' A function to coerce an object to a list and append the list to an environment
#' 
#' @importFrom magrittr %>%
#' @export
#' 
#' @param object An object which can be coerced to a list (e.g. an environment)
#' @param environment An environment
#' 
#' @return An appended environment
#' 
#' @examples
#' toEnv(data.frame(a = 1, b = 2), environment())
#' toEnv(list(a = 1, b = 2), environment())
#' toEnv(environment(), new.env())
toEnv = function(object, environment)
  object %>%
  as.list %>%
  list2env(environment)

#' Append an object to a dataframe.
#' 
#' Convert an object to a list, select only vector entries, coerce to a data.frame, 
#' and append to the given data frame.
#' 
#' @importFrom magrittr %>%
#' @export
#' 
#' @param object An object which can be coerced to a list (e.g. an environment)
#' @param dataframe A data frame
#' 
#' @return An appended dataframe
#' 
#' @examples
#' toDf(list(a = 1, b = 2), data.frame())
#' toDf(environment(), data.frame())
toDf = function(object, dataframe)
  append =
    object %>%
    as.list %>%
    rlist::list.filter(is.vector(.)) %>%
    as.data.frame %>%
    dplyr::bind_rows(dataframe, .)
  
  #' An implementation of a SAS datastep in a class
  #'  
  #' @docType class
  #' @format An \code{\link{R6Class}} generator object
  #' 
  #' @importFrom R6 R6Class
  #' @importFrom magrittr %>%
  #' @importFrom lazyeval lazy
  #' @export
  #' 
  #' @section Fields:
  #' 
  #' \describe{\item{\code{i}}{
  #'     \code{i} begins at 0 and is incremented for each iteration of the data step.}}
  #'   
  #' \describe{\item{\code{results}}{
  #'    The \code{results} frame is initialized as an empty data frame. 
  #' It is populated row-wise with each iteration of the data step.}}
  #'   
  #' \describe{\item{\code{continue}}{
  #'     \code{continue} is a marker which signals that the step should 
  #' continue repeating. When \code{continue} is 1, repetition will continue, and when 
  #' \code{continue} is 0, repitition will cease. It is initialized to 0.}}
  #'   
  #' \describe{\item{\code{eval}}{
  #'     \code{eval} is initialized as NULL, but will store a pointer to 
  #' the current evaluation environment. This pointer is helps pass the evaluation 
  #' environment from one iteration of the data step to the next.}}
  #'
  #' @section Methods:
  #' 
  #' \describe{\item{\code{begin(env)}}{
  #'     \code{begin} does three things: imports the environment of the previous 
  #' step to the current, stores the current environment (or the environment specified), 
  #' and increments \code{i} by 1. It takes one argument, \code{envir}, which should 
  #' typically be set to \code{environment()}.}}
  #' 
  #' \describe{\item{\code{set(dataframe, group_id)}}{
  #'     \code{set} takes two arguments: a data frame and an onptional unquoted \code{group_id} 
  #' variable. This \code{group_id} variable must contain a consecutive sequence of natural 
  #' numbers from 1 to some maximum. In each data step, rows where \code{i} matches the 
  #' \code{group_id} variable (or simply the ith row if no group_id variable is given) are selected, 
  #' and the slice is split into vectors and imported into the evaluation environment. 
  #' \code{continue} is set to 0 once \code{set} reaches the maximum value in the \code{group_id} 
  #' column, ceasing repetition of the datastep, else \code{continue} is set to 1.}}
  #' 
  #' \describe{\item{\code{set_(dataframe, group_id)}}{
  #'     A standard evaluation version of \code{set_}, in which the \code{group_id} 
  #' variable is included as a string, formula, or lazy object.}}
  #' 
  #' \describe{\item{\code{output}}{
  #'     \code{output} takes an optional list argument. Either the list, or, if none is given,
  #' all vectors in the evaluation environment are gathered into a data.frame, and this data.frame 
  #' appended to \code{results}.}}
  #' 
  #' \describe{\item{\code{end}}{
  #'     \code{end} will, if \code{continue} is 1, evaluate 
  #' the function given within the evaluation environment. Typically, the function 
  #' given will be the current function: that is, steps are joined recursively.}}
  #' 
  #' @examples
  #' step = dataStepClass$new()
  #' 
  #' frame = data.frame(x = 1:10)
  #' 
  #' stairs = function() {
  #'   step$begin(environment())
  #'   step$set(frame)
  #'   y = x + 1
  #'   step$output()
  #'   step$end(stairs)
  #' }
  #' 
  #' stairs()
  #' 
  #' step$results
  dataStepClass = R6::R6Class(
    "dataStepClass",
    public = list(
      i = 0,
      continue = 0,
      results = data.frame(),
      eval = NULL,
      
      begin = function(env) {
        # transfer in old environment
        self$eval %>% toEnv(env)
        # evaluate within new environment
        self$eval = env
        # iterate
        self$i = self$i + 1
      },
      
      set_ = function(dataframe, group_id) {
        
        if (missing("group_id")) {
          slice = dataframe %>% dplyr::slice(self$i)
          max = nrow(dataframe) 
        } else {
          #do some lazy interpretation
          filter = 
            lazyeval::lazy(group_id == self$i) %>%
            lazyeval::interp(group_id = group_id)
          
          #add dataframe slice to the environment
          slice =
            dataframe %>%
            dplyr::filter_(filter)
          
          max = 
            lazyeval::lazy(max(dataframe$group_id)) %>%
            lazyeval::interp(group_id = group_id) %>%
            lazyeval::lazy_eval()
        }
        
        slice %>% toEnv(self$eval)
        
        # if we've reached the end, mark for the end of recursion
        if (self$i != max) self$continue = 1 else self$continue = 0
      },
      
      set = function(dataframe, group_id) {
        if (missing("group_id")) self$set_(dataframe) else
          self$set_(
            dataframe, 
            lazyeval::lazy(group_id))
      },
      
      output = function(list = self$eval) {
        # add eval contents to results dataframe
        self$results = list %>% toDf(self$results)
      },
      
      end = function(FUN) {
        if (self$continue == 1) 
          do.call(FUN, args = list(), envir = self$eval)
        # if marked for recursion, recur within the eval environment
      }
    )
  )