################################################################################
# from GA package:
setClassUnion("numericOrchar", members = c("numeric", "character"))


################################################################################
#' An S4 class to return the results of using a GA to estimate a FSM with
#' \code{\link{evolve_model}}.
#'
#' @slot call Language from the call of the function \code{\link{evolve_model}}.
#' @slot actions Numeric vector with the number of actions.
#' @slot states Numeric vector with the number of states.
#' @slot GA S4 object created by ga() from the GA package.
#' @slot state_mat Numeric matrix with rows as states and columns as predictors.
#' @slot action_vec Numeric vector indicating what action to take for each
#'   state.
#' @slot predictive Numeric vector of length one with test data accuracy if test
#'   data was supplied; otherwise, a character vector with a message that the
#'   user should provide test data for better estimate of performance.
#' @slot varImp Numeric vector same length as number of columns of state matrix
#'   with relative importance scores for each predictor.
#' @slot timing Numeric vector length one with the total elapsed seconds it took
#'   \code{\link{evolve_model}} to execute.
#' @slot diagnostics Character vector length one, to be printed with base::cat().
#'
#' @export

#' @importClassesFrom GA ga
setClass("ga_fsm",
         slots = c(call = "language",
                   actions = "numeric",
                   states = "numeric",
                   GA = "ga", # from package "GA"
                   state_mat = "matrix",
                   action_vec = "numeric",
                   predictive = "numericOrchar",
                   varImp = "numeric",
                   timing = "numeric",
                   diagnostics = "character")
)

################################################################################
#' @describeIn ga_fsm An S4 method for printing a ga_fsm S4 object
#' @param x S4 ga_fsm object.
#'  @export

setMethod("print", "ga_fsm",
          function(x, ...) utils::str(x)
)

################################################################################
#' @describeIn ga_fsm An S4 method for showing a ga_fsm S4 object

setMethod("show", "ga_fsm",
          function(object) {
                  cat("An object of class \"ga_fsm\"\n")
                  cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
                  cat("Available slots:\n")
                  print(methods::slotNames(object))
          }
)

################################################################################
#' Turns ga_fsm S4 object into list of summaries for printing and then prints it.
#' @describeIn ga_fsm An S4 method for summarizing a ga_fsm S4 object
#'
#' @param object S4 ga_fsm object
#' @param digits Optional numeric vector length one for how many significant digits to
#' print, default is 3.
#'
#'  @export

setMethod("summary", "ga_fsm",
          function(object, digits = 3) {
                  x <- list(
                          # ga-related
                          popSize = object@GA@popSize,
                          maxiter = object@GA@maxiter,
                          elitism = object@GA@elitism,
                          pcrossover = object@GA@pcrossover,
                          pmutation = object@GA@pmutation,
                          iter = object@GA@iter,
                          fitness = object@GA@fitnessValue,
                          bit_string_solution = object@GA@solution,
                          # fsm-related
                          actions = object@actions,
                          states = object@states,
                          state_mat = object@state_mat,
                          action_vec = object@action_vec,
                          predictive = object@predictive,
                          varImp = object@varImp
                  )
                  cat("                                    \n")
                  cat("Gentic Algorithm Settings: \n")
                  cat(paste("Population size       = ", x$popSize, "\n"))
                  cat(paste("Number of generations = ", x$maxiter, "\n"))
                  cat(paste("Elitism               = ", x$elitism, "\n"))
                  cat(paste("Crossover probability = ", format(x$pcrossover, digits = digits), "\n"))
                  cat(paste("Mutation probability  = ", format(x$pmutation, digits = digits), "\n"))

                  cat("\nFinite State Machine Settings: \n")
                  cat(paste("Actions = ", x$actions, "\n"))
                  cat(paste("States  = ", x$states, "\n"))

                  cat("\nResults: \n\n")
                  cat(paste("Iterations For This Run              =", format(x$iter, digits = digits), "\n"))
                  cat(paste("Training Data Fitness Function Value =", format(x$fitness, digits = digits), "\n"))
                  cat(paste("Test Data Fitness Function Value     =", format(x$predictive, digits = digits), "\n"))

                  cat(paste("\n(Bit String Form) of Solution: \n"))
                  print(x$bit_string_solution[1, ], digits = digits)

                  cat("\nState Matrix of Solution: \n")
                  print(x$state_mat, digits = digits)

                  cat("\nAction Vector of Solution: \n")
                  print(x$action_vec, digits = digits)

                  cat("\nVariable Importance: \n")
                  print(x$varImp, digits = digits)

                  invisible(x)
          }
)

################################################################################
#' Plots ga_fsm S4 object's state transition matrix
#' @describeIn ga_fsm
#'
#' @aliases plot,ga_fsm-method
#' 
#' @param y not used.
#' @param maintitle optional character vector 
#' @param action_label optional character vector same length as action vector,
#'   where each ith element corresponds to what that ith element in the action
#'   vector represents. This will be used to fill in the states (circles) of the
#'   state transition matrix to be plotted.
#' @param transition_label optional character vector same length as number of
#'   columns of state transition matrix.
#' @param curvature optional numeric vector specifying the curvature of the
#'   lines for a diagram of 2 or more states.
#'   
#' @export

setMethod("plot", "ga_fsm",
          function(x, y, maintitle = "Transition Diagram",
                   action_label = NULL, 
                   transition_label = NULL,
                   curvature = c(0.3, 0.6, 0.8)) {
            actions <- x@actions
            states <- x@states
            state_mat <- x@state_mat
            action_vec <- x@action_vec
            predictive <- x@predictive
            varImp <- x@varImp
            
            s <- t(state_mat)
            
            if (missing(action_label))
              action_label <- paste0(action_vec)
            
            if (missing(transition_label))
              transition_label <- as.character(seq(nrow(s)))
            
            M <- as.data.frame(matrix(nrow = length(action_vec), ncol = length(action_vec), 
                                      byrow = TRUE, data = 0))
            
            for (i in seq(nrow(s))){
              for (j in seq(ncol(s))){
                if (M[s[i, j], j] == "0"){
                  M[s[i, j], j] <- transition_label[i] #row.names(s)[i]
                } else {
                  M[s[i, j], j] <- as.character(paste(M[s[i, j], j], transition_label[i], sep=" - ")) #paste0(M[s[i, j], j], "\n\n", row.names(s)[i])
                }
              }
            }
            
            diagram::plotmat(M, 
                             pos = length(action_vec), 
                             curve = curvature[(length(action_vec)-1)], 
                             name = action_label, 
                             lwd = 1, box.lwd = 2, 
                             cex.txt = 0.8, 
                             box.type = "circle", 
                             box.col = "lightblue",
                             box.prop = 1,
                             add = FALSE,
                             main  = maintitle)
          }
)

################################################################################
#' Plots ga_fsm S4 object's variable importances
#' @describeIn ga_fsm
#'
#' @param height ga_fsm S4 object
#' @param ... arguments to be passed to/from other methods.
#'   
#' @export

setMethod("barplot", "ga_fsm",
          function(height, ...) {
            height <- varImp(height)
            barplot(height, ...)
          }
)

################################################################################
#' Plots ga_fsm S4 object's variable importances
#' @describeIn ga_fsm
#' @param labels  vector of labels for each point. For vectors the default is to
#'   use names(x) and for matrices the row labels dimnames(x)[[1]].
#' @export

setMethod("dotchart", "ga_fsm",
          function(x, labels) {
            x <- varImp(x)
            if (missing(labels)) return(dotchart(x))
            if (!missing(labels)) return(dotchart(x, labels))
          }
)

################################################################################
#' Extracts slot relevant to estimating the fsm
#' @param x S4 ga_fsm object
#' @export

setGeneric("estimation_details", function(x){
  standardGeneric("estimation_details")
})

################################################################################
#' Extracts slot relevant to estimating the fsm
#' @describeIn ga_fsm
#'  @export

setMethod("estimation_details", "ga_fsm",
          function(x) methods::slot(x, "GA")
)

################################################################################
#' Extracts performance
#' @param x S4 ga_fsm object
#' @export

setGeneric("best_performance", function(x){
  standardGeneric("best_performance")
})

################################################################################
#' Extracts performance
#' @describeIn ga_fsm
#'  @export

setMethod("best_performance", "ga_fsm",
          function(x) estimation_details(x)@fitnessValue
)

################################################################################
#' Extracts slot of variable importances
#' @param x S4 ga_fsm object
#' @export

setGeneric(name = "varImp", function(x){
  standardGeneric("varImp")
})

################################################################################
#' Extracts slot of variable importances
#' @describeIn ga_fsm
#'  @export

setMethod("varImp", "ga_fsm",
          function(x) methods::slot(x, "varImp")
)

################################################################################
#' Extracts slot of action_vec
#' @param x S4 ga_fsm object
#' @export

setGeneric(name = "action_vec", function(x){
  standardGeneric("action_vec")
})

################################################################################
#' Extracts slot of action_vec
#' @describeIn ga_fsm
#'  @export

setMethod("action_vec", "ga_fsm",
          function(x) methods::slot(x, "action_vec")
)

################################################################################
#' Extracts number of states
#' @param x S4 ga_fsm object
#' @export

setGeneric(name = "states", function(x){
  standardGeneric("states")
})

################################################################################
#' Extracts number of states
#' @describeIn ga_fsm
#'  @export

setMethod("states", "ga_fsm",
          function(x) methods::slot(x, "states")
)

################################################################################
#' @describeIn ga_fsm
#' @param type Not currently used.
#' @param na.action Optional function.
#' @inheritParams evolve_model
#' 
#' @export
setMethod("predict", "ga_fsm",
          function(object, data,
                   type = "prob", na.action = stats::na.omit, ...){
            
            ## Data-related errors:
            if (missing(data)) 
              stop(paste("You must supply data. At the very least, you can supply data.",
                         "This should be a data.frame that has columns named 'period' and 'outcome' (period",
                         "is the time period that the outcome action was taken), and the rest of the",
                         "columns are predictors, ranging from one to three predictors. All of the",
                         "(3-5 columns) should be named. The period and outcome columns should be",
                         "integer vectors and the columns with the predictor variable data should be",
                         "logical vectors."))
            if(!is.data.frame(data)) {
              data <- as.data.frame(data)
              warning(paste("You did not supply a data.frame for 'data' argument of this function",
                            "so we converted it to one. To ensure this works right, run this again with a",
                            "data.frame that has columns named 'period' and 'outcome' (period",
                            "is the time period that the outcome action was taken), and the rest of the",
                            "columns are predictors, ranging from one to three predictors. All of the",
                            "(3-5 columns) should be named. The period and outcome columns should be",
                            "integer vectors and the columns with the predictor variable data should be",
                            "logical vectors."))
            }
            
            ## Data preparation:
            period <- data$period
            
            inputs <- 2^(ncol(data[ , -which(names(data) %in% c("period", "outcome")), drop = FALSE]))
            
            # change any non-logical predictor variable vectors to logical
            data[ , -which(names(data) %in% c("period", "outcome"))] <-
              data.frame(lapply(data[ , -which(names(data) %in% c("period", "outcome"))],
                                function(x) {
                                  if (class(x)!="logical") {
                                    as.logical(x)
                                  } else {
                                    x
                                  }}))
            # replace all NA's with 0 or 1 so these rows are not dropped
            # this works fine if the NAs are only for the first period play bc
            # then the predictor columns dont make a difference bc the FSM will initialize
            # with the same action regardless of the predictors at that time
            # but this would bias the results if NA's are occuring in predictors in other periods
            # so return an error for that:
            if (any(!stats::complete.cases(data) & !data$period == 1))
              stop(paste("Error: You have missing values in your training data somewhere other than the first period interactions.",
                         "You can only have missing values for predictor columns, AND these must be in rows where period==1."))
            data[is.na(data)] <- TRUE
            
            names <- colnames(data[ , -which(names(data) %in% c("period", "outcome")), drop = FALSE])
            
            if (length(names)==1){
              form <- paste("outcome ~ 0 +", names, sep=" ")
              data <- stats::model.matrix(eval(parse(text=form)), data)
            } else {
              predictors <- paste(names, collapse=":")
              form <- paste("outcome ~ 0 +", predictors, sep=" ")
              data <- stats::model.matrix(eval(parse(text=form)), data)
            }
            
            if (ncol(data) != inputs)
              stop(paste("Error: At least one of your predictor variables in your data",
                         "does not have exactly 2 levels."))
            
            ##########
            state_mat <- object@state_mat
            action_vec <-  object@action_vec
            results <- fitnessCPP(action_vec, state_mat, data, period)
            
            ##########
            if (anyNA(results) | length(results)==0){
              warning("Error: Results from fitness evaluation have missing values.")
            }
            
            results
          }
)
