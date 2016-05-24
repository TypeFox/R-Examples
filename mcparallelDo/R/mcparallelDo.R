#'A repository for a variety of useful functions.
#'
#' The primary function of this package is mcparallelDo().  
#' To use mcparallelDo(), simply invoke the function with a curly braced wrapped code and the character element name to which you want to assign the results.
#'
#' @name mcparallelDo-package
#' @docType package
#' @title mcparallelDo-package placeholder
#'
NULL
.onAttach <- function(libname,pkgname) {
  if (.Platform$OS.type != "unix") {
    warning("'mcparallelDo' only performs parallel processing on unix alikes; there will be no further warnings")
  }
}
NULL

#' The mcparallelDoManager Class and Object
#' @aliases mcparallelDoManager
#' @docType class
#' @importFrom R6 R6Class
mcparallelDoManagerClass <- R6::R6Class("mcparallelDoManager",
        public = list(
          h = taskCallbackManager()
          ,runningJobs = list()
          ,addJob = function(jobName, targetValue, verbose, targetEnvironment) {
            self$h$add(jobCompleteSelfDestructingHandler(jobName, targetValue, verbose, targetEnvironment))
            self$runningJobs[[jobName]] <- list(jobName=jobName, targetValue=targetValue, verbose=verbose, targetEnvironment=targetEnvironment)
            assign(targetValue, value = NULL, envir = targetEnvironment)
          }
          ,removeJob = function(x) {
            self$runningJobs <- self$runningJobs[names(self$runningJobs)!=x]
          }
          ,checkJobs = function() {
            sapply(names(self$runningJobs), function(x) {
              checkIfJobStillRunning(
                targetJob = x, 
                targetValue = self$runningJobs[[x]]$targetValue, 
                verbose = self$runningJobs[[x]]$verbose, 
                targetEnvironment = self$runningJobs[[x]]$targetEnvironment
              )
            })
          }
        )
)
.mcparallelDoManager <- mcparallelDoManagerClass$new()

#' mcparallelDoCheck
#'
#' Forces a check on all mcparallelDo() jobs and returns their values to the target environment if they are complete.
#' @return A named logical vector, TRUE if complete, FALSE if not complete, and an empty logical vector if on Windows
#' @export
mcparallelDoCheck <- function() {
  # Special handling for Windows
  if (.Platform$OS.type != "unix") {
    return(logical())
  }
  
  jobNames <- names(.mcparallelDoManager$runningJobs)
  jobStatus <- !.mcparallelDoManager$checkJobs()
  names(jobStatus) <- jobNames
  return(jobStatus)
}
NULL

#' checkIfJobStillRunning
#'
#' @param targetJob (character) The job name
#' @param targetValue (character) The return variable name
#' @param verbose (logical) Whether a message will be generated when complete
#' @param targetEnvironment (environment) Target environment
#'
#' @return logical; TRUE if still running; FALSE if not running
checkIfJobStillRunning <- function(targetJob, targetValue, verbose, targetEnvironment) {
  # Job is only still available for collection if it is in .mcparallelDoManager$runningJobs
  if (targetJob %in% names(.mcparallelDoManager$runningJobs)) {
    jobResult <- parallel::mccollect(get(targetJob, envir = targetEnvironment), wait=FALSE)
    if(is.null(jobResult)) {
      return(TRUE)
    } else {
      rm(list = targetJob, envir = targetEnvironment)
      assign(targetValue, jobResult[[1]], envir = targetEnvironment)
      if (verbose) message(targetValue, " has a new value")
      .mcparallelDoManager$removeJob(targetJob)
      return(FALSE)
    }
  } else {
    return(FALSE)
  }
}
NULL

#' jobCompleteDestructingHandler
#'
#' Creates a callback handler function that can be added via addTaskCallback().
#' These functions run at the end of each completed R statement.
#' This particular handler watches for the completion of the target job, which is created via mcparallel()
#' @param targetJob (character) Name of the mcparallel job variable that is waiting for a result
#' @param targetValue A character element indicating the variable that the result of that job should be assigned to targetEnvironment
#' @param verbose A boolean element; if TRUE the completion of the fork expr will be accompanied by a message
#' @param targetEnvironment The environment in which you want targetValue to be created
#'
#' @return callback handler function
jobCompleteSelfDestructingHandler <- function(targetJob, targetValue, verbose, targetEnvironment) {
  function(expr, value, ok, visible) {
    return(
      checkIfJobStillRunning(targetJob, targetValue, verbose, targetEnvironment)
    )
  }
}
NULL

#' mcparallelDo
#'
#' This function creates a fork, 
#' sets the variable named \code{targetValue} in the \code{targetEnvironment} to NULL,
#' evaluates a segment of code evaluated in the fork, 
#' and the result of the fork returned in a variable named \code{targetValue} in the \code{targetEnvironment} after the next top-level command completes.
#' If there is an error in the code, the returned variable will be a \code{try-error}.
#' These effects are accomplished via the automatic creation and destruction of a taskCallback and other functions inside the mcparallelDoManager.
#' If job results have to be collected before you return to the top level, use \link{mcparallelDoCheck}.
#' 
#' @param code The code to evaluate within a fork wrapped in {}
#' @param targetValue A character element indicating the variable that the result of that job should be assigned to targetEnvironment
#' @param verbose A boolean element; if TRUE the completion of the fork expr will be accompanied by a message
#' @param targetEnvironment The environment in which you want targetValue to be created
#'
#' @return The variable name of the job, this can be manually collected via mccollect or, if on Windows, an empty string
#'
#' @examples
#' ## Create data
#' data(ToothGrowth)
#' ## Trigger mcparallelDo to perform analysis on a fork
#' mcparallelDo({glm(len ~ supp * dose, data=ToothGrowth)},"interactionPredictorModel")
#' ## Do other things
#' binaryPredictorModel <- glm(len ~ supp, data=ToothGrowth)
#' gaussianPredictorModel <- glm(len ~ dose, data=ToothGrowth)
#' ## The result from mcparallelDo returns in your targetEnvironment, 
#' ## e.g. .GlobalEnv, when it is complete with a message (by default)
#' summary(interactionPredictorModel)
#' 
#' # Example of not returning a value until we return to the top level
#' for (i in 1:10) {
#'   if (i == 1) {
#'     mcparallelDo({2+2}, targetValue = "output")
#'   }
#'   if (exists("output")) print(i)
#' }
#' 
#' # Example of getting a value without returning to the top level
#' for (i in 1:10) {
#'   if (i == 1) {
#'     mcparallelDo({2+2}, targetValue = "output")
#'   }
#'   mcparallelDoCheck()
#'   if (exists("output")) print(i)
#' }
#' @importFrom ArgumentCheck addError finishArgCheck
#' @importFrom R.utils tempvar
#' @export
mcparallelDo <- function(code, targetValue, verbose = TRUE, targetEnvironment = .GlobalEnv) {
  Check <- ArgumentCheck::newArgCheck()
  if (!is.character(targetValue)) {
    ArgumentCheck::addError(
      msg = "targetValue must be a character",
      argcheck = Check
    )
  }
  if (length(targetValue) != 1) {
    ArgumentCheck::addError(
      msg = "targetValue must be a single element",
      argcheck = Check
    )
  }
  if (!is.environment(targetEnvironment)) {
    ArgumentCheck::addError(
      msg = "targetEnvironment must be an environment",
      argcheck = Check
    )
  }
  ArgumentCheck::finishArgCheck(Check)
  
  # Special handling for Windows
  if (!.Platform$OS.type=="unix") {
    assign(targetValue, value = {code}, envir = targetEnvironment)
    return("")
  }
  
  jobName <- R.utils::tempvar(".mcparallelDoJob", 
                              value = parallel::mcparallel({try(code)}),
                              envir = targetEnvironment)
  .mcparallelDoManager$addJob(jobName, targetValue, verbose, targetEnvironment)
  return(jobName)
}
NULL