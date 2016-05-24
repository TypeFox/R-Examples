
# calibrarDemo ------------------------------------------------------------

#' @title Demos for the calibrar package
#' @description Creates demo files able to be processed for a full calibration using
#' the calibrar package
#' 
#' @param path Path to create the demo files 
#' @param model Model to be used in the demo files, see details. 
#' @param \dots Additional parameters to be used in the construction of
#' the demo files.
#' @return A list with the following elements:
#' \item{path}{Path were the files were saved}
#' \item{par}{Real value of the parameters used in the demo} 
#' \item{constants}{Constants used in the demo} 
#' @author Ricardo Oliveros--Ramos
#' @references Oliveros-Ramos and Shin (2014)
#' @keywords demo calibration 
#' @examples
#' \dontrun{
#' require(calibrar)
#' set.seed(880820)
#' path = NULL # NULL to use the current directory
#' # create the demonstration files
#' demo = calibrarDemo(model="PoissonMixedModel", L=5, T=100) 
#' # get calibration information
#' calibrationInfo = getCalibrationInfo(path=demo$path)
#' # get observed data
#' observed = getObservedData(info=calibrationInfo, path=demo$path)
#' # read forcings for the model
#' forcing = read.csv(file.path(demo$path, "master", "environment.csv"), row.names=1)
#' # Defining 'runModel' function
#' runModel = function(par, forcing) {
#' output = calibrar:::.PoissonMixedModel(par=par, forcing=forcing)
#' # adding gamma parameters for penalties
#' output = c(output, list(gammas=par$gamma)) 
#' return(output)
#' }
#' # real parameters
#' cat("Real parameters used to simulate data\n")
#' print(demo$par)
#' # objective functions
#' obj  = createObjectiveFunction(runModel=runModel, info=calibrationInfo, 
#'                                observed=observed, forcing=forcing)
#' cat("Starting calibration...\n")
#' control = list(weights=calibrationInfo$weights, maxit=3.6e5) # control parameters
#' cat("Running optimization algorithms\n", "\t", date(), "\n")
#' cat("Running optim AHR-ES\n")
#' ahr = calibrate(par=demo$guess, fn=obj, lower=demo$lower, upper=demo$upper, control=control)
#' summary(ahr)
#' } 
#' @export 
calibrarDemo = function(path=NULL, model=NULL,  ...) {
  
  if(is.null(path)) path = getwd()
  if(is.null(model)) {
    model = "default"
    warning("Using default demo 'PoissonMixedModel'")
  }
  
  output = switch(model, 
                  PoissonMixedModel = .generatePoissonMixedModel(path=path, ...),
                  PredatorPrey      = .generatePredatorPreyModel(path=path, ...),
                  .generatePoissonMixedModel(path=path, ...)  
  )
  output$value = NA
  output$time = NA
  output$counts = c('function'=NA, gradient=NA)
  class(output) = c("calibrar.demo", "calibrar.results", class(output))
  return(output)                
  
}

#' @export
#' @method print calibrar.demo
print.calibrar.demo = function(x, ...) {
  print.default(x, ...)
}

