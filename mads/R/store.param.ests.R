#' Updates bootstrap.ddf.statistics 
#'
#' #' Updates bootstrap.ddf.statistics with the latest parameter estimates from
#'   the ddf object supplied.
#'
#' @param bootstrap.ddf.statistics a list containing various arrays and vectors
#'   storing model details.
#' @param species.name name of species the model relates to
#' @param model.name name of ddf object
#' @param ddf.model ddf object
#' @param rep.no iteration number of the bootstrap
#' @return bootstrap.ddf.statistics the updated list
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @seealso \code{create.param.arrays}
#' @keywords data manipulation
#'        
store.param.ests <- function(bootstrap.ddf.statistics, species.name, model.name, ddf.model, rep.no){
# store.param.ests function to non-parametrically resample the observations
#
# Arguments:
#   bootstrap.ddf.statistics   - list of arrays and vectors
#   species.name               - name of species model relates to
#   model.name                 - name of ddf object
#   rep.no                     - iteration number of the bootstrap
#
# Value: bootstrap.ddf.statistics  - the updated list
#
# Functions Used: none
#
  model.type <- ddf.model$method
  if(model.type%in%c("ds")){
    param.ests <- c(ddf.model$ds$aux$ddfob$scale$parameters, ddf.model$ds$aux$ddfob$shape$parameters, ddf.model$ds$aux$ddfob$adjustment$parameters)
    bootstrap.ddf.statistics[[species.name]][[model.name]]$ds.param[rep.no, 1:length(param.ests)] <- param.ests
  }
  if(model.type%in%c("trial", "io")){
    param.ests <- c(ddf.model$ds$ds$aux$ddfob$scale$parameters, ddf.model$ds$ds$aux$ddfob$shape$parameters, ddf.model$ds$ds$aux$ddfob$adjustment$parameters)
    bootstrap.ddf.statistics[[species.name]][[model.name]]$ds.param[rep.no, 1:length(param.ests)] <- param.ests
  }
  if(model.type%in%c("trial", "io")){
    param.ests <- ddf.model$mr$mr$coefficients
    bootstrap.ddf.statistics[[species.name]][[model.name]]$mr.param[rep.no, 1:length(param.ests)] <- param.ests
  }
  if(model.type%in%c("trial.fi", "io.fi")){
    param.ests <- ddf.model$mr$coefficients
    bootstrap.ddf.statistics[[species.name]][[model.name]]$mr.param[rep.no, 1:length(param.ests)] <- param.ests
  }
  return(bootstrap.ddf.statistics)
}


