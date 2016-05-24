#' Refits the detection functions to the resampled data
#'
#' Fits all the models named in model.names to the associated data supplied in 
#' ddf.dat.working. If more than one model is supplied for any species the 
#' model with the minimum selection crieteria will be selected.
#'  
#' @param ddf.dat.working list of dataframes containing the data to which the 
#'   models will be fitted
#' @param model.names list of unique character vectors giving the names of the 
#'   ddf objects for each species.
#' @param ddf.models a list of ddf objects
#' @param criterion character option specifying the model selection criteria - 
#'   "AIC", "AICc" or "BIC".
#' @param bootstrap.ddf.statistics array storing parameter estimates
#' @param rep.no numeric value indicating iteration number
#' @param MAE.warnings character vector of warning messages
#' @return list of ddf objects 
#' @note Internal function not intended to be called by user.
#' @author Laura Marshall
#' @importFrom stats na.omit
#'
fit.ddf.models <- function(ddf.dat.working, model.names, ddf.models, criterion, bootstrap.ddf.statistics, rep.no, MAE.warnings){
# fit.ddf.models function to refits the detection functions to the data provided
#
# Arguments:
#   ddf.dat.working list of dataframes containing the datasets
#   model.names list of unique character vectors specifying the model names
#   ddf.models a list of ddf objects
#   criterion character model selection option ("AIC", "AICc" or "BIC")
#
# Value: list of ddf objects 
#
# Function Calls: mrds::ddf
#
  species.name <- names(model.names)
  #create storage list
  ddf.results <- list()
  #for every species code
  for(sp in seq(along = species.name)){
    temp.results <- list()
    selection.criterion.values <- NULL
    #get dataset
    usedata <- ddf.dat.working[[species.name[sp]]]
    #for every model
    for(m in seq(along = model.names[[species.name[sp]]])){
      #get model name
      current.model.name <- model.names[[species.name[sp]]][m]  
      current.model <- ddf.models[[current.model.name]] 
      #get model call                        
      model.call <- current.model$call
      #point the call at the new data
      model.call$data <- as.name("usedata")
      #get parameter estimates from original model to act as start values to aid convergence
      model.type <- ddf.models[[current.model.name]]$method
      if(model.type%in%c("ds")){
        start.values <- list(scale=current.model$ds$aux$ddfob$scale$parameters, shape=current.model$ds$aux$ddfob$shape$parameters, adjustment=current.model$ds$aux$ddfob$adjustment$parameters)                           
        model.call$control <- call("list",initial=start.values)
      }
      if(model.type%in%c("trial", "io")){
        start.values <- list(scale=current.model$ds$ds$aux$ddfob$scale$parameters, shape=current.model$ds$ds$aux$ddfob$shape$parameters, adjustment=current.model$ds$ds$aux$ddfob$adjustment$parameters)                           
        model.call$control <- call("list",initial=start.values)
      }
      #refit ddf model
      options(show.error.messages = FALSE)
      temp.results[[m]] <- try(eval(model.call), silent = TRUE)
      options(show.error.messages = TRUE)
      if(any(class(temp.results[[m]]) == "try-error")){
        #Model fitting threw an error
        bootstrap.ddf.statistics[[species.name[sp]]]$convergence[2,current.model.name] <- bootstrap.ddf.statistics[[species.name[sp]]]$convergence[2,current.model.name] + 1
        MAE.warnings <- mae.warning(paste("Model was not successfull for species ",species.name[sp]," model ",current.model.name,".", sep = ""), warning.mode="store", MAE.warnings)
      }else if(check.convergence(temp.results[[m]])){
        #Model converged save info
        bootstrap.ddf.statistics[[species.name[sp]]]$convergence[1,current.model.name] <- bootstrap.ddf.statistics[[species.name[sp]]]$convergence[1,current.model.name] + 1
        bootstrap.ddf.statistics <- store.param.ests(bootstrap.ddf.statistics, species.name[sp], current.model.name, temp.results[[m]], rep.no)
        lnl <- temp.results[[m]]$lnl 
        k <- length(temp.results[[m]]$par)
        n <- nrow(temp.results[[m]]$data)
        selection.criterion.values[m] <- switch(criterion,
          AIC  = 2*k-2*lnl,
          AICc = 2*k-2*lnl+(2*k*(k+1))/(n-k-1),
          BIC  = k*log(n)-2*lnl)
        bootstrap.ddf.statistics[[species.name[sp]]][[current.model.name]][[criterion]][rep.no] <- selection.criterion.values[m]             
      }else{
        #Model failed to converge - this is an error so never gets here
        MAE.warnings <- mae.warning(paste("Model did not converge for species ",species.name[sp]," model ",current.model.name,". Convergence code was not zero.", sep = ""), warning.mode="store", MAE.warnings)
        bootstrap.ddf.statistics[[species.name[sp]]]$convergence[2,current.model.name] <- bootstrap.ddf.statistics[[species.name[sp]]]$convergence[2,current.model.name] + 1        
      }
    }#next model 
    #Check at least one converged                  
    if(!is.null(selection.criterion.values)){  
      selected.model <- which(selection.criterion.values == min(na.omit(selection.criterion.values)))  
      ddf.results[[species.name[sp]]] <- temp.results[[selected.model[1]]]
      selected.model.name <- model.names[[species.name[sp]]][selected.model[1]]
      bootstrap.ddf.statistics[[species.name[sp]]]$convergence[3,selected.model.name] <- bootstrap.ddf.statistics[[species.name[sp]]]$convergence[3,selected.model.name] + 1
      bootstrap.ddf.statistics[[species.name[sp]]][[selected.model.name]]$selected[rep.no] <- 1
    }else{
      #If none converged return NULL and exit function
      MAE.warnings <- mae.warning(paste("No models converged for species ",species.name[sp],", this bootstrap iteration is being skipped", sep=""), warning.mode="store", MAE.warnings)
      return(list(mae.warnings = MAE.warnings))
    }
  }#next species
  return(list(ddf.results = ddf.results, bootstrap.ddf.statistics = bootstrap.ddf.statistics, mae.warnings = MAE.warnings))
}



