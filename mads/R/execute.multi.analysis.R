#' Performs Multiple Analyses on Distance Data
#'
#' Analyses are performed for multiple species contained within the same
#' dataset. Individual detection function analyses of each species must have
#' already been completed using the \code{ddf} function in the \code{mrds}
#' library. This function may then perform additional tasks such as assessing
#' variance via a non-parametric bootstrap, including covariate variability via
#' a parametric bootstrap, including model uncertainty and dealing with species
#' codes which relate to unidentified sightings.
#'
#' This is a new package with limited testing on real data, please drop
#' me a line if you plan on using it (lhm[at]st-and.ac.uk).
#'
#' The model fitting code in this function obtains its data and the model
#' descriptions from the ddf objects passed in via the \code{ddf.models} argument.
#'
#' If you wish to include model uncertainty then each model which you wish to
#' be included in the analyses must have already been run and should be
#' provided in the \code{ddf.models} argument. The \code{model.names} argument
#' tells this function which \code{"ddf"} objects are
#' associated with which species code in the dataset. This object must be
#' constructed as a list of vectors. Each element in the list must be named
#' corresponding to one of the species codes in the dataset and contain a
#' character vector of object names.
#'
#' For the majority of analyses the variance will be estimated using a
#' non-parametric bootstrap, indicated by the \code{bootstrap} argument. You
#' may select options for the bootstrap using the \code{bootstrap.options}
#' argument. This is a list with elements specifying the number of repetitions
#' and whether to resample samples within strata (\code{$resample = "samples"})
#' or observations withing strata (\code{$resample = "observations"}). In
#' addition, the \code{bootstrap.covariates} is a boolean argument specifying
#' whether or not a parametric bootstrap should be performed on any of the
#' covariates. The details of which variables should be resampled and from
#' which distributions should be entered in the \code{covariate.uncertainty}
#' dataframe. This dataframe should contain 7 columns with the following names:
#' \code{variable.layer}, \code{variable.name},
#' \code{cor.factor.layer}, \code{cor.factor.name}, \code{uncertainty.layer},
#' \code{uncertainty.name}, \code{uncertainty.measure} and
#' \code{sampling.distribution}. [Currently this is only implemented for the
#' observation layer]. The \code{variable.name} and
#' \code{uncertainty.name} should be the names of the variable in the dataset
#' giving the covariate to be resampled and the variable containing the
#' uncertainty respectively. The \code{cor.factor.layer} specifies the data
#' layer which contains the correction factor variable, although alternatively
#' "numeric" can be entered. The \code{cor.factor.name} specifies the name of
#' the correction factor variable or the correction factor value if "numeric"
#' was specified for the correction factor layer.
#' The \code{uncertainty.name} should specify what
#' values the uncertainty variable contains and should be one of \code{"sd"},
#' \code{"var"} or \code{"CV"}. The \code{sampling.distribution} should specify
#' one of the following distributions to parametrically resample from
#' \code{"Normal"}, \code{"Normal.Absolute"}, \code{"Lognormal.BC"},
#' \code{"Poisson"} or \code{"TruncPoissonBC"}. The remaning column in this
#' dataset, \code{variable.correction.factos}, allows the user to specify a
#' value by which the variable should be scaled. If this is not reqied this
#' should be set to 1.
#'
#' If there are unidentified sightings in the dataset then the
#' \code{unidentified.sightings} argument should be \code{true} and a
#' \code{species.code.definitions} list should be provided. This list must
#' contain one element for every unidentified species code which should be
#' named according to this code. Each element will contain a vector of
#' identified species codes corresponding to those species which the
#' unidentified code could have potentially been. This function uses this
#' information to prorate the abundance estimated from the unidentified species
#' codes to the relevant abundances from the identified codes. The prorating is
#' done individually for each strata. The function can be forced not to prorate
#' to any given species in any selected strata using the \code{species.presence}
#' argument. This is a list containing one element for each strata, each must be
#' named using the appropriate strata name. Each element should contain a vector
#' of identified species codes corresponding to which species are present in
#' each strata.
#'
#' @param species.code vector of all the species codes to be included in
#'   the analysis
#' @param unidentified.sightings  a list with an element for each
#'   unidentified code which contains a vector of corresponding identified
#'   species codes or NULL if not required
#' @param species.presence must be specified if species.code.definitions is
#'  specified. A list with an element for each strata which contains the vector
#'  of species codes present in that strata
#' @param covariate.uncertainty a dataframe detailing the variables to be
#'   resampled - variable.layer, variable.name, cor.factor.layer,
#'   cor.factor.name , uncertainty.layer, uncertainty.name,
#'   uncertainty.measure, sampling.distribution. or NULL if not required
#' @param models.by.species.code a list of character vectors of model names
#'   with the elements named by species code
#' @param ddf.model.objects a list of all the ddf models named in models.by.species.code
#' @param ddf.model.options a list of options 1) selection.criterion either "AIC",
#'   "AICc" or "BIC" 2) species.field.name describing the field name in the ddf
#'   dataset containing species codes.
#' @param region.table dataframe of region records - Region.Label and Area
#' @param sample.table dataframe of sample records - Region.Label,
#'   Sample.Label, Effort
#' @param obs.table dataframe of observation records with fields object,
#'   Region.Label, and Sample.Label which give links to sample.table,
#'   region.table and the data records used in \code{model}
#' @param dht.options list containing option for dht: convert.units indicated
#'   if the distance measurement units are different from shapefile and transect
#'   coordinate units.
#' @param bootstrap if TRUE resamples data to obtain variance estimate
#' @param bootstrap.options a list of options that can be set 1) n: number of
#'   repetitions 2) resample: how to resample data ("samples", "observations")
#' @param silent boolean used to surpress some warning messages
#' @return object of class "ma" which consists of a list of objects of class
#'   "ma.element". Each "ma.element" consists of the following elements:
#'   \item{individuals}{Summary, N (abundance) and D (density) tables}
#'   \item{clusters}{Summary, N (abundance) and D (density) tables}
#'   \item{Expected.S}{Expected cluster size table}
#'   \item{ddf}{Model details including a summary of convergence and selection
#'     as well as parameter estimates for selected models.}
#' @export
#' @author Laura Marshall
#' @references
#'   Marques, F.F.C. and S.T. Buckland. 2004. Covariate models for the detection
#'     function. In: Advanced Distance Sampling, eds. S.T. Buckland,
#'     D.R.Anderson, K.P. Burnham, J.L. Laake, D.L. Borchers, and L. Thomas.
#'     Oxford University Press.
#'   Gerrodette, T. and Forcada, J. 2005 Non-recovery of two spotted and spinner
#'     dolphin populations in the eastern tropical Pacific Ocean. Marine Ecology
#'     Progress Series, 291:1-21.
#' @keywords ~distance sampling, unidentified sightings, covariate uncertainty, model uncertainty
#' @examples
#'
#' \dontrun{
#'
#' ex.filename<-system.file("testData/input_checks/ddf_dat.robj", package="mads")
#' load(ex.filename)
#' ex.filename<-system.file("testData/input_checks/obs_table.robj", package="mads")
#' load(ex.filename)
#' ex.filename<-system.file("testData/input_checks/region_table.robj", package="mads")
#' load(ex.filename)
#' ex.filename<-system.file("testData/input_checks/sample_table.robj", package="mads")
#' load(ex.filename)
#'
#' #run ddf analyses
#' ddf.1 <- ddf(dsmodel = ~mcds(key = "hn", formula = ~ size),
#'              method='ds', data=ddf.dat,meta.data=list(width=4))
#' ddf.2 <- ddf(dsmodel = ~mcds(key = "hr", formula = ~ size),
#'              method='ds', data=ddf.dat,meta.data=list(width=4))
#'
#' model.names <- list("CD"=c("ddf.1","ddf.2"), "WD"=c("ddf.1","ddf.2"),
#'                    "UnidDol"=c("ddf.1","ddf.2"))
#' ddf.models  <- list("ddf.1" = ddf.1, "ddf.2" = ddf.2)
#'
#' unidentified.code.definitions <- list("UnidDol" = c("CD","WD"))
#' bootstrap.options             <- list(resample="samples", n=10, quantile.type = 7)
#'
#' results<- execute.multi.analysis(
#'              species.code = names(model.names),
#'              unidentified.sightings = unidentified.code.definitions,
#'              models.by.species.code = model.names,
#'              ddf.model.objects = ddf.models,
#'              ddf.model.options = list(criterion="AIC"),
#'              region.table = region.table,
#'              sample.table = sample.table,
#'              obs.table = obs.table,
#'              bootstrap = TRUE,
#'              bootstrap.option = bootstrap.options)
#'
#' }
#'
execute.multi.analysis <- function(species.code, unidentified.sightings = NULL, species.presence = NULL, covariate.uncertainty = NULL, models.by.species.code, ddf.model.objects, ddf.model.options = list(criterion="AIC", species.field.name = "species"), region.table, sample.table, obs.table, dht.options = list(convert.units = 1), bootstrap, bootstrap.options = list(resample="samples", n = 1, quantile.type = 7), silent = FALSE){

  #create global variable to store error messages
  MAE.warnings <- NULL

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #INPUT CHECKS AND INPUT PROCESSING
  #ddf models and model options
  ddf.model.info           <- check.ddf.models(models.by.species.code, ddf.model.objects)
  clusters                 <- ddf.model.info$clusters
  double.observer          <- ddf.model.info$double.observer
  # If the user has not specified the criteria set it
  if(is.null(ddf.model.options$criterion)){
    ddf.model.options$criterion <- "AIC"
  }
  # If the user has not specified the species field name set it
  if(is.null(ddf.model.options$species.field.name)){
    ddf.model.options$species.field.name <- "species"
  }

  #Species codes and unidentified sightings
  if(!is.null(unidentified.sightings)){
    species.code.definitions <- check.species.code.definitions(unidentified.sightings, species.code)
  }else{
    temp <- list()
    temp[[species.code]] <- species.code
    species.code.definitions <- list(unidentified = FALSE, species.code.definitions = temp)
  }
  unidentified.species     <- species.code.definitions$unidentified
  species.code.definitions <- species.code.definitions$species.code.definitions

  #Species presence - if null it populates it
  species.presence         <- check.species.presence(species.presence, species.code, strata.name = as.character(region.table$Region.Label))

  #Covariate uncertainty
  if(!is.null(covariate.uncertainty)){
    covariate.uncertainty    <- check.covar.uncertainty(covariate.uncertainty)
  }

  #Bootstrap Options
  if(is.null(bootstrap.options$resample)){
    bootstrap.options$resample <- "samples"
  }
  if(is.null(bootstrap.options$n)){
    bootstrap.options$n <- 1
  }
  if(is.null(bootstrap.options$quantile.type)){
    bootstrap.options$quantile.type <- 7
  }
  check.bootstrap.options(bootstrap, bootstrap.options$resample, bootstrap.options$n, sample.table)
  bootstrap.options$n <- ifelse(bootstrap, bootstrap.options$n, 1)
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

  #Make master copies of all the datasets
  ddf.dat.master      <- get.datasets(models.by.species.code, ddf.model.objects)
  unique.model.names  <- ddf.dat.master$unique.model.names
  model.index         <- ddf.dat.master$model.index
  ddf.dat.master      <- ddf.dat.master$ddf.dat.master
  obs.table.master    <- obs.table
  sample.table.master <- sample.table

  #Create storage for results (only for the species codes not the unidentified codes)
  bootstrap.results <- create.result.arrays(species.code, species.code.definitions, region.table, clusters, bootstrap.options$n)
  bootstrap.ddf.statistics <- create.param.arrays(unique.model.names, ddf.model.objects, bootstrap.options$n, ddf.model.options$criterion)

  #Set up a loop
  for(n in 1:bootstrap.options$n){

      #if(!is.null(seed.array)){
      #  set.seet(seed.array[n])
      #}

      #Resample Data
      if(bootstrap){
        ddf.dat.working <- resample.data(resample=bootstrap.options$resample, obs.table.master, sample.table.master, ddf.dat.master, double.observer)
        obs.table       <- ddf.dat.working$obs.table
        sample.table    <- ddf.dat.working$sample.table
        ddf.dat.working <- ddf.dat.working$ddf.dat.working
      }else{
        ddf.dat.working <- ddf.dat.master
      }

      #Add uncertainty to covariates
      if(!is.null(covariate.uncertainty)){
        ddf.dat.working <- resample.covariates(ddf.dat.working, covariate.uncertainty, MAE.warnings)
        MAE.warnings <- ddf.dat.working$MAE.warnings
        ddf.dat.working <- ddf.dat.working$ddf.dat.working
      }

      #Fit ddf models to all species codes
      ddf.results <- fit.ddf.models(ddf.dat.working, unique.model.names, ddf.model.objects, ddf.model.options$criterion, bootstrap.ddf.statistics, n, MAE.warnings)
      if(length(ddf.results) > 1){
        bootstrap.ddf.statistics <- ddf.results$bootstrap.ddf.statistics
        ddf.results <- ddf.results$ddf.results
        MAE.warnings <- ddf.results$mae.warnings
      }else{
        #If the ddf results are not valid for all species move to next bootstrap iteration
        MAE.warnings <- ddf.results$mae.warnings
        next
      }

      #Calculate densities and abundance for all species codes
      dht.results <- calculate.dht(species.code, ddf.model.options$species.field.name, model.index, ddf.results, region.table, sample.table, obs.table, dht.options)

      #Deal with unidentified sightings if present or format dht results if not
      if(unidentified.species){
        formatted.dht.results <- prorate.unidentified(dht.results, species.code.definitions, species.presence, clusters)
      }else{
        formatted.dht.results <- format.dht.results(dht.results, species.code, clusters)
      }

      #Format / Record results
      bootstrap.results <- accumulate.results(n, bootstrap.results, formatted.dht.results, clusters)

  }#next iteration

  #Debugging
  #save(bootstrap.results, file = "/Users/laura/mads/bootstrap.results.robj")
  #cat("saving bootstrap results")

  #process results
  results <- process.bootstrap.results(bootstrap.results, model.index, clusters, bootstrap.ddf.statistics, bootstrap.options$quantile.type, analysis.options = list(bootstrap = bootstrap, n = bootstrap.options$n, covariate.uncertainty = covariate.uncertainty, clusters = clusters, double.observer = double.observer, unidentified.species = unidentified.species, species.code.definitions = species.code.definitions, model.names = models.by.species.code))
  #process warning messages
  process.warnings(MAE.warnings)

  #return results
  class(results) <- "ma"
  class(results$analysis.options) <- "ma.analysis"
  class(results$species) <- "ma.allspecies"
  for(sp in seq(along = results$species)){
    class(results$species[[sp]]) <- "ma.species"
  }
  if(!is.null(results$unidentified)){
    class(results$unidentified) <- "ma.allunid"
    for(sp in seq(along = results$unidentified)){
      class(results$unidentified[[sp]]) <- "ma.unid"
    }
  }

  return(results)
}












