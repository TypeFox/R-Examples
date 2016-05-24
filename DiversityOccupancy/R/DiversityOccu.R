globalVariables(c("data2"))
#' Fits occupancy models for multiple species detection history
#'
#' This function takes a data.frame with multiple detection history from
#' various species in different sites, covariates of each site to calculate
#' occupancy, variables specific to sampling days to calculate probability of
#' detection. It features an automatic model selection when dredge = TRUE.
#'
#' @param pres a data.frame where rows are the sites and columns are a series of
#' presence-absence observation from multiple species, every species needs to
#' have the same number of observations.
#' @param sitecov a data.frame where every row is a site, and every column is a
#' measurement of that site, such as elevation or slope, this covariates are
#' usually more constant.
#' @param obscov a list where every element is a data frame with the daily
#' covariates for each site, that is a measurement for each day, such as average
#' temperature of a day, this covariates are usually very .
#' @param spp the number of species in the pres data.frame
#' @param form a formula in the format ~ obscov ~ sitcov, the first arguments
#' will be used to calculate probability of detection and the second part the
#' occupancy.
#' @param dredge default = FALSE, if TRUE, for each species, the best occupancy
#' model will be determined by fitting all possible models and ranking by AICc.
#' @param pos where to do the removal. By default, uses the current environment.
#' @param envir the environment to use.
#' @return A list with the fitted models for each species and the calculated
#' Alpha diversity for each site.
#' @details
#' This function fits the single season occupancy model of MacKenzie et al (2002),
#' for multiple species and it can automatically select the best model for each
#' specie based on AICc.
#' @examples
#' \dontrun{
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' BatOccupancy <-batchoccu(pres = BatOccu, sitecov = sampling.cov, obscov =
#' Dailycov,spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum + sdtemp ~
#' Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy +
#' I(Burn.intensity.Canopy^2) + Burn.intensity.basal + I(Burn.intensity.basal^2))
#' #plot the response of occupancy to individual variables for species 4, 11
#' #and 15
#'
#' responseplot.occu(batch = BatOccupancy, spp = 4, variable = Burn.intensity.soil)
#'
#' responseplot.occu(batch = BatOccupancy, spp = 11, variable = Burn.intensity.soil)
#'
#' responseplot.occu(batch = BatOccupancy, spp = 15, variable = Burn.intensity.soil)
#' }
#' #Dredge for 2 species
#' A <- batchoccu(pres = BatOccu[,1:6], sitecov = sampling.cov, obscov = Dailycov,
#' spp = 2, form = ~ Meanhum + Meantemp ~  Burn.intensity.basal +
#' I(Burn.intensity.basal^2), dredge = TRUE)
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @export
#' @importFrom unmarked occu
#' @importFrom unmarked unmarkedFrameOccu
#' @importMethodsFrom unmarked predict
#' @importFrom MuMIn dredge
#' @importFrom MuMIn get.models
#' @importFrom MuMIn AICc
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

batchoccu<- function(pres, sitecov, obscov, spp, form, dredge = FALSE,  pos = 1, envir = as.environment(pos)) {
  secuencia <- c(1:spp)*(ncol(pres)/spp)
  secuencia2<-secuencia-(secuencia[1]-1)
  models <- list()
  data<- list()
  fit <- list()
  dredged <- list ()
  if (dredge == FALSE){

    for(i in 1:length(secuencia)) {
      data[[i]] <-c(secuencia2[i]:secuencia[i])
      data[[i]] <- pres[, data[[i]]]
      models[[i]] <- unmarkedFrameOccu(y = data[[i]], siteCovs = sitecov, obsCovs = obscov)
      models[[i]] <- occu(form, models[[i]])
      fit[[i]] <- predict(models[[i]], type = "state", newdata = sitecov)$Predicted
      fit <- as.data.frame(fit)
      colnames(fit) = paste("species",c(1:ncol(fit)), sep =".")
    }
  }

  else if (dredge==TRUE) {
    for(i in 1:length(secuencia)) {
      data[[i]] <-c(secuencia2[i]:secuencia[i])
      data[[i]] <- pres[, data[[i]]]
      #data is a list of class unmarkedFrames from package unmarked.
      # NM: write to the global environment so he data won't be "lost"
      assign("data2",  unmarkedFrameOccu(y = data[[i]], siteCovs = sitecov, obsCovs = obscov), envir = envir)
      models[[i]] <- occu(form, data2)
      #selects models
      # NM: saved this to dredged object rather than overwriting models object
      dredged[[i]] <- dredge(models[[i]], data2)
      #select the first model
      models[[i]] <- get.models(dredged[[i]], 1)[[1]]
      #predictions for the best model
      fit[[i]] <- predict(models[[i]], type = "state", newdata = sitecov)$Predicted
      fit<- as.data.frame(fit)
      colnames(fit) = paste("species",c(1:ncol(fit)), sep =".")
      data[[i]] <- data2
    }
    # remove temporary data file from the global environment
    rm(data2, pos= envir)
  }

  result <- list(Covs = sitecov, models = models, fit = fit)
  class(result)<- "batchoccupancy"
  return(result)
}


#' Calculates alpha diversity from multiple species occupancy data
#'
#' This function takes a data.frame with multiple presence absence-data from
#' various species in different sites, covariates of each site to calculate
#' occupancy, variables specific to sampling days to calculate probability of
#' detection, and it calculates the alpha diversity for each site.
#'
#' @param pres a data.frame where rows are the sites and columns are a series of
#' presence-absence observation from multiple species, every species needs to
#' have the same number of observations.
#' @param sitecov a data.frame where every row is a site, and every column is a
#' measurement of that site, such as elevation or slope, this covariates are
#' usually more constant.
#' @param obscov a list where every element is a data frame with the daily
#' covariates for each site, that is a measurement for each day, such as average
#' temperature of a day, this covariates are usually very .
#' @param spp the number of species in the pres data.frame
#' @param form a formula in the format ~ obscov ~ sitcov, the first arguments
#' will be used to calculate probability of detection and the second part the
#' occupancy.
#' @param index Diversity index, one of "shannon", "simpson" or "invsimpson".
#' @param dredge default = FALSE, if TRUE, for each species, the best occupancy
#' model will be determined by fitting all possible models and ranking by AICc.
#' @param pos where to do the removal. By default, uses the current environment.
#' @param envir the environment to use.
#' @return A list with the fitted models for each species, the calculated
#' Alpha diversity for each site, and a dataframe with the abundance of each
#' species and diversity.
#' @details
#' This function fits the latent abundance mixture model described in Royle and
#' Nichols (2003), to calculate the abundance of every species in each site, the
#' using that abundance it calculates the alpha diversity index for each site
#' based on that abundance.
#' @examples
#' \dontrun{
#' #Load the data
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#'
#' #Model the abundance for 17 bat species and calculate alpha diversity from that
#'
#' BatDiversity <-diversityoccu(pres = BatOccu, sitecov = sampling.cov,
#' obscov = Dailycov,spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum +
#' sdtemp ~ Burn.intensity.soil + I(Burn.intensity.soil^2) +
#' Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#'
#' #To see the estimates and p values for each model:
#'
#' BatDiversity$models
#' }
#' @seealso \code{\link[vegan]{diversity}}
#' @seealso \code{\link[DiversityOccupancy]{model.diversity}}
#' @export
#' @importFrom vegan diversity
#' @importFrom unmarked occuRN
#' @importFrom unmarked unmarkedFrameOccu
#' @importMethodsFrom unmarked predict
#' @importFrom MuMIn dredge
#' @importFrom MuMIn get.models
#' @importFrom MuMIn AICc

#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>
#' @author Nicole L. Michel

diversityoccu<- function(pres, sitecov, obscov, spp, form, index = "shannon", dredge = FALSE,  pos = 1, envir = as.environment(pos)) {
  secuencia <- c(1:spp)*(ncol(pres)/spp)
  secuencia2<-secuencia-(secuencia[1]-1)

  models <- list()
  data<- list()
  div <- list()
  dredged <- list ()
  if (dredge == FALSE){

    for(i in 1:length(secuencia)) {
      data[[i]] <-c(secuencia2[i]:secuencia[i])
      data[[i]] <- pres[, data[[i]]]
      models[[i]] <- unmarkedFrameOccu(y = data[[i]], siteCovs = sitecov, obsCovs = obscov)
      models[[i]] <- occuRN(form, models[[i]])
      div[[i]] <- predict(models[[i]], type = "state", newdata = sitecov)$Predicted
      div<- as.data.frame(div)
      colnames(div) = paste("species",c(1:ncol(div)), sep =".")
      h<- diversity(div, index)
      DF <- cbind(h, div)
    }
  }

  else if (dredge==TRUE) {
    for(i in 1:length(secuencia)) {
      data[[i]] <-c(secuencia2[i]:secuencia[i])
      data[[i]] <- pres[, data[[i]]]
      #data is a list of class unmarkedFrames from package unmarked.
      # NM: write to the global environment so he data won't be "lost"
      assign("data2",  unmarkedFrameOccu(y = data[[i]], siteCovs = sitecov, obsCovs = obscov), envir = envir)
      models[[i]] <- occuRN(form, data2)
      #selects models
      # NM: saved this to dredged object rather than overwriting models object
      dredged[[i]] <- dredge(models[[i]], data2)
      #select the first model
      models[[i]] <- get.models(dredged[[i]], 1)[[1]]
      #predictions for the best model
      div[[i]] <- predict(models[[i]], type = "state", newdata = sitecov)$Predicted
      div<- as.data.frame(div)
      colnames(div) = paste("species",c(1:ncol(div)), sep =".")
      data[[i]] <- data2
      h <- diversity(div, index)
      DF <- cbind(h, div)
    }
    # remove temporary data file from the global environment
    rm(data2, pos=envir)
  }

  result <- list(Covs = sitecov, models = models, Diversity = h, species = DF)
  class(result)<- "diversityoccupancy"

  return(result)
}

#' Find the best GLM model explaining the alpha divesity of the species
#'
#' This function takes a diversityoccu object and heuristically searches for the
#' glm that best explains the alpha diversity of the modelled species.
#'
#' @param DivOcc is an object returned by the divesityoccu function of this
#' package
#' @param method The method to be used to explore the candidate set of models.
#' If "h" an exhaustive screening is undertaken. If "g" the genetic algorithm is
#' employed (recommended for large candidate sets). If "l", a very fast
#' exhaustive branch-and-bound algorithm is used. Package leaps must then be
#' loaded, and this can only be applied to linear models with covariates and no
#' interactions.
#' @param delta	The number of models that will be returned will be the ones that
#' have a maximum AICc difference with the top model equal to delta.
#' @param squared, if FALSE (Default), only GLMs with linear components will be
#' evaluated; If TRUE, GLMs with both linear and quadratic components will be evaluated.
#' WARNING if squared is TRUE, the number of parameters duplicates and the models
#' grow exponentially, this may result in to many variables for a CPU to compute.
#' @return An object with the best fitted model, the coefficients of that model,
#' a table with the top 5 fitted models ranked by AICc and the data used for the
#' model
#' @details
#' This function fits every first order glm possible and ranks them by AICc.
#' @examples
#' \dontrun{
#' #Load the data
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#'
#' #Model the abundance for 17 bat species and calculate alpha diversity from that
#'
#' BatDiversity <-diversityoccu(pres = BatOccu, sitecov = sampling.cov,
#' obscov = Dailycov,spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum +
#' sdtemp ~ Burn.intensity.soil + I(Burn.intensity.soil^2) +
#' Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#'
#' #Select the best model that explains diversity using genetic algorithms
#' set.seed(123)
#' glm.Batdiversity <- model.diversity(BatDiversity, method = "g")
#'
#' #see the best models
#'
#' glm.Batdiversity$Best.model
#'
#' #plot the response of diversity to individual variables
#'
#' plot(glm.Batdiversity, Burn.intensity.soil)
#'
#' #To add the quadratic components of models
#'
#' batdiversity <-diversityoccu(pres = BatOccu, sitecov = sampling.cov,
#' obscov = Dailycov, spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum +
#' sdtemp ~Burn.intensity.soil + I(Burn.intensity.soil^2) + Burn.intensity.Canopy +
#' I(Burn.intensity.Canopy^2) + Burn.intensity.basal +I(Burn.intensity.basal^2))
#' set.seed(123)
#' glm.batdiversity <- model.diversity(batdiversity , method = "g", squared = TRUE)
#'
#' responseplot.diver(glm.batdiversity, Burn.intensity.Canopy)
#' }
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @export
#' @importFrom stats glm
#' @importFrom stats as.formula
#' @importFrom glmulti glmulti
#' @importFrom glmulti weightable
#' @importFrom dplyr filter
#' @importFrom qpcR akaike.weights
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>


model.diversity <- function(DivOcc, method = "h", delta = 2, squared = FALSE){
  Delta.AICc <- NULL
  A <- cbind(DivOcc$Diversity, DivOcc$Covs)
  colnames(A)[1]<-"Diversity"
  B <- paste(names(DivOcc$Covs), "+")
  B <- toString(B)
  B <- gsub(",", " ", B)
  if (squared == TRUE) {
    C <- paste("I(", names(DivOcc$Covs), "^2) +")
    C <- toString(C)
    C <- gsub (",", " ", C)
    B <- paste("Diversity ~", B, C, collapse = " ")
  }
  else if (squared == FALSE) {
    B <- paste("Diversity ~", B, collapse = " ")
  }

  B <- as.formula(substr(B, 1, nchar(B)-1))
  B <- glm(B, data = A)
  D <- glmulti(B, level = 1, crit = "aicc", plotty = FALSE, method = method)
  Best.model <- D@formulas[[1]]
  Table <- weightable(D)
  Table$Delta.AICc <- Table[,2]-Table[1,2]
  Table <- filter(Table, Delta.AICc < delta)
  Table$weights <- akaike.weights(Table$aicc)$weights
  d<-summary(glm(Best.model, data = A))
  result <- list(Best_model = Best.model, Table = Table, coeff = d, dataset= A)
  class(result) <- "modeldiversity"
  return(result)
}

#' Makes a spacially explicit prediction of the occupancy of multiple species
#' and alpha diversity, and select the area where
#'
#' This function takes an deiversityoccu object and predicts occupancy for all species
#' in new data, either a data.frame or a rasterstack. It can also return a subset
#' of the total area of a rasterstack, where diversity and occupancy/abundance are
#' higher than the nth quantile.
#' @param model A result from diversityoccu
#' @param diverse A result from the model.diversity function.
#' @param new.data a rasterstack, or a dataframe containing the same variables as
#' the siteCovs variable in diversityoccu or batchoccu
#' @param quantile.nth the nth quantile, over which is a goal to keep both diversity
#' and selected species. default = NULL
#' @param species, a boolean vector of the species to take into acount
#' @return a data frame with predicted values, or a raster stack with predictions
#' for each species, a raster for diversity and a raster with the area meeting the
#' quantile criteria.
#' @examples
#' \dontrun{
#' #Load the data
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' data("plumas.stack")
#'
#' #Model the abundance for 17 bat species and calculate alpha diversity from that
#'
#' BatDiversity <-diversityoccu(pres = BatOccu, sitecov = sampling.cov,
#' obscov = Dailycov,spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum +
#' sdtemp ~ Burn.intensity.soil + I(Burn.intensity.soil^2) +
#' Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#'
#' #Select the best model that explains diversity using genetic algorithms
#' set.seed(123)
#' glm.Batdiversity <- model.diversity(BatDiversity, method = "g")
#'
#' # get the area where the first two bat species Myyu and Myca are most abundant
#' # and the diversity is most abundant
#'
#' Selected.area <- diversity.predict(model = BatDiversity, diverse = glm.Batdiversity,
#' new.data = plumas.stack, quantile.nth = 0.85, species =
#' c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,
#' FALSE, FALSE, FALSE, FALSE,FALSE,FALSE))
#'
#' Selected.area
#' }
#' @export
#' @seealso \code{\link[DiversityOccupancy]{diversityoccu}}
#' @seealso \code{\link[DiversityOccupancy]{batchoccu}}
#' @seealso \code{\link[DiversityOccupancy]{model.diversity}}
#' @importFrom graphics plot
#' @importFrom stats glm
#' @importFrom raster addLayer
#' @importFrom raster KML
#' @importFrom raster quantile
#' @importFrom raster reclassify
#' @importFrom raster stack
#' @importFrom raster subset
#' @importFrom raster unstack
#' @importFrom raster writeRaster
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>

diversity.predict<- function(model, diverse, new.data, quantile.nth = 0.8 , species) {
  models <- model$models[species]
  layers <- list()
  for (i in 1:length(models)){
    layers [[i]] <- predict(models[[i]], new.data, type = "state")$Predicted
  }
  glm.model <- glm(diverse$Best_model, data = diverse$dataset)
  diversity.raster<- predict(object = new.data, model = glm.model)
  layers <- stack(unlist(layers))
  desition <- addLayer(layers, diversity.raster)
  nths <- quantile(desition, quantile.nth)
  desition <- unstack(desition)
  rc<-list()
  for (i in 1:length(nths)){
    m <- c(-Inf, nths[i], NA,  nths[i], Inf, 1)
    rclmat <- matrix(m, ncol=3, byrow=TRUE)
    rc[[i]]<-reclassify(desition[[i]] , rcl = rclmat)
  }
  rc<- stack(unlist(rc))
  priority.area <- prod(rc)
  plot(priority.area, colNA="black", legend = FALSE)
  KML(priority.area, file='priority_area.kml', overwrite = TRUE, col = "red")
  result <- list(species = layers, diversity.raster = diversity.raster, priority.area = priority.area)
  return(result)
}

#' Predicts occupancy for all the species in a batchoccupancy class object
#'
#' This function takes an batchoccupancy object and predicts occupancy for all species
#' in new data, either a data.frame or a rasterstack.
#' @param batch A result from the batchoccu
#' @param new.data a rasterstack, or a dataframe containing the same variables as
#' the siteCovs variable in batchoccu
#' @return a raster stack with predictions
#' for each species.
#' @examples
#' \dontrun{
#' #Load the data
#' data("BatOccu")
#' data("Dailycov")
#' data("sampling.cov")
#' data("plumas.stack")
#'
#' #Model the abundance for 17 bat species and calculate alpha diversity from that
#'
#' BatOccupancy <-batchoccu(pres = BatOccu, sitecov = sampling.cov,
#' obscov = Dailycov,spp = 17, form = ~ Julian + Meanhum + Meantemp + sdhum +
#' sdtemp ~ Burn.intensity.soil + I(Burn.intensity.soil^2) +
#' Burn.intensity.Canopy + I(Burn.intensity.Canopy^2) + Burn.intensity.basal +
#' I(Burn.intensity.basal^2))
#'
#' Occupancy.stack <- occupancy.predict(batch = BatOccupancy, new.data =
#' plumas.stack)
#' }
#' @export
#' @seealso \code{\link[DiversityOccupancy]{batchoccu}}
#' @importFrom raster plot
#' @importFrom raster addLayer
#' @importFrom raster stack
#' @importFrom raster subset
#' @importFrom graphics plot
#' @author Derek Corcoran <derek.corcoran.barrios@gmail.com>


occupancy.predict<- function(batch, new.data) {
  models <- batch$models
  layers <- list()
  for (i in 1:length(models)){
    layers [[i]] <- predict(models[[i]], new.data, type = "state")
    layers [[i]] <- subset(layers[[i]], 1)
  }
  layers <- stack (unlist(layers))
  return(layers)
}
