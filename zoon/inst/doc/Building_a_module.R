## ---- echo = FALSE-------------------------------------------------------
require(zoon, quietly = TRUE)
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", fig.width=6, fig.height=6)

## ---- eval = FALSE-------------------------------------------------------
#  # Load zoon
#  library(zoon)
#  
#  # Start building our function
#  Lorem_ipsum_UK <- function(){

## ---- eval = FALSE-------------------------------------------------------
#  # First I retrieve the data from figshare
#  # Here is the URL
#  URL <- "https://ndownloader.figshare.com/files/2519918"
#  
#  # Here is the data
#  out <- read.csv(URL)
#  head(out)

## ---- echo = FALSE, message = FALSE--------------------------------------
URL <- "https://ndownloader.figshare.com/files/2519918"
head(read.csv(URL))

## ---- eval = FALSE-------------------------------------------------------
#  # Keep only Lat Long columns
#  out <- out[, c("latitude", "longitude")]
#  
#  # Add in the columns we dont have
#  out$value <- 1 # all our data are presences
#  out$type <- 'presence'
#  out$fold <- 1 # we don't add any folds
#  
#  # Now the data is in the correct format we can return it
#  return(out)

## ------------------------------------------------------------------------
Lorem_ipsum_UK <- function(){
  
  # Get data
  URL <- "https://ndownloader.figshare.com/files/2519918"
  out <- read.csv(URL)
  out <- out[, c("latitude", "longitude")]
  
  # Add in the columns we dont have
  out$value <- 1 # all our data are presences
  out$type <- 'presence'
  out$fold <- 1 # we wont add any folds
  
  return(out)
}

## ---- fig.width = 5, fig.height = 5, warning = FALSE---------------------
workl1 <- workflow(occurrence = Lorem_ipsum_UK,
                   covariate = UKBioclim,
                   process = OneHundredBackground,
                   model = LogisticRegression,
                   output = PrintMap)

## ------------------------------------------------------------------------
# Let's build our module
BuildModule(Lorem_ipsum_UK,
            type = 'occurrence',
            title = 'A dataset of Lorem ipsum occurrences',
            description = paste0('The module retrieves a dataset of',
            'Lorem ipsum records from figshare. This dataset contains',
            'precence only data and was collected between 1990 and',
            '2000 by members of to Lorem ipsum appreciation society'),
            details = 'This dataset is fake, Lorem ipsum does not exist',
            author = 'A.B. Ceidi',
            email = 'ABCD@anemail.com',
            dataType = 'presence-only')

## ---- eval = FALSE-------------------------------------------------------
#  # First we remove the function from our workspace
#  rm(list = 'Lorem_ipsum_UK')
#  
#  # This is how you would use a module that a colleague has sent you
#  LoadModule(module = 'Lorem_ipsum_UK.R')
#  
#  work2 <- workflow(occurrence = Lorem_ipsum_UK,
#                    covariate = UKBioclim,
#                    process = OneHundredBackground,
#                    model = LogisticRegression,
#                    output = PrintMap)

## ---- eval = FALSE-------------------------------------------------------
#  # Our function will take an argument to set the variable
#  # the user wants returned
#  AustraliaAir <- function(variable = 'rhum'){

## ---- echo = FALSE-------------------------------------------------------
URL <- "http://files.figshare.com/2527274/aus_air.rdata"
load(url(URL, method = 'libcurl'))
print(ras)

## ---- eval = FALSE-------------------------------------------------------
#  # Load in the data
#  URL <- "http://files.figshare.com/2527274/aus_air.rdata"
#  load(url(URL, method = 'libcurl')) # The object is called 'ras'
#  
#  # Subset the data according the the variable parameter
#  ras <- subset(ras, variables)
#  
#  return(ras)

## ---- fig.width = 5, fig.height = 5, message = FALSE, warning = FALSE----
AustraliaAir <- function(variables = 'rhum'){

  URL <- "http://files.figshare.com/2527274/aus_air.rdata"
  load(url(URL, method = 'libcurl')) # The object is called 'ras'
  ras <- subset(ras, variables)
  return(ras)
  
}

# Select the variables we want
myVariables <- c('air','hgt','rhum','shum','omega','uwnd','vwnd')

work3 <- workflow(occurrence = SpOcc(extent = c(111, 157, -46, -6),
                                     species = 'Varanus varius',
                                     limit = 500),
                  covariate = AustraliaAir(variables = myVariables),
                  process = OneHundredBackground,
                  model = LogisticRegression,
                  output = PrintMap)

## ------------------------------------------------------------------------
# Build our module
BuildModule(AustraliaAir,
            type = 'covariate',
            title = 'Australia Air data from NCEP',
            description = paste('This modules provides access to the',
                                'NCEP air data for austrlia provided by',
                                'NCEP and should be attributed to Climatic',
                                'Research Unit, University of East Anglia'),
            details = paste('These data are redistributed under the terms of',
                            'the Open Database License',
                            'http://opendatacommons.org/licenses/odbl/1.0/'),
            author = 'Z.O. Onn',
            email = 'zoon@zoon-zoon.com',
            paras = list(variables = paste('A character vector of air variables',
                         'you wish to return. This can include any number of',
                         "the following: 'air','hgt','rhum','shum','omega',",
                         "'uwnd','vwnd'")))

## ---- eval = FALSE-------------------------------------------------------
#  # remove the original function from our environment
#  rm(list = 'AustraliaAir')
#  
#  # Load the module script
#  LoadModule('AustraliaAir.R')
#  
#  work4 <- workflow(occurrence = SpOcc(extent = c(111, 157, -46, -6),
#                                       species = 'Varanus varius',
#                                       limit = 500),
#                    covariate = AustraliaAir,
#                    process = OneHundredBackground,
#                    model = LogisticRegression,
#                    output = PrintMap)

## ---- warning = FALSE----------------------------------------------------
# We run a very simple workflow so that we can get example input
# for our module
work5 <- workflow(occurrence = UKAnophelesPlumbeus,
                  covariate  = UKAir,
                  process    = NoProcess,
                  model      = LogisticRegression,
                  output     = PrintMap)

# The output from a process module is in the same format as the 
# input, so we can use the output of NoProcess as the testing
# input for our module. Note that this object should be called
# .data
.data <- work5$process.output[[1]]

str(.data, 2)

## ---- eval = FALSE-------------------------------------------------------
#  # Start writing our module
#  ClipOccurence <- function(.data, extent = c(-180, 180, -180, 180)){

## ---- eval = FALSE-------------------------------------------------------
#  # Write the body of our function
#  # extract the occurrence data from the .data object
#  occDF <- .data$df
#  
#  # Subset by longitude
#  occSub <- occDF[occDF$longitude >= extent[1] &
#                  occDF$longitude <= extent[2], ]
#  
#  # Subset by latitude
#  occSub <- occSub[occSub$latitude >= extent[3] &
#                   occSub$latitude <= extent[4], ]
#  
#  # assign this data.frame back to the .data object
#  .data$df <- occSub

## ------------------------------------------------------------------------
ClipOccurrence <- function(.data, extent = c(-180, 180, -180, 180)){
  
  # Write the body of our function
  # extract the occurrence data from the .data object
  occDF <- .data$df
  
  occSub <- occDF[occDF$longitude >= extent[1] &
                  occDF$longitude <= extent[2], ]
 
  occSub <- occSub[occSub$latitude >= extent[3] &
                   occSub$latitude <= extent[4], ]
  
  .data$df <- occSub
  
  return(.data)
  
}

## ---- message = FALSE, warning = FALSE-----------------------------------
# Run a workflow with our new process
# In this example we first add background points, then clip the data
work6 <- workflow(occurrence = UKAnophelesPlumbeus,
                  covariate  = UKAir,
                  process    = Chain(OneHundredBackground,
                                     ClipOccurrence(extent = c(-3, 2, 50, 53))),
                  model      = LogisticRegression,
                  output     = PrintMap)

## ------------------------------------------------------------------------
# Build our module
BuildModule(ClipOccurrence,
            type = 'process',
            title = 'Clip occurrence data to extent',
            description = paste('This process module clips the occurrence',
                                'data that is returned from the occurrence',
                                'module to a user defined extent'),
            details = paste('The extent is a square region which denotes the',
                            'area within which observations will be kept.',
                            'All data that falls outside of the extent will',
                            'be removed and will be not be used in the',
                            'modelling process'),
            author = 'Z.O. Onn',
            email = 'zoon@zoon-zoon.com',
            paras = list(extent = paste('A numeric vector of length for',
                                        'giving (in this order) the minimum',
                                        'longitude, maximum longitude, minimum',
                                        'latitude, maximum latitude.')),
            dataType = c('presence-only', 'presence/absence',
                         'presence/background', 'abundance',
                         'proportion'))

## ---- warning = FALSE----------------------------------------------------
# remove the original function from our environment
rm(list = 'ClipOccurrence')

# Load the module script
LoadModule('ClipOccurrence.R')

work7 <- workflow(occurrence = CWBZimbabwe,
                  covariate = Bioclim(extent = c(31, 34, -22, -18)),
                  process = ClipOccurrence(extent = c(32, 33, -21, -19)),
                  model = LogisticRegression,
                  output = PrintMap)

## ---- eval = FALSE-------------------------------------------------------
#  GamGam <- function(.df){

## ---- eval = FALSE-------------------------------------------------------
#  # Specify the packages we need using the function
#  # GetPackage
#  zoon::GetPackage("gam")

## ---- eval = FALSE-------------------------------------------------------
#  # Create a data.frame of covariate data
#  covs <- as.data.frame(.df[, 6:ncol(.df)])
#  names(covs) <- names(.df)[6:ncol(.df)]
#  
#  # do a bit of copy-pasting to define smooth terms for each covariate
#  f <- sprintf('.df$value ~ s(%s)',
#                      paste(colnames(covs),
#                            collapse = ') + s('))
#  
#  # Run our gam model
#  m <- gam::gam(formula = formula(f),
#                data = covs,
#                family = binomial)

## ---- eval = FALSE-------------------------------------------------------
#  # Create a ZoonModel object to return.
#  # this includes our model, predict method
#  # and the packages we need.
#  ZoonModel(model = m,
#            code = {
#  
#            # create empty vector of predictions
#            p <- rep(NA, nrow(newdata))
#  
#            # omit NAs in new data
#            newdata_clean <- na.omit(newdata)
#  
#            # get NA indices
#            na_idx <- attr(newdata_clean, 'na.action')
#  
#            # if there are no NAs then the index should
#            # include all rows, else it should name the
#            # rows to ignore
#            if (is.null(na_idx)){
#              idx <- 1:nrow(newdata)
#            } else {
#              idx <- -na_idx
#            }
#  
#            # Use the predict function in gam to predict
#            # our new values
#            p[idx] <- gam::predict.gam(model,
#                                       newdata_clean,
#                                       type = 'response')
#            return (p)
#          },
#          packages = 'gam')

## ------------------------------------------------------------------------
GamGam <- function(.df){

  # Specify the packages we need using the function
  # GetPackage
  zoon::GetPackage("gam")
  
  # Create a data.frame of covariate data
  covs <- as.data.frame(.df[, 6:ncol(.df)])
  names(covs) <- names(.df)[6:ncol(.df)]
  
  # do a bit of copy-pasting to define smooth terms for each covariate
  f <- sprintf('.df$value ~ s(%s)',
                      paste(colnames(covs),
                            collapse = ') + s('))
  
  # Run our gam model
  m <- gam::gam(formula = formula(f),
                data = covs,
                family = binomial)
  
  # Create a ZoonModel object to return.
  # this includes our model, predict method
  # and the packages we need.
  ZoonModel(model = m,
            code = {
            
            # create empty vector of predictions
            p <- rep(NA, nrow(newdata))
            
            # omit NAs in new data
            newdata_clean <- na.omit(newdata)
            
            # get their indices
            na_idx <- attr(newdata_clean, 'na.action')
            
            # if there are no NAs then the index should 
            # include all rows, else it should name the 
            # rows to ignore
            if (is.null(na_idx)){
              idx <- 1:nrow(newdata)
            } else {
              idx <- -na_idx
            }
            
            # Use the predict function in gam to predict
            # our new values
            p[idx] <- gam::predict.gam(model,
                                       newdata_clean,
                                       type = 'response')
            return (p)
          },
          packages = 'gam')
  
}

## ----BuildMod------------------------------------------------------------
BuildModule(object = GamGam,
            type = 'model',
            title = 'GAM sdm model',
            description = 'This is my mega cool new model.',
            details = paste('This module performs GAMs (Generalised Additive',
                            'Models) using the gam function from the package gam.'),
            author = 'Z. Oon',
            email = 'zoon@zoon.com',
            dataType = c('presence-only', 'presence/absence'))

## ---- warning=FALSE, message=FALSE---------------------------------------
# remove the function in our workspace else
# this will cause problems
rm(GamGam)

# Load in teh module we just built
LoadModule('GamGam.R')

# Run a workflow using our module
work8 <- workflow(occurrence = UKAnophelesPlumbeus,
                  covariate = UKAir,
                  process  = OneHundredBackground,
                  model = GamGam,
                  output = PrintMap)

## ---- warning = FALSE, eval = FALSE--------------------------------------
#  # We run a very simple workflow so that we can get example input
#  # for our module
#  work9 <- workflow(occurrence = UKAnophelesPlumbeus,
#                    covariate  = UKAir,
#                    process    = OneHundredBackground,
#                    model      = LogisticRegression,
#                    output     = PrintMap)
#  
#  # The input to an output module is a combination of the output
#  # from the model module and the covariate module. We can recreate
#  # it for this work flow like this
#  .model <- work9$model.output[[1]]
#  .ras <- work9$covariate.output[[1]]

## ---- eval = FALSE-------------------------------------------------------
#  # Our output module takes the default parameters and a user-defined
#  # Raster* object that has the same structure as the raster layer output
#  # by the covariate module
#  PredictNewRasterMap <- function(.model, .ras, raster = .ras){

## ---- eval = FALSE-------------------------------------------------------
#  # The first step is to load in the packages we need
#  zoon::GetPackage(raster)
#  
#  # Then extract the covariate values
#  # from the user provided raster
#  vals <- data.frame(getValues(raster))
#  colnames(vals) <- names(raster)

## ---- eval = FALSE-------------------------------------------------------
#  # Make predictions to the new values
#  pred <- ZoonPredict(.model$model,
#                      newdata = vals)
#  
#  # Create a copy of the users' raster...
#  # (just a single layer)
#  pred_ras <- raster[[1]]
#  
#  # ... and assign the predicted values to it
#  pred_ras <- setValues(pred_ras, pred)

## ---- eval = FALSE-------------------------------------------------------
#  # Plot the predictions as a map
#  plot(pred_ras)
#  
#  # Return the raster of predictions
#  return (pred_ras)

## ------------------------------------------------------------------------
PredictNewRasterMap <- function(.model, .ras, raster = .ras){
  
  zoon::GetPackage(raster)
  
  # Extract the values from the user provided raster
  vals <- data.frame(getValues(raster))
  colnames(vals) <- names(raster)
  
  # Make predictions to the new values
  pred <- ZoonPredict(.model$model,
                      newdata = vals)
  
  pred_ras <- raster[[1]]
  pred_ras <- setValues(pred_ras, pred)
  
  # Print the predictions as a map
  plot(pred_ras)
  
  return(pred_ras)
}

## ---- message = FALSE, warning = FALSE-----------------------------------
# Run it with the defaults
work10 <- workflow(occurrence = UKAnophelesPlumbeus,
                   covariate  = UKBioclim,
                   process    = OneHundredBackground,
                   model      = LogisticRegression,
                   output     = PredictNewRasterMap)

# Now I'm going to run it with a different raster
library(raster)

# Get Bioclim data (using the getData function in the raster package,
# which zoon loads) ...
BioclimData <- getData('worldclim', var = 'bio', res = 5)
BioclimData <- BioclimData[[1:19]]

# ... and crop to Australia
cropped <- crop(BioclimData,
                c(109,155,-46,-7))

# Run it with my new raster
work11 <- workflow(occurrence = UKAnophelesPlumbeus,
                   covariate  = UKBioclim,
                   process    = OneHundredBackground,
                   model      = LogisticRegression,
                   output     = PredictNewRasterMap(raster = cropped))

# The prediction map should also be returned as a raster
str(work11$report, 2)

## ------------------------------------------------------------------------
# Build our module
BuildModule(PredictNewRasterMap,
            type = 'output',
            title = 'Predict to a new raster and map',
            description = paste('This output module predicts the species',
                                'distribution in a new area given a new',
                                'raster'),
            details = paste('The results are printed as a map and a raster is',
                            'returned with the predicted values. It is important',
                            'that the new raster has the same structure as the',
                            'raster provided by the covariate module.',
                            'It must have the same covariate columns in the',
                            'same order.'),
            author = 'Z.O. On',
            email = 'zoon@zoon-zoon.com',
            paras = list(raster = paste('A RasterBrick, RasterLayer or RasterStack in',
                                        'the same format as the raster provided',
                                        'by the covariate module. Predicted values',
                                        'will be estimated for this raster using',
                                        'the results from the model module')),
            dataType = c('presence-only', 'presence/absence', 'abundance',
                         'proportion'))

## ---- warning = FALSE----------------------------------------------------
# remove the original function from our environment
rm(list = 'PredictNewRasterMap')

# Load the module script
LoadModule('PredictNewRasterMap.R')

# Now I model a crop pest from Zimbabwe in its home
# range and in Australia by chaining together
# output modules
work12 <- workflow(occurrence = CWBZimbabwe,
                   covariate = Bioclim(extent = c(28, 38, -24, -16)),
                   process = NoProcess,
                   model = RandomForest,
                   output = Chain(PrintMap,
                                  PredictNewRasterMap(raster = cropped)))

## ---- echo = FALSE-------------------------------------------------------
# Clean up
unlink('AustraliaAir.R')
unlink('ClipOccurrence.R')
unlink('GamGam.R')
unlink('Lorem_ipsum_UK.R')
unlink('PredictNewRasterMap.R')

