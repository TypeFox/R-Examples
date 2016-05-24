#' Create and write heat wave projections
#'
#' This function creates datasets of identified and characterized heat waves
#' for all ensemble members in all climate models in a directory of climate
#' projections for a user-specified set of locations.The resulting heat wave
#' projections are written out to a specified directory on the user's local
#' computer.
#'
#' @param out Character string with pathway to directory to which
#'    heat wave files will be written. This should be a pathname to a directory
#'    on the user's local computer. If the directory already exists, it will
#'    be overwritten by this function, so the user should either specify a
#'    pathname for a directory that does not yet exist or be willing to
#'    overwrite the existing directory. The parent directory of the specified
#'    directory must exist.
#' @param dataFolder Character string with pathway to a directory with
#'    climate projection data. This directory must have a specific structure--
#'    see the \code{futureheatwaves} vignette for guidance on setting up this
#'    directory.
#' @param citycsv Character string giving the filepath to a
#'    comma-separated (.csv) file with, for each study city, a unique city
#'    identifier, latitude, and longitude. These values must be specified with
#'    the column names \code{city}, \code{lat}, and \code{lon}. See the
#'    \code{futureheatwaves} vignette for guidance on setting up this
#'    file.
#' @param coordinateFilenames Character string the with filename of each
#'    grid point location file. This filename should be identical for all
#'    ensemble member subdirectories included in the \code{dataFolder} directory.
#'    See the package vignette for an example of the required structure for this
#'    file.
#' @param tasFilenames Character string the with filename of each climate
#'    projection file. This filename should be identical for all ensemble
#'    member subdirectories included in the \code{dataFolder} directory. See the
#'    package vignette for an example of the required structure for this file.
#' @param timeFilenames Character string the with filename of each projection
#'    dates file. This filename should be identical for all ensemble
#'    member  subdirectories included in the \code{dataFolder} directory. See
#'    the package vignette for an example of the required structure for this
#'    file.
#' @param IDheatwavesFunction A character string with the name of the R function
#'    to use to identify heat waves. This function may be a user-specified custom
#'    function, but it must be loaded into the current R session. The
#'    function name must be put in quotation marks. For more guidance on how to
#'    write and use a custom function to identify heat waves, see the package
#'    vignette for \code{futureheatwaves}.
#' @param thresholdBoundaries A numeric vector with the custom time boundaries
#'    to be used to determine the threshold temperatures for the heat wave
#'    definition. The required format for this vector is c(start year, end
#'    year), with the restriction that bounds must be contained within the
#'    time boundaries of one of the two experiment subdirectories specified
#'    by the \code{dataDirectories} argument in \code{\link{gen_hw_set}}.
#'    The default value is 1990 to 1999.
#' @param projectionBoundaries A numeric vector with the custom time boundaries
#'    for which the user wants to create heat wave projections. The required
#'    format for this vector is c(start year, end year), with the restriction
#'    that bounds must be contained within the time boundaries of one of the two
#'    experiment subdirectories specified by the \code{dataDirectories} argument
#'    in \code{\link{gen_hw_set}}. The default value is 2070 to 2079.
#' @param referenceBoundaries A numeric vector with the custom time boundaries
#'    to use in calculating relative characteristics for heat waves (i.e., to use
#'    when exploring the role of adaptation in projections). For more
#'    information on how reference temperatures are used, see the package
#'    vignette for \code{futureheatwaves}. The required format for this vector
#'    is c(start year, end year), with the restriction that bounds must be
#'    contained within the time boundaries of one of the two experiment
#'    subdirectories specified by the \code{dataDirectories} argument in
#'    \code{\link{gen_hw_set}}. The default value is 2070 to 2079. If the
#'    time bounds used differ from those used for projections, these reference
#'    temperatures will be pulled from the ensemble member for each climate
#'    model specified by \code{threshold_ensemble}.
#' @param probThreshold Numerical value between 0 and 1 specifying the
#'    percentile to be used for the threshold temperature used to define heat
#'    waves. The default value is 0.98 (i.e., a heat wave is a certain number of
#'    days above the city's 98th percentile temperature).
#' @param numDays Integer greater than 0 giving the number of days to
#'    use in the heat wave definition (e.g., \code{numDays = 2} would define a
#'    heat wave as two or more days above the threshold temperature).
#' @param printWarning TRUE / FALSE, specifies whether to print out a warning
#'    informing the user that the function will write out results to the local
#'    directory specified by the user with \code{out}. This warning prints out
#'    by default; the user must opt out of this warning by specifying FALSE
#'    for this argument, for example if running this function within a loop.
#' @param lat_lon_colnames A character vector of length two with the column names
#'    in the \code{citycsv} dataframe for latitude (first vector element) and
#'    longitude (second vector element)
#' @param models_to_run A character vector with either "all" (the default),
#'    in which case the function runs through all models in the \code{dataFolder}
#'    directory, or the names of the models to run, using the names of each
#'    model's subdirectory within the data directory (e.g.,
#'    \code{c("bcc1", "ccsm")}).
#' @param dataDirectories A list object, with two elements, one for each of the
#'    two subdirectories included in the main directory. Typically, these will
#'    be separate directories of historical and projection experiments from
#'    climate models. Each element of the list should be named with the name of
#'    the subdirectory and should provide a numeric vector with the starting and
#'    ending years of the data within each of the two subdirectories (e.g.,
#'    \code{list("historical" = c(1990, 1999), "rcp85" = c(2060, 2079))}
#'    for a \code{dataFolder} with historical experiment data for 1990 to 1999
#'    and RCP8.5 projections for 2060 to 2079).
#' @param threshold_ensemble A character vector giving the name of the ensemble
#'    member that should be used when determining the city-specific threshold
#'    temperatures  for each climate model (e.g., \code{"r1i1p1"}). This
#'    threshold is used for relative heat wave definitions. See the
#'    \code{futureheatwaves} vignette for more on heat wave definitions.
#'    If any climate model lacks that ensemble member for the specified
#'    dates for calculating the threshold, it will be excluded from the
#'    processing.
#' @param input_metric A character string indicating the temperature metric
#'    of the climate projection data being processed. Choices are "kelvin",
#'    "fahrenheit", and "celsius".
#'
#' @return This function creates, and writes to the user's computer, files with
#'    the heat waves and their characteristics for the specified climate
#'    projections and dates. For more information on customizing this function,
#'    see the \code{futureheatwaves} vignette. This function also returns a
#'    dataframe listing the name of each climate model processed, as well as the
#'    number of historical and future projection ensemble members for each
#'    model. This output can be used as a check that the function processed
#'    through the directory of input files specified using the \code{dataFolder}
#'    argument.
#'
#' @examples
#' \dontrun{
#' projection_dir_location <- system.file("extdata/cmip5",
#'                                       package = "futureheatwaves")
#' city_file_location <- system.file("extdata/cities.csv",
#'                                  package = "futureheatwaves")
#' gen_hw_set(out = "example_results",
#'            dataFolder = projection_dir_location ,
#'            dataDirectories = list("historical" = c(1990, 1999),
#'                                   "rcp85" = c(2060, 2079)),
#'            citycsv = city_file_location,
#'            coordinateFilenames = "latitude_longitude_NorthAmerica_12mo.csv",
#'            tasFilenames = "tas_NorthAmerica_12mo.csv",
#'            timeFilenames = "time_NorthAmerica_12mo.csv")
#' }
#'
#' @export
#'
#' @importFrom dplyr %>%
gen_hw_set <- function(out,
                       dataFolder,
                       dataDirectories = list("historical" = c(1980, 2004),
                                              "rcp85" = c(2006, 2099)),
                       citycsv,
                       coordinateFilenames,
                       tasFilenames,
                       timeFilenames,
                       IDheatwavesFunction = "IDHeatwavesCPPwrapper",
                       thresholdBoundaries = c(1990, 1999),
                       projectionBoundaries = c(2070, 2079),
                       referenceBoundaries = c(2070, 2079),
                       models_to_run = "all",
                       probThreshold = 0.98,
                       numDays = 2,
                       printWarning = TRUE,
                       threshold_ensemble = "r1i1p1",
                       lat_lon_colnames = c("lat", "lon"),
                       input_metric = "kelvin"){

        # If `dataFolder` does not end in "/", add it.
        split_dataFolder <- unlist(strsplit(dataFolder, split = ""))
        last_char <- split_dataFolder[length(split_dataFolder)]
        if(last_char != "/"){
                dataFolder <- paste0(dataFolder, "/")
        }

        # If `out` does not end in "/", add it.
        split_out <- unlist(strsplit(out, split = ""))
        last_char <- split_out[length(split_out)]
        if(last_char != "/"){
                out <- paste0(out, "/")
        }

        input_metric <- tolower(input_metric)

        # Check the parameters for errors
        check_params(out = out,
                     dataFolder = dataFolder,
                     dataDirectories = dataDirectories,
                     citycsv = citycsv,
                     coordinateFilenames = coordinateFilenames,
                     tasFilenames = tasFilenames,
                     timeFilenames = timeFilenames,
                     IDheatwavesFunction = IDheatwavesFunction,
                     thresholdBoundaries = thresholdBoundaries,
                     projectionBoundaries = projectionBoundaries,
                     referenceBoundaries = referenceBoundaries,
                     input_metric = input_metric,
                     numDays = numDays)

  # Add warning for user that this will write new files
        if(printWarning){
                cat("\n", "Warning: This function will write new files",
                    "to your computer in the \n", out,
                    "directory of your computer. If that directory already",
                    "exists,\nrunning this function will write over it. \n",
                    "Do you want to continue? (y / n): \n")
                user_prompt <- scan(n = 1, what = "character")
                user_prompt <- tolower(user_prompt)
                if(!(user_prompt %in% c("y", "yes"))){
                        stop("User chose to exit at prompt.")
                }
        }

        # Put the directories into nested list form
        models <- acquireDirectoryStructure(dataFolder = dataFolder,
                                            coordinateFilenames = coordinateFilenames,
                                            tasFilenames = tasFilenames,
                                            timeFilenames = timeFilenames,
                                            models_to_run = models_to_run,
                                            dataDirectories = dataDirectories,
                                            threshold_ensemble = threshold_ensemble,
                                            thresholdBoundaries = thresholdBoundaries)

        # Read the cities data file
        cities <- utils::read.csv(citycsv) %>%
                process_cities_file(lat_lon_colnames = lat_lon_colnames)

        # Create "global" list object that will hold variables that all
        # functions that need then will have access to
        global <- list("output" = out,
                       "data" = dataFolder,
                       "dataDirectories" = dataDirectories,
                       "cities" = cities,
                       "coordinateFilenames" = coordinateFilenames,
                       "tasFilenames" = tasFilenames,
                       "timeFilenames" = timeFilenames,
                       "threshold_ensemble" = threshold_ensemble,
                       "input_metric" = input_metric)

        # Create the "custom" list object that will hold all of the user's
        # custom settings.
        custom <- list("IDheatwaves" = IDheatwavesFunction,
                       "getBounds" = c(thresholdBoundaries,
                                       projectionBoundaries),
                       "processModel" = referenceBoundaries,
                       "createHwDataframe" = !identical(projectionBoundaries,
                                                       referenceBoundaries),
                       "probThreshold" = probThreshold,
                       "numDays" = as.integer(floor(numDays)))

        # Create accumulator closure
        accumulators <- createAccumulators()

        # Process the entire dataset
        referenceEnsembles <- sapply(models, processModel,
                                     global = global,
                                     custom = custom,
                                     accumulators = accumulators,
                                     dataDirectories = dataDirectories)

        # Write the model information from the model information accumulator
        out <- accumulators("return model information")
        writeAccumulators(modelInfoAccumulator = accumulators("return model information"),
                          locationList = accumulators("return locations"),
                          global = global)

        cat("All operations completed. Exiting.", "\n\n")
        return(out)
}

#' Check for input parameter errors
#'
#' This function goes through all parameter inputs for \code{\link{gen_hw_set}}
#' and makes sure all parameter entries are in the appropriate format.
#' If any parameters are in an incorrect format, the function stops
#' and returns an error describing the problem.
#'
#' @inheritParams gen_hw_set
#'
#' @return Only stops and returns an error if any parameters are incorrect.
#'
#' @note This function does not check if the data is organized in the proper
#'    structure or if any data exists within the directory at all, so a
#'    call to \code{\link{gen_hw_set}} could still pass through this check and
#'    make it further through the function code with those mistakes.This
#'    function also does not check if the three ensemble final .csv data
#'    files exist, only if they have the .csv extension if they do exist.
check_params <- function(out,
                         dataFolder,
                         dataDirectories,
                         citycsv,
                         coordinateFilenames,
                         tasFilenames,
                         timeFilenames,
                         IDheatwavesFunction,
                         thresholdBoundaries,
                         projectionBoundaries,
                         referenceBoundaries,
                         input_metric,
                         numDays){

        # Check to see if the folder that holds the climate data exists.
        tryCatch(
                dir.exists(dataFolder),
                error = function(){
                        stop(paste("The pathway to the projection data directory",
                                   "(`dataFolder`) is invalid."))
                },
                finally = {}
        )

        # Check if the city information .csv can be opened.
        # Note: Does not check if the city information is valid.
        tryCatch(
                utils::read.csv(citycsv, header = TRUE),
                error = function(x){
                        stop(paste("The community location file ",
                                   "(`citycsv`) cannot be opened."))
                }
        )

        # Check 'Filenames' parameters for .csv extension.
        if(!grepl(".csv", coordinateFilenames)){
                stop("The `coordinateFilenames` is not a .csv file.")
        }
        if(!grepl(".csv", tasFilenames)){
                stop("The `tasFilenames` is not a .csv file.")
        }
        if(!grepl(".csv", timeFilenames)){
                stop("The `timeFilenames` is not a .csv file.")
        }

        if(!(input_metric %in% c("kelvin", "fahrenheit", "celsius"))){
                stop("`input_metric` must be `kelvin`, `fahrenheit`, or `celsius`.")
        }

        if(numDays <= 0){
                stop("`numDays` must be 1 or larger.")
        } else if(!is.numeric(numDays)){
                stop("`numDays` must be a numeric value.")
        }

        checkCustomBounds(boundList = thresholdBoundaries,
                          dataDirectories = dataDirectories)
        checkCustomBounds(boundList = projectionBoundaries,
                          dataDirectories = dataDirectories)
        checkCustomBounds(boundList = referenceBoundaries,
                          dataDirectories = dataDirectories)
}

#' Check year boundaries for errors
#'
#' This function inputs the boundary lists specified in
#' \code{\link{gen_hw_set}}, \code{thresholdBoundaries},
#' \code{projectionBoundaries}, and
#' \code{referenceBoundaries}, and checks them for errors in structure of the
#' input or in the years selected.
#'
#' @param boundList A set of boundary years in the format
#'    c(start year, end year).
#' @inheritParams gen_hw_set
checkCustomBounds <- function(boundList, dataDirectories){

        if(class(boundList) != "numeric"){
                stop("All date boundaries must have the class `numeric`.")
        }

        if(boundList[1] > boundList[2]){
                stop(paste("In date boundaries, the first value must equal",
                           "or be lower than the second value."))
        }

        if(boundList[1] < dataDirectories[[1]][1] |
           boundList[2] > dataDirectories[[2]][2]){
                stop(paste0("Date boundaries must be within the years ",
                           dataDirectories[[1]][1],
                           " and ",
                           dataDirectories[[2]][2], ", the years you",
                           " have specified as being covered by your ",
                           " `dataDirectories`."))
        }

        if(boundList[1] <= dataDirectories[[1]][2]){
                if(boundList[2] >= dataDirectories[[2]][1]){
                        stop(paste0("Date boundaries cannot span between ",
                                   "the first (",
                                   dataDirectories[[1]][1],
                                   "-",
                                   dataDirectories[[1]][2],
                                   ") and second (",
                                   dataDirectories[[2]][1],
                                   "-",
                                   dataDirectories[[2]][2],
                                   ") directory time spans."))
                }
        }
}

#' Create accumulator closure
#'
#' This function creates a closure that holds, adds to, and returns data
#' structures that the user wishes to grow at various points in the execution of
#' the package (e.g., location and model information dataframes).
#'
#' As an example, when the generated
#' closure is used with the command "append location list", it will add
#' information on the cities and closest grid point locations based
#' on the climate model it has just completed analyzing to a growing
#' dataframe with this information for all climate models. After the function
#' run to generate the heat wave projections is completed, this closure can
#' be used with the command "return locations" to output the completed
#' dataframe of this location information.The closure can be used in a
#' similar manner to aggregate and then return meta-data on the models
#' analyzed based on their inclusion in the user-specified projections
#' directory.
#'
#' @return A closure that accepts commands to access and append new data onto
#'    data structures as the program executes. The closure created by this
#'    function accepts two arguments: (1) the command and (2) an element to
#'    be appended to the end of the data structure of the command. These two
#'    arguments must be entered in this exact order. The first argument (the
#'    command) can be one of the following options:
#'    \code{return model information}, \code{append model information},
#'    \code{return locations}, and \code{append location list}. The second
#'    argument for the created closure should only be used
#'    in conjunction with the two "append" commands for the closure.
createAccumulators <- function(){
        modelInfoAccumulator <- data.frame(c(), c(), c())
        locationList <- data.frame(c(), c(), c(), c(), c(), c())

        function(command, newElement = FALSE){

                # Commands for model information accumulator
                if(command == "return model information"){
                        return(modelInfoAccumulator)

                } else if(command == "append model information"){
                        modelInfoAccumulator <<- rbind(modelInfoAccumulator,
                                                       newElement)

                # Commands for location list accumulator
                } else if(command == "return locations"){
                        return(locationList)

                } else if(command == "append location list"){
                        locationList <<- rbind(locationList, newElement)
                }

                # If user passes an invalid command, halt the program.
                else{
                        stop("Accumulator closure: Bad command. Exiting.")
                }
        }
}
