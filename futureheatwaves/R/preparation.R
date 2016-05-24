#' Acquire structure of input directory
#'
#' This function parses the structure of the user-provided
#' directory of climate projection files to create a list of the
#' climate models and ensemble members included.
#'
#' @inheritParams gen_hw_set
#'
#' @return The function returns a list object outlining the file structure of
#'    the \code{dataFolder} directory.
#'
#' @note Projection, grid locations, and projection time files within the
#'    \code{dataFolder} directory must be comma-separated (.csv) files and
#'    must be named using the names specified by the arguments
#'    \code{coordinateFilenames}, \code{tasFilenames}, and \code{timeFilenames}.
#'    See the \code{futureheatwaves} vignette for more information about
#'    setting up the \code{dataFolder} directory.
#'
#' @importFrom dplyr %>%
acquireDirectoryStructure <- function(dataFolder, coordinateFilenames,
                                      tasFilenames, timeFilenames,
                                      models_to_run, dataDirectories,
                                      threshold_ensemble,
                                      thresholdBoundaries){

        # Acquire all pathnames to csv files rooted at dataPath
        all <- list.files(dataFolder, recursive = TRUE,
                          pattern = "\\.csv$")

        # Convert list to a dataframe
        split_all <- strsplit(all, "/")
        split_all <- split_all[sapply(split_all, length) == 4]
        df_all <- as.data.frame(matrix(unlist(split_all),
                                        ncol = 4, byrow = TRUE))
        colnames(df_all) <- c("exp", "model", "ens", "type")

        # Determine which experiment subdirectory to check for the
        # threshold ensemble member
        if(thresholdBoundaries[1] <= dataDirectories[[1]][2]){
                threshold_experiment <- 1
        } else {
                threshold_experiment <- 2
        }
        # Only get climate models with (a) both historical and
        # rcp85 results and (b) r1i1p1 ensemble historical results
        df_all <- dplyr::group_by_(df_all, ~ model) %>%
                dplyr::summarize_(check_1 = ~ names(dataDirectories)[1] %in% exp,
                          check_2 = ~ names(dataDirectories)[2] %in% exp,
                          check_3 = ~ paste(
                                  names(dataDirectories)[threshold_experiment],
                                          threshold_ensemble) %in%
                                  paste(exp, ens)) %>%
                dplyr::left_join(df_all, by = "model") %>%
                dplyr::filter_(~ check_1 & check_2 & check_3 &
                               type %in% c(coordinateFilenames,
                                           tasFilenames,
                                           timeFilenames)) %>%
                dplyr::select_(~ exp, ~ model, ~ ens, ~ type)

        models <- as.character(unique(df_all$model))
        experiments <- as.character(unique(df_all$exp))

        # Generate the nested lists that will be used for the processing step
        # Structure: model -> experiment -> ensemble
        finalList <- lapply(models, buildStructureModels,
                            experiments = experiments,
                            dataFolder = dataFolder,
                            coordinateFilenames = coordinateFilenames,
                            tasFilenames = tasFilenames,
                            timeFilenames = timeFilenames,
                            dataDirectories = dataDirectories)

        # Limit to only the models the user wants to run (if specified)
        if(models_to_run[1] != "all"){
                names(finalList) <- sapply(finalList, function(x) x[[1]])
                finalList <- finalList[models_to_run]
        }

        return(finalList)
}

#' Generate list of file structure
#'
#' This function takes input from \code{\link{acquireDirectoryStructure}} and
#' uses it to generate a list object with the projection directory
#' file structure. This parsed file structure is later used to lead other
#' code through all climate models and ensemble members in the input
#' projection directory.
#'
#' @param modelName Character string of climate model name (e.g., "bcc1"). This
#'    name is generated from the subdirectory name for the climate model within
#'    \code{dataFolder}.
#' @param experiments Character string of the experiment(s). Possible variables
#'    are the names of elements in the list object specified by the
#'    \code{dataDirectories} argument in \code{\link{gen_hw_set}}.
#' @inheritParams gen_hw_set
#'
#' @return A list of length 3. The first element is the name of the model
#'    whose structure was being built. The second element is, for this climate
#'    model, the hierarchy of the first subdirectory specified by
#'    \code{dataDirectories}. The third element is the hierarchy of the second
#'    subdirectory specified by \code{dataDirectories}. The second and third
#'    elements are return values of \code{\link{buildStructureExperiments}}.
buildStructureModels <- function(modelName, experiments,
                                 dataFolder,
                                 coordinateFilenames, tasFilenames,
                                 timeFilenames, dataDirectories){
        return(list(modelName,
                    buildStructureExperiments(modelName = modelName,
                                              experiment = experiments[1],
                                              dataPath = dataFolder,
                                              coordinateFilenames =
                                                      coordinateFilenames,
                                              tasFilenames = tasFilenames,
                                              timeFilenames = timeFilenames,
                                              dataDirectories = dataDirectories),
                    buildStructureExperiments(modelName = modelName,
                                              experiment = experiments[2],
                                              dataPath = dataFolder,
                                              coordinateFilenames =
                                                      coordinateFilenames,
                                              tasFilenames = tasFilenames,
                                              timeFilenames = timeFilenames,
                                              dataDirectories =
                                                      dataDirectories)))
}

#' Generate file structure for an experiment
#'
#' This function generates a list object with the file structure of files
#' in the \code{dataFolder} directory for a single experiment
#' (e.g., "historical" or "rcp85").
#'
#' @param experiment Character string of the experiment. Possible variables are
#'    the names of elements in the list object specified by the
#'    \code{dataDirectories} argument in \code{\link{gen_hw_set}}.
#' @param dataPath Character string of the file path to \code{dataFolder}.
#'    Must include the final `/`.
#' @inheritParams buildStructureModels
#' @inheritParams gen_hw_set
#'
#' @return A list with one element for each ensemble member. Each element
#'    is a return value of the \code{\link{buildStructureEnsembles}} function for
#'    one ensemble member in the experiment and climate model.
buildStructureExperiments <- function(modelName, experiment,
                                      dataPath,
                                      coordinateFilenames, tasFilenames,
                                      timeFilenames,
                                      dataDirectories){

        # List all ensembles in the given experiment
        ensembles <- list.files(paste0(dataPath, experiment, "/", modelName),
                                full.names = TRUE)

        # Build the directory structure of each ensemble
        ret <- lapply(ensembles, buildStructureEnsembles,
                      coordinateFilenames = coordinateFilenames,
                      tasFilenames = tasFilenames,
                      timeFilenames = timeFilenames,
                      dataDirectories = dataDirectories)
        return(ret)
}

#' List files for a single ensemble member
#'
#' This function reads through the user-specified \code{dataFolder} directory
#' and creates a list with pathnames to all three files (projection times, grid
#' points, and projections) for a single ensemble member.
#'
#' @param ensemblePath A character string that gives the absolute file path
#'    for files for an ensemble member within the user-specified projection
#'    directory (\code{dataFolder}).
#' @inheritParams gen_hw_set
#'
#' @return A list of length 2. The first element is the name of the ensemble
#' that was processed. The second element is a vector with the filepaths of
#' the three files.
buildStructureEnsembles <- function(ensemblePath, coordinateFilenames,
                                    tasFilenames, timeFilenames,
                                    dataDirectories){

        # Extract name of the ensemble (two directories below the
        # one named for the experiment, "historical" or "rcp85").
        splist = strsplit(ensemblePath, "/")
        ensemble_name_index <- which(sapply(splist,
                                         function(x) x %in%
                                                 names(dataDirectories)))
        # If the user has some earlier directories that match subdirectory
        # names, only take the last match (after that, should just be
        # climate models and projection files).
        if(length(ensemble_name_index) > 1){
                ensemble_name_index <- utils::tail(ensemble_name_index, 1)
        }

        ensemble_name_index <- ensemble_name_index + 2
        ensembleName <- splist[[1]][ensemble_name_index]

        # remove any irrelevant files from the file structure
        ens_files <- list.files(ensemblePath)
        coor <- ens_files[grep(coordinateFilenames, ens_files)]
        tas <- ens_files[grep(tasFilenames, ens_files)]
        time_file <- ens_files[grep(timeFilenames, ens_files)]
        ens_files <- c(coor, tas, time_file)
        directories <- unlist(lapply(ens_files, function(x){
                return(paste(ensemblePath, x, sep = "/"))
        }))
        return(c(ensembleName, directories))
}
