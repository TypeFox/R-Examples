
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2015-04-03 12:24:30 +0200 (Fri, 03 Apr 2015) $
# $Rev: 338 $
# ----------------------------------------------------------------

## ----------------------------------------------------------------
## Process of reading in climate data and storing them to a data.frame
##  object will take place here.
##
## - models2wux, the 'main' WUX data input routine.
## - Helperfunctions for models2wux.
## ----------------------------------------------------------------


models2wux <- function(userinput,
                       modelinput = NULL
                       ) {
  ## Creates a dataframe containing climate change signals of climate models
  ## listed in user.input.
  ##
  ## Args:
  ##   userinput: Filename ot list object from input created by the user.
  ##                  The file contains a
  ##                   general.input section with tuning parameters and
  ##                   a model.input section listing the climate models
  ##                   to be processed.
  ##   modelinput: filepath to model.input file. default is NULL,
  ##                which reads in the internal IniModelsDictionary file.
  ##
  ## Returns:
  ##   Dataframe containing climate change signals for each model.
  ##
  ## History:
  ##   2010-10-20 | Original code (thm)
  ##   2010-11-03 | user.input sourced from Models2Wux
  ##                model.input inserted in user.input
  ##                renaming of several variables (thm)
  ##   2010-11-17 | model.info read in after reading in parameter data (thm)
  ##   2010-11-29 | some additional user.input - model.input keywords added
  ##   2011-02-11 | user.input handling changes, renaming from 'Modes2Wux'
  ##                to 'models2wux'
  ##   2011-04-29 | structural changes for generating WUX dataframe
  ##   2011-09-14 | added display.subregions keyword to models2wux function
  ##   2011-09-22 | added case filenames = NA (thus not reading in any data)
  ##              | plus: return value added (it got lost somewhere...)
  ##   2014-11-12 | alternative model.input pathway introduced (thm)
  ##   2014-11-14 | rename argument "input" to "userinput" (thm)
  ##   2014-11-19 | now omitting argument does.plot.subregions (thm)
  ##   2014-11-21 | platform independent (thm)
  ##   2015-03-31 | new handling of subregion plotting

  cat("\n")

  ## get user.input list from 'input.filename'
  if (!is.list(userinput)) {
    user.input.sourced <- source(userinput, local = TRUE)
    user.input <- user.input.sourced$value
  } else {
    user.input <- userinput
  }

  ## -------------------------------------------
  ## ++++++++++ READING GENERAL INPUT ++++++++++
  ## -------------------------------------------

############## EXTRACTING INFORMATION FROM user.input FILE ##############

  ## create outdir for data if path not existing
  save.as.data <- user.input$save.as.data
  if ( !file.exists(dirname(save.as.data)) )
    dir.create(dirname(save.as.data), recursive = TRUE)

  ## create outdir for plots if path not existing
  if ( !is.null(user.input$plot.subregions$save.subregions.plots) ) {
    save.subregions.plots <- user.input$plot.subregions$save.subregions.plots
    if ( !file.exists(save.subregions.plots) )
      dir.create(save.subregions.plots, recursive = TRUE)
  } else {
    save.subregions.plots <- NULL
  }

  ## user input for area fraction
  if ( !is.null(user.input$area.fraction) ) {
    does.area.fraction <- user.input$area.fraction
  } else {
    does.area.fraction <- FALSE
  }

  ## user input for land mask
  if ( !is.null(user.input$use.land.mask) ) {
    use.land.mask <- user.input$use.land.mask
    print("USING LAND-SEA MASK")
  } else {
    use.land.mask <- FALSE
  }

  ##  user input for area weighting
  if ( !is.null(user.input$spatial.weighting) ) {
    does.spatial.weighting <- user.input$spatial.weighting
  } else {
    does.spatial.weighting <- FALSE
  }

  reference.period <- user.input$reference.period
  scenario.period <- user.input$scenario.period
  periods <- c("reference.period" = reference.period,
               "scenario.period" = scenario.period)

  ##  user input for area weighting
  if ( !is.null(user.input$na.rm) ) {
    na.rm <- user.input$na.rm
  } else {
    na.rm <- FALSE
  }

  ## subregion specification
  subregions <- user.input$subregions

  ## temporal aggregation specification
  temporal.aggr <- user.input$temporal.aggregation

  ## generate list for dataframes for different parameters
  parameter.names <- user.input$parameter.names
  wux.df.list <- vector("list", length(parameter.names))
  names(wux.df.list) <- parameter.names

  ## -------------------------------------------
  ## +++++++++++ READING MODEL DATA ++++++++++++
  ## -------------------------------------------

############## LOOP OVER PARMS, MODELS, PERIODS AND EXTRACT DATA ##############

  ## first loop over specified parameters
  for ( parameter.longname in parameter.names ) {

    cat("\n******************************************************************")
    cat("****\n", "PROCESSING PARAMETER \"", parameter.longname,
        "\"\n", sep = "")
    cat("******************************************************************")
    cat("****\n")

    ## initialize WUX dataset
    wux.df <- data.frame()

    ## climate model names from user.input
    user.model.names <- user.input$climate.models
    ## load Model Dictionary, i.e. a named list containing all the information
    ## for known models as filenames, directories and further metainformation
    model.input <- ReadFromModelDictionary(user.model.names, modelinput)
    model.names <- names(model.input)

    ## second loop over specified models
    for (modelname.counter in model.names) {
      ## Need parameter short name (i.e. variable name from NetCDF file)
      ## to extract data.
      ## If "parameters" tag is provided in InitModelDictionary, use its value
      ## as NetCDF variable name.
      model.pars <- model.input[[modelname.counter]]$parameters
      parameter.shortname <- GetParameterShortName(parameter.longname,
                                                   model.pars)

      ## getting general model info for data to be read in
      institute.name <- model.input[[modelname.counter]]$institute
      global.model.name <- model.input[[modelname.counter]]$gcm
      gcm.run <- model.input[[modelname.counter]]$gcm.run
      rcm.name <- model.input[[modelname.counter]]$rcm
      emission.scenario <- model.input[[modelname.counter]]$emission.scenario
      is.corrected <- model.input[[modelname.counter]]$corrected.data
      resolution <-  model.input[[modelname.counter]]$resolution
      what.timesteps <- model.input[[modelname.counter]]$what.timesteps
      ## optional input for models having other periods
      model.ref.per <- model.input[[modelname.counter]]$reference.period
      model.scn.per <- model.input[[modelname.counter]]$scenario.period
      ## optional input for model having different period lengths per file
      years.in.ref.per.files <-
        model.input[[modelname.counter]]$years.in.ref.per.files
      years.in.scn.per.files <-
        model.input[[modelname.counter]]$years.in.scn.per.files
      ## optional input for respecification of NetCDF entries
      time.units <- model.input[[modelname.counter]]$time.units
      calendar <- model.input[[modelname.counter]]$calendar
      count.first.time.value <-
        model.input[[modelname.counter]]$count.first.time.value
      lonlat.var.name <-
        model.input[[modelname.counter]]$lonlat.var.name

      ## location of gridfile
      ## model grid
      grid.filename <-
        file.path(model.input[[modelname.counter]]$gridfile.path,
              model.input[[modelname.counter]]$gridfile.filename)
      ## model subregion
      subregion.filename <-
        file.path(model.input[[modelname.counter]]$subregion.path,
              model.input[[modelname.counter]]$subregion.filename)
      ## in case no subregionfiles are passed, we flag filename with "FALSE"
      if (length(subregion.filename) == 0)
        subregion.filename <- FALSE
      ## land mask
      land.mask <- model.input[[modelname.counter]]$land.mask
      land.mask.name <- model.input[[modelname.counter]]$land.mask.name

      ## third loop over periods (reference period and scenario period)
      for (periods.counter in periods) {
        cat("\n-------------------------------------------------------------")
        cat("---------\n\n  PROCESSING MODEL \"",  modelname.counter,
            "\" FOR PERIOD \"", periods.counter, "\"\n", sep = "")

        ## convert periods to POSIXct dates
        period.begin <- ISOdate(strsplit(periods.counter, "-")[[1]][1],
                                1, 1, 0)
        period.end <-   ISOdate(strsplit(periods.counter, "-")[[1]][2],
                                12, 31, 23)

        ## flag whether in reference period or not
        current.period.char = names(which(periods == periods.counter))
        if (current.period.char == "reference.period") {
          is.refperiod <- "yes"
          ## in case model has specific reference period, use it here
          if (!is.null(model.ref.per))
            periods.counter <- model.ref.per
        } else {
          is.refperiod <- "no"
          ## in case model has specific scenario period, use it here
          if (!is.null(model.scn.per))
            periods.counter <- model.scn.per
        }

        ## retrieving filenames whose data will be read in and processed
        ## get required files to cover the period
        cat("    GETTING REQUIRED FILENAMES FOR PERIOD", "\n", sep = "")
        current.model.input <- model.input[[modelname.counter]]
        filenames <- GetFileNames(current.model.input,
                                  period.begin, period.end,
                                  parameter.shortname,
                                  period = current.period.char,
                                  calendar = calendar,
                                  count.first.time.value =
                                  count.first.time.value,
                                  time.units = time.units,
                                  what.timesteps = what.timesteps)
        cat(paste("      ", filenames, "\n", sep=""))

        ## read in data, corresponding grids and aggregate over time to
        ## get corresponding seasonal values
        if ( !all(is.na(filenames)) ) {

          data.aggr.list <-
          GetAggregatedSpatioTemporalData(filenames,
                                          model.name = modelname.counter,
                                          plot.subregion = user.input$plot.subregion,
                                          grid.filenames = grid.filename,
                                          subregions = subregions,
                                          parameter.name = parameter.shortname,
                                          interpolate.to.grid = FALSE,
                                          temporal.aggr = temporal.aggr,
                                          startdate = period.begin,
                                          enddate = period.end,
                                          ## '...' for  ReadNetCdfTimeData
                                          time.units = time.units,
                                          calendar = calendar,
                                          count.first.time.value =
                                          count.first.time.value,
                                          lonlat.var.name = lonlat.var.name,
                                          ## for AggregateTemporal(...)
                                          what.timesteps = what.timesteps,
                                          na.rm = na.rm,
                                          ## '...' for GetSubregionShapes
                                          area.fraction = does.area.fraction,
                                          spatial.weighting = does.spatial.weighting,
                                          ## '...' for GetAgg***Data.nointerpol
                                          use.land.mask = use.land.mask,
                                          land.mask = land.mask,
                                          land.mask.name = land.mask.name)

          ## get number of datapoints read in each subregion (not used further)
          n.of.pixels <-
            lapply(data.aggr.list, function(z) {
              lapply(z, function(x) {lapply(x, function(y) length(!is.na(y)))})})

          ## aggregate over spatial component
          data.aggr.list <- lapply(data.aggr.list,
                                         function(x)
                                         lapply(x$data,
                                                AggregateWeightedSubregions,
                                                weight=x$weight)
                                         )

          ## check if dealing with a time series
          if (sum(which(unlist(sapply(temporal.aggr,
                                      "[", "time.series")) == FALSE)) == 0) {
            is.time.series <- TRUE
          } else {
            is.time.series <- FALSE
          }

          ## get dataframe for particular model & period & parameter
          model.data <- MakeModelDataFrame(modelname.counter,
                                           institute.name,
                                           rcm.name,
                                           global.model.name,
                                           gcm.run = gcm.run,
                                           is.reference.period = is.refperiod,
                                           time.step = periods.counter,
                                           subregion = names(data.aggr.list),
                                           is.corrected = is.corrected,
                                           resolution = resolution,
                                           season = sort(names(data.aggr.list[[1]])),
                                           is.time.series = is.time.series,
                                           data = data.aggr.list,
                                           parameter.name = parameter.longname,
                                           emission.scenario = emission.scenario)

          ## append to dataset of previos models
          wux.df <- rbind(wux.df, model.data)
        }
      }
    }

    ## now we have all models read in for this specific parameter, so we plug
    ## the data.frame to a list to merged to a dataframe again later
    wux.df.list[[parameter.longname]] <- wux.df
  }

############## MERGING ALL CREATED data.frames AND SAVE THEM ON DISK  ##########

  ## merge datasets with different parmers to single data.frame
  wux.df <- MergeParameterDataframes(wux.df.list, parameter.names)

  ## getting climate change signal calculating the diff of scen.per and ref.per.
  if (!is.time.series) {
    ## wux.df has still seperated data for ref.per and scen.per. To get the
    ## climate change signal we calculate the difference of these two periods
    ## (and for precipitation the perentual change as well)
    wux.diff.df <- WuxDataDiff(wux.df, parameter.names)
    write.table(wux.diff.df,
                file = paste(save.as.data, "_diff.csv", sep = ""), sep=";")
  }

  ## storing WUX dataset to harddisk
  cat("\n-------------------------------------------------------------------")
  cat("---\n")
  cat("\n", "SAVING DATA", "\n", sep = "")
  write.table(wux.df, file = paste(save.as.data, ".csv", sep = ""), sep=";")

  ## time to say goodbye
  cat("\n*******************************************************************")
  cat("***\n")
  cat("         YOUR WUX DATA FRAME HAS BEEN CREATED SUCCESSFULLY :-)      ")
  cat("    \n")
  cat("********************************************************************")
  cat("**\n")

  if (!is.time.series) {
    class(wux.df) <- append("wux.df", class(wux.df))
    return(wux.diff.df)
  } else {
    class(wux.df) <- append("wux.ts.df", class(wux.df))
    return(wux.df)
  }

}


ReadFromModelDictionary <- function(modelnames, modelinput = NULL){
  ## Extracts model information from InitModelDictionary for given models.
  ##
  ## Args:
  ##   modelnames: Character vector of model names (acronyms) to be extracted
  ##               from InitModelDictionary, or model ensemble name for multiple
  ##               climate models, which then has to be defined here
  ##               (eg. "ensembles").
  ##   modelinput: filepath to model.input file. default is NULL,
  ##                which reads in the internal IniModelsDictionary file.
  ##
  ## Returns:
  ##   Named list of model information needed for models2wux.
  ##
  ## History:
  ##   2010-10-27 | Original code (thm)
  ##       ...    | many things have happened
  ##   2012-08-29 | extended CMIP5 ensemble, alphabetized vectors (msu)
  ##   2014-04-14 | extended CMIP5 ensemble. up to date ensemble! (thm)
  ##   2014-11-12 | alternative model.input possible (thm)

  ## check wheteher user specified a specific model input (NULL means no)
  if (is.null(modelinput)){
    ## read in whole InitModelDictionary model.list
    model.input <- InitModelDictionary()
  } else {
    cat(is.list(modelinput))
    ## check whether modelinput is a directory or a list object
    if (is.list(modelinput)) {
      ## this is a list object
      model.input <- modelinput
    } else {
      ## this is a file to be read in
      model.input.sourced <- source(modelinput, local = TRUE)
      model.input <- model.input.sourced$value
    }
    if (!is.list(model.input))
      stop("YOUR model.input FILE IS WRONG. BE SURE YOU HAVE SPECIFIED A LIST CALLED \"model.input\". SEE ?models2wux FOR DETAILS. ")
  }

  ## known model ensembles from the CMIP3 project and the "ENSEMBLES" project
  ##-------------------------------------------------------------------
  ## CMIP3 A1B scenario:
  cmip3.sresa1b.labels <- c("cmip3-sresa1b", "CMIP3-SRESA1B")
  cmip3.sresa1b.modelnames <- c("bccr_bcm2_0-r1",
                                "cccma_cgcm3_1-r1",
                                "cccma_cgcm3_1-r2",
                                "cccma_cgcm3_1-r3",
                                "cccma_cgcm3_1-r4",
                                "cccma_cgcm3_1-r5",
                                "cccma_cgcm3_1_t63-r1",
                                "cnrm_cm3-r1",
                                "csiro_mk3_0-r1",
                                "csiro_mk3_5-r1",
                                "gfdl_cm2_0-r1",
                                "gfdl_cm2_1-r1",
                                "giss_aom-r1",
                                "giss_aom-r2",
                                "giss_model_e_h-r1",
                                "giss_model_e_h-r2",
                                "giss_model_e_h-r3",
                                "giss_model_e_r-r1",
                                "giss_model_e_r-r2",
                                "giss_model_e_r-r3",
                                "giss_model_e_r-r4",
                                "giss_model_e_r-r5",
                                "iap_fgoals1_0_g-r1",
                                "iap_fgoals1_0_g-r2",
                                "iap_fgoals1_0_g-r3",
                                "ingv_echam4-r1",
                                "inmcm3_0-r1",
                                "ipsl_cm4-r1",
                                "miroc3_2_hires-r1",
                                "miroc3_2_medres-r1",
                                "miroc3_2_medres-r2",
                                "miroc3_2_medres-r3",
                                "miub_echo_g-r1",
                                "miub_echo_g-r2",
                                "miub_echo_g-r3",
                                "mpi_echam5-r1",
                                "mpi_echam5-r2",
                                "mpi_echam5-r3",
                                "mpi_echam5-r4",
                                "mri_cgcm2_3_2a-r1",
                                "mri_cgcm2_3_2a-r2",
                                "mri_cgcm2_3_2a-r3",
                                "mri_cgcm2_3_2a-r4",
                                "mri_cgcm2_3_2a-r5",
                                "ncar_ccsm3_0-r1",
                                "ncar_ccsm3_0-r2",
                                "ncar_ccsm3_0-r3",
                                "ncar_ccsm3_0-r5",
                                "ncar_ccsm3_0-r6",
                                "ncar_ccsm3_0-r7",
                                "ncar_ccsm3_0-r9",
                                "ukmo_hadcm3-r1",
                                "ukmo_hadgem1-r1")
  ##-------------------------------------------------------------------


  ##-------------------------------------------------------------------
  ## CMIP3 A2 scenario:
  cmip3.sresa2.labels <- c("cmip3-sresa2", "CMIP3-SRESA2")
  cmip3.sresa2.modelnames <- c("bccr_bcm2_0-sresa2-r1",
                               "cccma_cgcm3_1-sresa2-r1",
                               "cccma_cgcm3_1-sresa2-r2",
                               "cccma_cgcm3_1-sresa2-r3",
                               "cccma_cgcm3_1-sresa2-r4",
                               "cccma_cgcm3_1-sresa2-r5",
                               "cnrm_cm3-sresa2-r1",
                               "csiro_mk3_0-sresa2-r1",
                               "csiro_mk3_5-sresa2-r1",
                               "gfdl_cm2_0-sresa2-r1",
                               "gfdl_cm2_1-sresa2-r1",
                               "giss_model_e_r-sresa2-r1",
                               "ingv_echam4-sresa2-r1",
                               "inmcm3_0-sresa2-r1",
                               "ipsl_cm4-sresa2-r1",
                               "miroc3_2_medres-sresa2-r1",
                               "miroc3_2_medres-sresa2-r2",
                               "miroc3_2_medres-sresa2-r3",
                               "miub_echo_g-sresa2-r1",
                               "miub_echo_g-sresa2-r2",
                               "miub_echo_g-sresa2-r3",
                               "mpi_echam5-sresa2-r1",
                               "mpi_echam5-sresa2-r2",
                               "mpi_echam5-sresa2-r3",
                               "mri_cgcm2_3_2a-sresa2-r1",
                               "mri_cgcm2_3_2a-sresa2-r2",
                               "mri_cgcm2_3_2a-sresa2-r3",
                               "mri_cgcm2_3_2a-sresa2-r4",
                               "mri_cgcm2_3_2a-sresa2-r5",
                               "ncar_ccsm3_0-sresa2-r1",
                               "ncar_ccsm3_0-sresa2-r2",
                               "ncar_ccsm3_0-sresa2-r3",
                               "ncar_ccsm3_0-sresa2-r4",
                               "ncar_ccsm3_0-sresa2-r5",
                               "ukmo_hadcm3-sresa2-r1",
                               "ukmo_hadgem1-sresa2-r1")
  ##-------------------------------------------------------------------


  ##-------------------------------------------------------------------
  ## CMIP3 B1 scenario:
  cmip3.sresb1.labels <- c("cmip3-sresb1", "CMIP3-SRESB1")
  cmip3.sresb1.modelnames <- c("bccr_bcm2_0-sresb1-r1",
                               "cccma_cgcm3_1-sresb1-r1",
                               "cccma_cgcm3_1-sresb1-r2",
                               "cccma_cgcm3_1-sresb1-r3",
                               "cccma_cgcm3_1-sresb1-r4",
                               "cccma_cgcm3_1-sresb1-r5",
                               "cccma_cgcm3_1_t63-sresb1-r1",
                               "cnrm_cm3-sresb1-r1",
                               "csiro_mk3_0-sresb1-r1",
                               "csiro_mk3_5-sresb1-r1",
                               "gfdl_cm2_0-sresb1-r1",
                               "gfdl_cm2_1-sresb1-r1",
                               "giss_aom-sresb1-r1",
                               "giss_aom-sresb1-r2",
                               "giss_model_e_r-sresb1-r1" ,
                               "iap_fgoals1_0_g-sresb1-r1",
                               "iap_fgoals1_0_g-sresb1-r2",
                               "iap_fgoals1_0_g-sresb1-r3",
                               "inmcm3_0-sresb1-r1",
                               "ipsl_cm4-sresb1-r1",
                               "miroc3_2_hires-sresb1-r1",
                               "miroc3_2_medres-sresb1-r1",
                               "miroc3_2_medres-sresb1-r2",
                               "miroc3_2_medres-sresb1-r3",
                               "miub_echo_g-sresb1-r1",
                               "miub_echo_g-sresb1-r2",
                               "miub_echo_g-sresb1-r3",
                               "mpi_echam5-sresb1-r1",
                               "mpi_echam5-sresb1-r2",
                               "mpi_echam5-sresb1-r3",
                               "mri_cgcm2_3_2a-sresb1-r1",
                               "mri_cgcm2_3_2a-sresb1-r2",
                               "mri_cgcm2_3_2a-sresb1-r3",
                               "mri_cgcm2_3_2a-sresb1-r4",
                               "mri_cgcm2_3_2a-sresb1-r5",
                               "ncar_ccsm3_0-sresb1-r1",
                               "ncar_ccsm3_0-sresb1-r2",
                               "ncar_ccsm3_0-sresb1-r3",
                               "ncar_ccsm3_0-sresb1-r4",
                               "ncar_ccsm3_0-sresb1-r5",
                               "ncar_ccsm3_0-sresb1-r6",
                               "ncar_ccsm3_0-sresb1-r7",
                               "ncar_ccsm3_0-sresb1-r9",
                               "ukmo_hadcm3-sresb1-r1")
  ##-------------------------------------------------------------------


  ##-------------------------------------------------------------------
  ## CMIP5
  ## rcp26
  cmip5.rcp26.labels <- c("cmip5-rcp26", "CMIP5-RCP26")
  cmip5.rcp26.modelnames <- c("BCC-CSM1-1-r1i1p1_rcp26", 
                              "BCC-CSM1-1-m-r1i1p1_rcp26", 
                              "BNU-ESM-r1i1p1_rcp26", 
                              "CanESM2-r1i1p1_rcp26", 
                              "CanESM2-r2i1p1_rcp26", 
                              "CanESM2-r3i1p1_rcp26", 
                              "CanESM2-r4i1p1_rcp26", 
                              "CanESM2-r5i1p1_rcp26", 
                              "CCSM4-r1i1p1_rcp26", 
                              "CCSM4-r2i1p1_rcp26", 
                              "CCSM4-r3i1p1_rcp26", 
                              "CCSM4-r4i1p1_rcp26", 
                              "CCSM4-r5i1p1_rcp26", 
                              "CCSM4-r6i1p1_rcp26", 
                              "CESM1-CAM5-r1i1p1_rcp26", 
                              "CESM1-CAM5-r2i1p1_rcp26", 
                              "CESM1-CAM5-r3i1p1_rcp26", 
                              "CNRM-CM5-r1i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r10i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r1i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r2i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r3i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r4i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r5i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r6i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r7i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r8i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r9i1p1_rcp26", 
                              "EC-EARTH-r12i1p1_rcp26", 
                              "EC-EARTH-r8i1p1_rcp26", 
                              "FGOALS-g2-r1i1p1_rcp26", 
                              "FGOALS-s2-r1i1p1_rcp26", 
                              "FIO-ESM-r1i1p1_rcp26", 
                              "FIO-ESM-r2i1p1_rcp26", 
                              "FIO-ESM-r3i1p1_rcp26", 
                              "GFDL-CM3-r1i1p1_rcp26", 
                              "GFDL-ESM2G-r1i1p1_rcp26", 
                              "GFDL-ESM2M-r1i1p1_rcp26", 
                              "GISS-E2-H-r1i1p1_rcp26", 
                              "GISS-E2-H-r1i1p2_rcp26", 
                              "GISS-E2-H-r1i1p3_rcp26", 
                              "GISS-E2-R-r1i1p1_rcp26", 
                              "GISS-E2-R-r1i1p2_rcp26", 
                              "GISS-E2-R-r1i1p3_rcp26", 
                              "HadGEM2-AO-r1i1p1_rcp26", 
                              "HadGEM2-ES-r1i1p1_rcp26", 
                              ## "HadGEM2-ES-r2i1p1_rcp26", making touble with dates in NetCDF file
                              ## "HadGEM2-ES-r3i1p1_rcp26",  making touble with dates in NetCDF file
                              "HadGEM2-ES-r4i1p1_rcp26", 
                              "IPSL-CM5A-LR-r1i1p1_rcp26", 
                              "IPSL-CM5A-LR-r2i1p1_rcp26", 
                              "IPSL-CM5A-LR-r3i1p1_rcp26", 
                              "IPSL-CM5A-LR-r4i1p1_rcp26", 
                              "IPSL-CM5A-MR-r1i1p1_rcp26", 
                              "MIROC5-r1i1p1_rcp26", 
                              "MIROC5-r2i1p1_rcp26", 
                              "MIROC5-r3i1p1_rcp26", 
                              "MIROC-ESM-r1i1p1_rcp26", 
                              "MIROC-ESM-CHEM-r1i1p1_rcp26", 
                              "MPI-ESM-LR-r1i1p1_rcp26", 
                              "MPI-ESM-LR-r2i1p1_rcp26", 
                              "MPI-ESM-LR-r3i1p1_rcp26", 
                              "MPI-ESM-MR-r1i1p1_rcp26", 
                              "MRI-CGCM3-r1i1p1_rcp26", 
                              "NorESM1-M-r1i1p1_rcp26", 
                              "NorESM1-ME-r1i1p1_rcp26")
  ##-------------------------------------------------------------------
  ## rcp26 including models with shorter simulations (usually until 2035)
  ## with MIROC5_r4i1p1_rcp26
  ##      MIROC5_r5i1p1_rcp26
  cmip5.rcp26.shortperiods.labels <- c("cmip5-rcp26-shortperiods", "CMIP5-RCP26-SHORTPERIODS")
  cmip5.rcp26.shortperiods.modelnames <- c("BCC-CSM1-1-r1i1p1_rcp26", 
                              "BCC-CSM1-1-m-r1i1p1_rcp26", 
                              "BNU-ESM-r1i1p1_rcp26", 
                              "CanESM2-r1i1p1_rcp26", 
                              "CanESM2-r2i1p1_rcp26", 
                              "CanESM2-r3i1p1_rcp26", 
                              "CanESM2-r4i1p1_rcp26", 
                              "CanESM2-r5i1p1_rcp26", 
                              "CCSM4-r1i1p1_rcp26", 
                              "CCSM4-r2i1p1_rcp26", 
                              "CCSM4-r3i1p1_rcp26", 
                              "CCSM4-r4i1p1_rcp26", 
                              "CCSM4-r5i1p1_rcp26", 
                              "CCSM4-r6i1p1_rcp26", 
                              "CESM1-CAM5-r1i1p1_rcp26", 
                              "CESM1-CAM5-r2i1p1_rcp26", 
                              "CESM1-CAM5-r3i1p1_rcp26", 
                              "CNRM-CM5-r1i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r10i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r1i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r2i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r3i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r4i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r5i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r6i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r7i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r8i1p1_rcp26", 
                              "CSIRO-Mk3-6-0-r9i1p1_rcp26", 
                              "EC-EARTH-r12i1p1_rcp26", 
                              "EC-EARTH-r8i1p1_rcp26", 
                              "FGOALS-g2-r1i1p1_rcp26", 
                              "FGOALS-s2-r1i1p1_rcp26", 
                              "FIO-ESM-r1i1p1_rcp26", 
                              "FIO-ESM-r2i1p1_rcp26", 
                              "FIO-ESM-r3i1p1_rcp26", 
                              "GFDL-CM3-r1i1p1_rcp26", 
                              "GFDL-ESM2G-r1i1p1_rcp26", 
                              "GFDL-ESM2M-r1i1p1_rcp26", 
                              "GISS-E2-H-r1i1p1_rcp26", 
                              "GISS-E2-H-r1i1p2_rcp26", 
                              "GISS-E2-H-r1i1p3_rcp26", 
                              "GISS-E2-R-r1i1p1_rcp26", 
                              "GISS-E2-R-r1i1p2_rcp26", 
                              "GISS-E2-R-r1i1p3_rcp26", 
                              "HadGEM2-AO-r1i1p1_rcp26", 
                              "HadGEM2-ES-r1i1p1_rcp26", 
                              "HadGEM2-ES-r2i1p1_rcp26", 
                              "HadGEM2-ES-r3i1p1_rcp26", 
                              "HadGEM2-ES-r4i1p1_rcp26", 
                              "IPSL-CM5A-LR-r1i1p1_rcp26", 
                              "IPSL-CM5A-LR-r2i1p1_rcp26", 
                              "IPSL-CM5A-LR-r3i1p1_rcp26", 
                              "IPSL-CM5A-LR-r4i1p1_rcp26", 
                              "IPSL-CM5A-MR-r1i1p1_rcp26", 
                              "MIROC5-r1i1p1_rcp26", 
                              "MIROC5-r2i1p1_rcp26", 
                              "MIROC5-r3i1p1_rcp26", 
                              "MIROC5-r4i1p1_rcp26", 
                              "MIROC5-r5i1p1_rcp26", 
                              "MIROC-ESM-r1i1p1_rcp26", 
                              "MIROC-ESM-CHEM-r1i1p1_rcp26", 
                              "MPI-ESM-LR-r1i1p1_rcp26", 
                              "MPI-ESM-LR-r2i1p1_rcp26", 
                              "MPI-ESM-LR-r3i1p1_rcp26", 
                              "MPI-ESM-MR-r1i1p1_rcp26", 
                              "MRI-CGCM3-r1i1p1_rcp26", 
                              "NorESM1-M-r1i1p1_rcp26", 
                              "NorESM1-ME-r1i1p1_rcp26")
  ##-------------------------------------------------------------------


  ##-------------------------------------------------------------------
  ## rcp45
  cmip5.rcp45.labels <- c("cmip5-rcp45", "CMIP5-RCP45")
  cmip5.rcp45.modelnames <- c("ACCESS1-0-r1i1p1_rcp45", 
                              "ACCESS1-3-r1i1p1_rcp45", 
                              "BCC-CSM1-1-r1i1p1_rcp45", 
                              "BCC-CSM1-1-m-r1i1p1_rcp45", 
                              "BNU-ESM-r1i1p1_rcp45", 
                              ## "CanCM4-r10i1p1_rcp45", too short timeseries...
                              ## "CanCM4-r1i1p1_rcp45", 
                              ## "CanCM4-r2i1p1_rcp45", 
                              ## "CanCM4-r3i1p1_rcp45", 
                              ## "CanCM4-r4i1p1_rcp45", 
                              ## "CanCM4-r5i1p1_rcp45", 
                              ## "CanCM4-r6i1p1_rcp45", 
                              ## "CanCM4-r7i1p1_rcp45", 
                              ## "CanCM4-r8i1p1_rcp45", 
                              ## "CanCM4-r9i1p1_rcp45", 
                              "CanESM2-r1i1p1_rcp45", 
                              "CanESM2-r2i1p1_rcp45", 
                              "CanESM2-r3i1p1_rcp45", 
                              "CanESM2-r4i1p1_rcp45", 
                              "CanESM2-r5i1p1_rcp45", 
                              "CCSM4-r1i1p1_rcp45", 
                              "CCSM4-r2i1p1_rcp45", 
                              "CCSM4-r3i1p1_rcp45", 
                              "CCSM4-r4i1p1_rcp45", 
                              "CCSM4-r5i1p1_rcp45", 
                              "CCSM4-r6i1p1_rcp45", 
                              "CESM1-BGC-r1i1p1_rcp45", 
                              "CESM1-CAM5-r1i1p1_rcp45", 
                              "CESM1-CAM5-r2i1p1_rcp45", 
                              "CESM1-CAM5-r3i1p1_rcp45", 
                              "CESM1-WACCM-r2i1p1_rcp45", 
                              ## "CESM1-WACCM-r3i1p1_rcp45",   too short timeseries...
                              ## "CESM1-WACCM-r4i1p1_rcp45", 
                              "CMCC-CM-r1i1p1_rcp45", 
                              "CMCC-CMS-r1i1p1_rcp45", 
                              "CNRM-CM5-r1i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r10i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r1i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r2i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r3i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r4i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r5i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r6i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r7i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r8i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r9i1p1_rcp45", 
                              "EC-EARTH-r12i1p1_rcp45", 
                              "EC-EARTH-r13i1p1_rcp45", 
                              "EC-EARTH-r14i1p1_rcp45", 
                              ## "EC-EARTH-r1i1p1_rcp45",  too short timeseries...
                              "EC-EARTH-r2i1p1_rcp45", 
                              "EC-EARTH-r6i1p1_rcp45", 
                              ## "EC-EARTH-r7i1p1_rcp45", ## missing historical runs
                              "EC-EARTH-r8i1p1_rcp45", 
                              "EC-EARTH-r9i1p1_rcp45", 
                              "FGOALS-g2-r1i1p1_rcp45", 
                              "FIO-ESM-r1i1p1_rcp45", 
                              "FIO-ESM-r2i1p1_rcp45", 
                              "FIO-ESM-r3i1p1_rcp45", 
                              ## "GFDL-CM2-1-r10i1p1_rcp45", too short timeseries...
                              ## "GFDL-CM2-1-r1i1p1_rcp45", 
                              ## "GFDL-CM2-1-r2i1p1_rcp45", 
                              ## "GFDL-CM2-1-r3i1p1_rcp45", 
                              ## "GFDL-CM2-1-r4i1p1_rcp45", 
                              ## "GFDL-CM2-1-r5i1p1_rcp45", 
                              ## "GFDL-CM2-1-r6i1p1_rcp45", 
                              ## "GFDL-CM2-1-r7i1p1_rcp45", 
                              ## "GFDL-CM2-1-r8i1p1_rcp45", 
                              ## "GFDL-CM2-1-r9i1p1_rcp45", 
                              "GFDL-CM3-r1i1p1_rcp45", 
                              "GFDL-ESM2G-r1i1p1_rcp45", 
                              "GFDL-ESM2M-r1i1p1_rcp45", 
                              "GISS-E2-H-r1i1p1_rcp45", 
                              "GISS-E2-H-r1i1p2_rcp45", 
                              "GISS-E2-H-r1i1p3_rcp45", 
                              "GISS-E2-H-r2i1p1_rcp45", 
                              "GISS-E2-H-r2i1p2_rcp45", 
                              "GISS-E2-H-r2i1p3_rcp45", 
                              "GISS-E2-H-r3i1p1_rcp45", 
                              "GISS-E2-H-r3i1p2_rcp45", 
                              "GISS-E2-H-r3i1p3_rcp45", 
                              "GISS-E2-H-r4i1p1_rcp45", 
                              "GISS-E2-H-r4i1p2_rcp45", 
                              "GISS-E2-H-r4i1p3_rcp45", 
                              "GISS-E2-H-r5i1p1_rcp45", 
                              "GISS-E2-H-r5i1p2_rcp45", 
                              "GISS-E2-H-r5i1p3_rcp45", 
                              "GISS-E2-H-r6i1p3_rcp45", 
                              "GISS-E2-H-CC-r1i1p1_rcp45", 
                              "GISS-E2-R-r1i1p1_rcp45", 
                              "GISS-E2-R-r1i1p2_rcp45", 
                              "GISS-E2-R-r1i1p3_rcp45", 
                              "GISS-E2-R-r2i1p1_rcp45", 
                              "GISS-E2-R-r2i1p2_rcp45", 
                              "GISS-E2-R-r2i1p3_rcp45", 
                              "GISS-E2-R-r3i1p1_rcp45", 
                              "GISS-E2-R-r3i1p2_rcp45", 
                              "GISS-E2-R-r3i1p3_rcp45", 
                              "GISS-E2-R-r4i1p1_rcp45", 
                              "GISS-E2-R-r4i1p2_rcp45", 
                              "GISS-E2-R-r4i1p3_rcp45", 
                              "GISS-E2-R-r5i1p1_rcp45", 
                              "GISS-E2-R-r5i1p2_rcp45", 
                              "GISS-E2-R-r5i1p3_rcp45", 
                              "GISS-E2-R-r6i1p1_rcp45", 
                              "GISS-E2-R-r6i1p3_rcp45", 
                              "GISS-E2-R-CC-r1i1p1_rcp45", 
                              ## "HadCM3-r10i1p1_rcp45", too short timeseries...
                              ## "HadCM3-r1i1p1_rcp45", 
                              ## "HadCM3-r2i1p1_rcp45", 
                              ## "HadCM3-r3i1p1_rcp45", 
                              ## "HadCM3-r4i1p1_rcp45", 
                              ## "HadCM3-r5i1p1_rcp45", 
                              ## "HadCM3-r6i1p1_rcp45", 
                              ## "HadCM3-r7i1p1_rcp45", 
                              ## "HadCM3-r8i1p1_rcp45", 
                              ## "HadCM3-r9i1p1_rcp45", 
                              "HadGEM2-AO-r1i1p1_rcp45", 
                              "HadGEM2-CC-r1i1p1_rcp45", 
                              "HadGEM2-ES-r1i1p1_rcp45", 
                              ## "HadGEM2-ES-r2i1p1_rcp45", making touble with dates in NetCDF file
                              ## "HadGEM2-ES-r3i1p1_rcp45", making touble with dates in NetCDF file
                              "HadGEM2-ES-r4i1p1_rcp45", 
                              "INM-CM4-r1i1p1_rcp45", 
                              "IPSL-CM5A-LR-r1i1p1_rcp45", 
                              "IPSL-CM5A-LR-r2i1p1_rcp45", 
                              "IPSL-CM5A-LR-r3i1p1_rcp45", 
                              "IPSL-CM5A-LR-r4i1p1_rcp45", 
                              "IPSL-CM5A-MR-r1i1p1_rcp45", 
                              "IPSL-CM5B-LR-r1i1p1_rcp45", 
                              ## "MIROC4h-r1i1p1_rcp45", too short timeseries...
                              ## "MIROC4h-r2i1p1_rcp45", 
                              ## "MIROC4h-r3i1p1_rcp45", 
                              "MIROC5-r1i1p1_rcp45", 
                              "MIROC5-r2i1p1_rcp45", 
                              "MIROC5-r3i1p1_rcp45", 
                              ## "MIROC5-r4i1p1_rcp45", too short timeseries...
                              ## "MIROC5-r5i1p1_rcp45", 
                              "MIROC-ESM-r1i1p1_rcp45", 
                              "MIROC-ESM-CHEM-r1i1p1_rcp45", 
                              "MPI-ESM-LR-r1i1p1_rcp45", 
                              "MPI-ESM-LR-r2i1p1_rcp45", 
                              "MPI-ESM-LR-r3i1p1_rcp45", 
                              "MPI-ESM-MR-r1i1p1_rcp45", 
                              "MPI-ESM-MR-r2i1p1_rcp45", 
                              "MPI-ESM-MR-r3i1p1_rcp45", 
                              "MRI-CGCM3-r1i1p1_rcp45", 
                              "NorESM1-M-r1i1p1_rcp45", 
                              "NorESM1-ME-r1i1p1_rcp45"
                              ## "CSIRO-Mk3L-1-2-r1i2p1_rcp45", ## missing historical runs
                              ## "CSIRO-Mk3L-1-2-r2i2p1_rcp45", 
                              ## "CSIRO-Mk3L-1-2-r3i2p1_rcp45",
                              )
  ## rcp45
  cmip5.rcp45.shortperiods.labels <- c("cmip5-rcp45-shortperiods", "CMIP5-RCP45-SHORTPERIODS")
  cmip5.rcp45.shortperiods.modelnames <- c("ACCESS1-0-r1i1p1_rcp45", 
                              "ACCESS1-3-r1i1p1_rcp45", 
                              "BCC-CSM1-1-r1i1p1_rcp45", 
                              "BCC-CSM1-1-m-r1i1p1_rcp45", 
                              "BNU-ESM-r1i1p1_rcp45", 
                              "CanCM4-r10i1p1_rcp45", 
                              "CanCM4-r1i1p1_rcp45", 
                              "CanCM4-r2i1p1_rcp45", 
                              "CanCM4-r3i1p1_rcp45", 
                              "CanCM4-r4i1p1_rcp45", 
                              "CanCM4-r5i1p1_rcp45", 
                              "CanCM4-r6i1p1_rcp45", 
                              "CanCM4-r7i1p1_rcp45", 
                              "CanCM4-r8i1p1_rcp45", 
                              "CanCM4-r9i1p1_rcp45", 
                              "CanESM2-r1i1p1_rcp45", 
                              "CanESM2-r2i1p1_rcp45", 
                              "CanESM2-r3i1p1_rcp45", 
                              "CanESM2-r4i1p1_rcp45", 
                              "CanESM2-r5i1p1_rcp45", 
                              "CCSM4-r1i1p1_rcp45", 
                              "CCSM4-r2i1p1_rcp45", 
                              "CCSM4-r3i1p1_rcp45", 
                              "CCSM4-r4i1p1_rcp45", 
                              "CCSM4-r5i1p1_rcp45", 
                              "CCSM4-r6i1p1_rcp45", 
                              "CESM1-BGC-r1i1p1_rcp45", 
                              "CESM1-CAM5-r1i1p1_rcp45", 
                              "CESM1-CAM5-r2i1p1_rcp45", 
                              "CESM1-CAM5-r3i1p1_rcp45", 
                              "CESM1-WACCM-r2i1p1_rcp45", 
                              "CESM1-WACCM-r3i1p1_rcp45", 
                              "CESM1-WACCM-r4i1p1_rcp45", 
                              "CMCC-CM-r1i1p1_rcp45", 
                              "CMCC-CMS-r1i1p1_rcp45", 
                              "CNRM-CM5-r1i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r10i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r1i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r2i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r3i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r4i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r5i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r6i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r7i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r8i1p1_rcp45", 
                              "CSIRO-Mk3-6-0-r9i1p1_rcp45", 
                              "EC-EARTH-r12i1p1_rcp45", 
                              "EC-EARTH-r13i1p1_rcp45", 
                              "EC-EARTH-r14i1p1_rcp45", 
                              "EC-EARTH-r1i1p1_rcp45", 
                              "EC-EARTH-r2i1p1_rcp45", 
                              "EC-EARTH-r6i1p1_rcp45", 
                              ## "EC-EARTH-r7i1p1_rcp45", ## missing historical runs
                              "EC-EARTH-r8i1p1_rcp45", 
                              "EC-EARTH-r9i1p1_rcp45", 
                              "FGOALS-g2-r1i1p1_rcp45", 
                              "FIO-ESM-r1i1p1_rcp45", 
                              "FIO-ESM-r2i1p1_rcp45", 
                              "FIO-ESM-r3i1p1_rcp45", 
                              "GFDL-CM2-1-r10i1p1_rcp45", 
                              "GFDL-CM2-1-r1i1p1_rcp45", 
                              "GFDL-CM2-1-r2i1p1_rcp45", 
                              "GFDL-CM2-1-r3i1p1_rcp45", 
                              "GFDL-CM2-1-r4i1p1_rcp45", 
                              "GFDL-CM2-1-r5i1p1_rcp45", 
                              "GFDL-CM2-1-r6i1p1_rcp45", 
                              "GFDL-CM2-1-r7i1p1_rcp45", 
                              "GFDL-CM2-1-r8i1p1_rcp45", 
                              "GFDL-CM2-1-r9i1p1_rcp45", 
                              "GFDL-CM3-r1i1p1_rcp45", 
                              "GFDL-ESM2G-r1i1p1_rcp45", 
                              "GFDL-ESM2M-r1i1p1_rcp45", 
                              "GISS-E2-H-r1i1p1_rcp45", 
                              "GISS-E2-H-r1i1p2_rcp45", 
                              "GISS-E2-H-r1i1p3_rcp45", 
                              "GISS-E2-H-r2i1p1_rcp45", 
                              "GISS-E2-H-r2i1p2_rcp45", 
                              "GISS-E2-H-r2i1p3_rcp45", 
                              "GISS-E2-H-r3i1p1_rcp45", 
                              "GISS-E2-H-r3i1p2_rcp45", 
                              "GISS-E2-H-r3i1p3_rcp45", 
                              "GISS-E2-H-r4i1p1_rcp45", 
                              "GISS-E2-H-r4i1p2_rcp45", 
                              "GISS-E2-H-r4i1p3_rcp45", 
                              "GISS-E2-H-r5i1p1_rcp45", 
                              "GISS-E2-H-r5i1p2_rcp45", 
                              "GISS-E2-H-r5i1p3_rcp45", 
                              "GISS-E2-H-r6i1p3_rcp45", 
                              "GISS-E2-H-CC-r1i1p1_rcp45", 
                              "GISS-E2-R-r1i1p1_rcp45", 
                              "GISS-E2-R-r1i1p2_rcp45", 
                              "GISS-E2-R-r1i1p3_rcp45", 
                              "GISS-E2-R-r2i1p1_rcp45", 
                              "GISS-E2-R-r2i1p2_rcp45", 
                              "GISS-E2-R-r2i1p3_rcp45", 
                              "GISS-E2-R-r3i1p1_rcp45", 
                              "GISS-E2-R-r3i1p2_rcp45", 
                              "GISS-E2-R-r3i1p3_rcp45", 
                              "GISS-E2-R-r4i1p1_rcp45", 
                              "GISS-E2-R-r4i1p2_rcp45", 
                              "GISS-E2-R-r4i1p3_rcp45", 
                              "GISS-E2-R-r5i1p1_rcp45", 
                              "GISS-E2-R-r5i1p2_rcp45", 
                              "GISS-E2-R-r5i1p3_rcp45", 
                              "GISS-E2-R-r6i1p1_rcp45", 
                              "GISS-E2-R-r6i1p3_rcp45", 
                              "GISS-E2-R-CC-r1i1p1_rcp45", 
                              "HadCM3-r10i1p1_rcp45", 
                              "HadCM3-r1i1p1_rcp45", 
                              "HadCM3-r2i1p1_rcp45", 
                              "HadCM3-r3i1p1_rcp45", 
                              "HadCM3-r4i1p1_rcp45", 
                              "HadCM3-r5i1p1_rcp45", 
                              "HadCM3-r6i1p1_rcp45", 
                              "HadCM3-r7i1p1_rcp45", 
                              "HadCM3-r8i1p1_rcp45", 
                              "HadCM3-r9i1p1_rcp45", 
                              "HadGEM2-AO-r1i1p1_rcp45", 
                              "HadGEM2-CC-r1i1p1_rcp45", 
                              "HadGEM2-ES-r1i1p1_rcp45", 
                              "HadGEM2-ES-r2i1p1_rcp45", 
                              "HadGEM2-ES-r3i1p1_rcp45", 
                              "HadGEM2-ES-r4i1p1_rcp45", 
                              "INM-CM4-r1i1p1_rcp45", 
                              "IPSL-CM5A-LR-r1i1p1_rcp45", 
                              "IPSL-CM5A-LR-r2i1p1_rcp45", 
                              "IPSL-CM5A-LR-r3i1p1_rcp45", 
                              "IPSL-CM5A-LR-r4i1p1_rcp45", 
                              "IPSL-CM5A-MR-r1i1p1_rcp45", 
                              "IPSL-CM5B-LR-r1i1p1_rcp45", 
                              "MIROC4h-r1i1p1_rcp45", 
                              "MIROC4h-r2i1p1_rcp45", 
                              "MIROC4h-r3i1p1_rcp45", 
                              "MIROC5-r1i1p1_rcp45", 
                              "MIROC5-r2i1p1_rcp45", 
                              "MIROC5-r3i1p1_rcp45", 
                              "MIROC5-r4i1p1_rcp45", 
                              "MIROC5-r5i1p1_rcp45", 
                              "MIROC-ESM-r1i1p1_rcp45", 
                              "MIROC-ESM-CHEM-r1i1p1_rcp45", 
                              "MPI-ESM-LR-r1i1p1_rcp45", 
                              "MPI-ESM-LR-r2i1p1_rcp45", 
                              "MPI-ESM-LR-r3i1p1_rcp45", 
                              "MPI-ESM-MR-r1i1p1_rcp45", 
                              "MPI-ESM-MR-r2i1p1_rcp45", 
                              "MPI-ESM-MR-r3i1p1_rcp45", 
                              "MRI-CGCM3-r1i1p1_rcp45", 
                              "NorESM1-M-r1i1p1_rcp45", 
                              "NorESM1-ME-r1i1p1_rcp45"
                              ## "CSIRO-Mk3L-1-2-r1i2p1_rcp45", ## missing historical runs
                              ## "CSIRO-Mk3L-1-2-r2i2p1_rcp45", 
                              ## "CSIRO-Mk3L-1-2-r3i2p1_rcp45",
                              )
  ##-------------------------------------------------------------------

  ##-------------------------------------------------------------------
  ## rcp60
  cmip5.rcp60.labels <- c("cmip5-rcp60", "CMIP5-RCP60")
  cmip5.rcp60.modelnames <- c("BCC-CSM1-1-r1i1p1_rcp60", 
                              "BCC-CSM1-1-m-r1i1p1_rcp60", 
                              "CCSM4-r1i1p1_rcp60", 
                              "CCSM4-r2i1p1_rcp60", 
                              "CCSM4-r3i1p1_rcp60", 
                              "CCSM4-r4i1p1_rcp60", 
                              "CCSM4-r5i1p1_rcp60", 
                              "CCSM4-r6i1p1_rcp60", 
                              "CESM1-CAM5-r1i1p1_rcp60", 
                              "CESM1-CAM5-r2i1p1_rcp60", 
                              "CESM1-CAM5-r3i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r10i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r1i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r2i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r3i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r4i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r5i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r6i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r7i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r8i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r9i1p1_rcp60", 
                              "FGOALS-s2-r1i1p1_rcp60", 
                              "FIO-ESM-r1i1p1_rcp60", 
                              "FIO-ESM-r2i1p1_rcp60", 
                              "FIO-ESM-r3i1p1_rcp60", 
                              "GFDL-CM3-r1i1p1_rcp60", 
                              "GFDL-ESM2G-r1i1p1_rcp60", 
                              "GFDL-ESM2M-r1i1p1_rcp60", 
                              "GISS-E2-H-r1i1p1_rcp60", 
                              "GISS-E2-H-r1i1p2_rcp60", 
                              "GISS-E2-H-r1i1p3_rcp60", 
                              "GISS-E2-R-r1i1p1_rcp60", 
                              "GISS-E2-R-r1i1p2_rcp60", 
                              "GISS-E2-R-r1i1p3_rcp60", 
                              "HadGEM2-AO-r1i1p1_rcp60", 
                              "HadGEM2-ES-r1i1p1_rcp60", 
                              ## "HadGEM2-ES-r2i1p1_rcp60",     making touble with dates in NetCDF file
                              ## "HadGEM2-ES-r3i1p1_rcp60",  making touble with dates in NetCDF file
                              "HadGEM2-ES-r4i1p1_rcp60", 
                              "IPSL-CM5A-LR-r1i1p1_rcp60", 
                              "IPSL-CM5A-MR-r1i1p1_rcp60", 
                              "MIROC5-r1i1p1_rcp60", 
                              "MIROC5-r2i1p1_rcp60", 
                              "MIROC5-r3i1p1_rcp60", 
                              ## "MIROC5-r4i1p1_rcp60",  too short timeseries...
                              ## "MIROC5-r5i1p1_rcp60", 
                              "MIROC-ESM-r1i1p1_rcp60", 
                              "MIROC-ESM-CHEM-r1i1p1_rcp60", 
                              "MRI-CGCM3-r1i1p1_rcp60", 
                              "NorESM1-M-r1i1p1_rcp60", 
                              "NorESM1-ME-r1i1p1_rcp60")
  ##-------------------------------------------------------------------
  ##-------------------------------------------------------------------
  ## rcp60 including all siulations with shorter periods
  cmip5.rcp60.shortperiods.labels <- c("cmip5-rcp60-shortperiods", "CMIP5-RCP60-SHORTPERIODS")
  cmip5.rcp60.shortperiods.modelnames <- c("BCC-CSM1-1-r1i1p1_rcp60", 
                              "BCC-CSM1-1-m-r1i1p1_rcp60", 
                              "CCSM4-r1i1p1_rcp60", 
                              "CCSM4-r2i1p1_rcp60", 
                              "CCSM4-r3i1p1_rcp60", 
                              "CCSM4-r4i1p1_rcp60", 
                              "CCSM4-r5i1p1_rcp60", 
                              "CCSM4-r6i1p1_rcp60", 
                              "CESM1-CAM5-r1i1p1_rcp60", 
                              "CESM1-CAM5-r2i1p1_rcp60", 
                              "CESM1-CAM5-r3i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r10i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r1i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r2i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r3i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r4i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r5i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r6i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r7i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r8i1p1_rcp60", 
                              "CSIRO-Mk3-6-0-r9i1p1_rcp60", 
                              "FGOALS-s2-r1i1p1_rcp60", 
                              "FIO-ESM-r1i1p1_rcp60", 
                              "FIO-ESM-r2i1p1_rcp60", 
                              "FIO-ESM-r3i1p1_rcp60", 
                              "GFDL-CM3-r1i1p1_rcp60", 
                              "GFDL-ESM2G-r1i1p1_rcp60", 
                              "GFDL-ESM2M-r1i1p1_rcp60", 
                              "GISS-E2-H-r1i1p1_rcp60", 
                              "GISS-E2-H-r1i1p2_rcp60", 
                              "GISS-E2-H-r1i1p3_rcp60", 
                              "GISS-E2-R-r1i1p1_rcp60", 
                              "GISS-E2-R-r1i1p2_rcp60", 
                              "GISS-E2-R-r1i1p3_rcp60", 
                              "HadGEM2-AO-r1i1p1_rcp60", 
                              "HadGEM2-ES-r1i1p1_rcp60", 
                              "HadGEM2-ES-r2i1p1_rcp60", 
                              "HadGEM2-ES-r3i1p1_rcp60", 
                              "HadGEM2-ES-r4i1p1_rcp60", 
                              "IPSL-CM5A-LR-r1i1p1_rcp60", 
                              "IPSL-CM5A-MR-r1i1p1_rcp60", 
                              "MIROC5-r1i1p1_rcp60", 
                              "MIROC5-r2i1p1_rcp60", 
                              "MIROC5-r3i1p1_rcp60", 
                              "MIROC5-r4i1p1_rcp60", 
                              "MIROC5-r5i1p1_rcp60", 
                              "MIROC-ESM-r1i1p1_rcp60", 
                              "MIROC-ESM-CHEM-r1i1p1_rcp60", 
                              "MRI-CGCM3-r1i1p1_rcp60", 
                              "NorESM1-M-r1i1p1_rcp60", 
                              "NorESM1-ME-r1i1p1_rcp60")
  ##-------------------------------------------------------------------


  ##-------------------------------------------------------------------
  ## rcp85
  cmip5.rcp85.labels <- c("cmip5-rcp85", "CMIP5-RCP85")
  cmip5.rcp85.modelnames <- c("ACCESS1-0-r1i1p1_rcp85", 
                              "ACCESS1-3-r1i1p1_rcp85", 
                              "BCC-CSM1-1-r1i1p1_rcp85", 
                              "BCC-CSM1-1-m-r1i1p1_rcp85", 
                              "BNU-ESM-r1i1p1_rcp85", 
                              "CanESM2-r1i1p1_rcp85", 
                              "CanESM2-r2i1p1_rcp85", 
                              "CanESM2-r3i1p1_rcp85", 
                              "CanESM2-r4i1p1_rcp85", 
                              "CanESM2-r5i1p1_rcp85", 
                              "CCSM4-r1i1p1_rcp85", 
                              "CCSM4-r2i1p1_rcp85", 
                              "CCSM4-r3i1p1_rcp85", 
                              "CCSM4-r4i1p1_rcp85", 
                              "CCSM4-r5i1p1_rcp85", 
                              "CCSM4-r6i1p1_rcp85", 
                              "CESM1-BGC-r1i1p1_rcp85", 
                              "CESM1-CAM5-r1i1p1_rcp85", 
                              "CESM1-CAM5-r2i1p1_rcp85", 
                              "CESM1-CAM5-r3i1p1_rcp85", 
                              "CESM1-WACCM-r2i1p1_rcp85", 
                              ## "CESM1-WACCM-r3i1p1_rcp85",too short timeseries
                              ## "CESM1-WACCM-r4i1p1_rcp85", 
                              "CMCC-CESM-r1i1p1_rcp85", 
                              "CMCC-CM-r1i1p1_rcp85", 
                              "CMCC-CMS-r1i1p1_rcp85", 
                              "CNRM-CM5-r10i1p1_rcp85", 
                              "CNRM-CM5-r1i1p1_rcp85", 
                              "CNRM-CM5-r2i1p1_rcp85", 
                              "CNRM-CM5-r4i1p1_rcp85", 
                              "CNRM-CM5-r6i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r10i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r1i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r2i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r3i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r4i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r5i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r6i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r7i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r8i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r9i1p1_rcp85", 
                              "EC-EARTH-r12i1p1_rcp85", 
                              "EC-EARTH-r13i1p1_rcp85", 
                              "EC-EARTH-r14i1p1_rcp85", 
                              ## "EC-EARTH-r1i1p1_rcp85", ## corrupt scen tas files
                              "EC-EARTH-r2i1p1_rcp85", 
                              "EC-EARTH-r6i1p1_rcp85", 
                              "EC-EARTH-r8i1p1_rcp85", 
                              "EC-EARTH-r9i1p1_rcp85", 
                              "FGOALS-g2-r1i1p1_rcp85", 
                              "FIO-ESM-r1i1p1_rcp85", 
                              "FIO-ESM-r2i1p1_rcp85", 
                              "FIO-ESM-r3i1p1_rcp85", 
                              "GFDL-CM3-r1i1p1_rcp85", 
                              "GFDL-ESM2G-r1i1p1_rcp85", 
                              "GFDL-ESM2M-r1i1p1_rcp85", 
                              "GISS-E2-H-r1i1p1_rcp85", 
                              "GISS-E2-H-r1i1p2_rcp85", 
                              "GISS-E2-H-r1i1p3_rcp85", 
                              "GISS-E2-H-r2i1p1_rcp85", 
                              "GISS-E2-H-r2i1p3_rcp85", 
                              "GISS-E2-H-CC-r1i1p1_rcp85", 
                              "GISS-E2-R-r1i1p1_rcp85", 
                              "GISS-E2-R-r1i1p2_rcp85", 
                              "GISS-E2-R-r1i1p3_rcp85", 
                              "GISS-E2-R-r2i1p1_rcp85", 
                              "GISS-E2-R-r2i1p3_rcp85", 
                              "GISS-E2-R-CC-r1i1p1_rcp85", 
                              "HadGEM2-AO-r1i1p1_rcp85", 
                              "HadGEM2-CC-r1i1p1_rcp85", 
                              ## "HadGEM2-CC-r2i1p1_rcp85",  making touble with dates in NetCDF file
                              ##  "HadGEM2-CC-r3i1p1_rcp85", making touble with dates in NetCDF file
                              "HadGEM2-ES-r1i1p1_rcp85", 
                              ## "HadGEM2-ES-r2i1p1_rcp85",  making touble with dates in NetCDF file
                              ## "HadGEM2-ES-r3i1p1_rcp85", making touble with dates in NetCDF file
                              "HadGEM2-ES-r4i1p1_rcp85", 
                              "INM-CM4-r1i1p1_rcp85", 
                              "IPSL-CM5A-LR-r1i1p1_rcp85", 
                              "IPSL-CM5A-LR-r2i1p1_rcp85", 
                              "IPSL-CM5A-LR-r3i1p1_rcp85", 
                              "IPSL-CM5A-LR-r4i1p1_rcp85", 
                              "IPSL-CM5A-MR-r1i1p1_rcp85", 
                              "IPSL-CM5B-LR-r1i1p1_rcp85", 
                              "MIROC5-r1i1p1_rcp85", 
                              "MIROC5-r2i1p1_rcp85", 
                              "MIROC5-r3i1p1_rcp85", 
                              ## "MIROC5-r4i1p1_rcp85",  too short timeseries...
                              ## "MIROC5-r5i1p1_rcp85", 
                              "MIROC-ESM-r1i1p1_rcp85", 
                              "MIROC-ESM-CHEM-r1i1p1_rcp85", 
                              "MPI-ESM-LR-r1i1p1_rcp85", 
                              "MPI-ESM-LR-r2i1p1_rcp85", 
                              "MPI-ESM-LR-r3i1p1_rcp85", 
                              "MPI-ESM-MR-r1i1p1_rcp85", 
                              "MRI-CGCM3-r1i1p1_rcp85", 
                              "NorESM1-M-r1i1p1_rcp85", 
                              "NorESM1-ME-r1i1p1_rcp85", 
                              "MRI-ESM1-r1i1p1_rcp85")
  ##-------------------------------------------------------------------
  ## rcp85
  cmip5.rcp85.shortperiods.labels <- c("cmip5-rcp85-shortperiods", "CMIP5-RCP85-SHORTPERIODS")
  cmip5.rcp85.shortperiods.modelnames <- c("ACCESS1-0-r1i1p1_rcp85", 
                              "ACCESS1-3-r1i1p1_rcp85", 
                              "BCC-CSM1-1-r1i1p1_rcp85", 
                              "BCC-CSM1-1-m-r1i1p1_rcp85", 
                              "BNU-ESM-r1i1p1_rcp85", 
                              "CanESM2-r1i1p1_rcp85", 
                              "CanESM2-r2i1p1_rcp85", 
                              "CanESM2-r3i1p1_rcp85", 
                              "CanESM2-r4i1p1_rcp85", 
                              "CanESM2-r5i1p1_rcp85", 
                              "CCSM4-r1i1p1_rcp85", 
                              "CCSM4-r2i1p1_rcp85", 
                              "CCSM4-r3i1p1_rcp85", 
                              "CCSM4-r4i1p1_rcp85", 
                              "CCSM4-r5i1p1_rcp85", 
                              "CCSM4-r6i1p1_rcp85", 
                              "CESM1-BGC-r1i1p1_rcp85", 
                              "CESM1-CAM5-r1i1p1_rcp85", 
                              "CESM1-CAM5-r2i1p1_rcp85", 
                              "CESM1-CAM5-r3i1p1_rcp85", 
                              "CESM1-WACCM-r2i1p1_rcp85", 
                              "CESM1-WACCM-r3i1p1_rcp85", 
                              "CESM1-WACCM-r4i1p1_rcp85", 
                              "CMCC-CESM-r1i1p1_rcp85", 
                              "CMCC-CM-r1i1p1_rcp85", 
                              "CMCC-CMS-r1i1p1_rcp85", 
                              "CNRM-CM5-r10i1p1_rcp85", 
                              "CNRM-CM5-r1i1p1_rcp85", 
                              "CNRM-CM5-r2i1p1_rcp85", 
                              "CNRM-CM5-r4i1p1_rcp85", 
                              "CNRM-CM5-r6i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r10i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r1i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r2i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r3i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r4i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r5i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r6i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r7i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r8i1p1_rcp85", 
                              "CSIRO-Mk3-6-0-r9i1p1_rcp85", 
                              "EC-EARTH-r12i1p1_rcp85", 
                              "EC-EARTH-r13i1p1_rcp85", 
                              "EC-EARTH-r14i1p1_rcp85", 
                              ## "EC-EARTH-r1i1p1_rcp85", ## corrupt scen tas files
                              "EC-EARTH-r2i1p1_rcp85", 
                              "EC-EARTH-r6i1p1_rcp85", 
                              "EC-EARTH-r8i1p1_rcp85", 
                              "EC-EARTH-r9i1p1_rcp85", 
                              "FGOALS-g2-r1i1p1_rcp85", 
                              "FIO-ESM-r1i1p1_rcp85", 
                              "FIO-ESM-r2i1p1_rcp85", 
                              "FIO-ESM-r3i1p1_rcp85", 
                              "GFDL-CM3-r1i1p1_rcp85", 
                              "GFDL-ESM2G-r1i1p1_rcp85", 
                              "GFDL-ESM2M-r1i1p1_rcp85", 
                              "GISS-E2-H-r1i1p1_rcp85", 
                              "GISS-E2-H-r1i1p2_rcp85", 
                              "GISS-E2-H-r1i1p3_rcp85", 
                              "GISS-E2-H-r2i1p1_rcp85", 
                              "GISS-E2-H-r2i1p3_rcp85", 
                              "GISS-E2-H-CC-r1i1p1_rcp85", 
                              "GISS-E2-R-r1i1p1_rcp85", 
                              "GISS-E2-R-r1i1p2_rcp85", 
                              "GISS-E2-R-r1i1p3_rcp85", 
                              "GISS-E2-R-r2i1p1_rcp85", 
                              "GISS-E2-R-r2i1p3_rcp85", 
                              "GISS-E2-R-CC-r1i1p1_rcp85", 
                              "HadGEM2-AO-r1i1p1_rcp85", 
                              "HadGEM2-CC-r1i1p1_rcp85", 
                              "HadGEM2-CC-r2i1p1_rcp85", 
                              "HadGEM2-CC-r3i1p1_rcp85", 
                              "HadGEM2-ES-r1i1p1_rcp85", 
                              "HadGEM2-ES-r2i1p1_rcp85", 
                              "HadGEM2-ES-r3i1p1_rcp85", 
                              "HadGEM2-ES-r4i1p1_rcp85", 
                              "INM-CM4-r1i1p1_rcp85", 
                              "IPSL-CM5A-LR-r1i1p1_rcp85", 
                              "IPSL-CM5A-LR-r2i1p1_rcp85", 
                              "IPSL-CM5A-LR-r3i1p1_rcp85", 
                              "IPSL-CM5A-LR-r4i1p1_rcp85", 
                              "IPSL-CM5A-MR-r1i1p1_rcp85", 
                              "IPSL-CM5B-LR-r1i1p1_rcp85", 
                              "MIROC5-r1i1p1_rcp85", 
                              "MIROC5-r2i1p1_rcp85", 
                              "MIROC5-r3i1p1_rcp85", 
                              "MIROC5-r4i1p1_rcp85", 
                              "MIROC5-r5i1p1_rcp85", 
                              "MIROC-ESM-r1i1p1_rcp85", 
                              "MIROC-ESM-CHEM-r1i1p1_rcp85", 
                              "MPI-ESM-LR-r1i1p1_rcp85", 
                              "MPI-ESM-LR-r2i1p1_rcp85", 
                              "MPI-ESM-LR-r3i1p1_rcp85", 
                              "MPI-ESM-MR-r1i1p1_rcp85", 
                              "MRI-CGCM3-r1i1p1_rcp85", 
                              "NorESM1-M-r1i1p1_rcp85", 
                              "NorESM1-ME-r1i1p1_rcp85", 
                              "MRI-ESM1-r1i1p1_rcp85")
  ##-------------------------------------------------------------------


  ##-------------------------------------------------------------------
  ## ENSEMBLES RCMs
  ensembles.labels <- c("ENSEMBLES", "ensembles")
  ensembles.modelnames <- c("METO-HC_HadRM3Q0",
                            "METO-HC_HadRM3Q16",
                            "METO-HC_HadRM3Q3",
                            "ETHZ-CLM",
                            "METNOHIRHAM_HadCM3Q0",
                            "METNOHIRHAM_BCM",
                            "ICTP-REGCM3",
                            "MPI-M-REMO",
                            "C4IRCA3",
                            "CNRM-RM5.1",
                            "CNRM-RM4.5",
                            "DMI-HIRHAM5_ARPEGE",
                            "DMI-HIRHAM5_ECHAM5",
                            "DMI-HIRHAM5_BCM",
                            "GKSS-CCLM4.8",
                            "KNMI-RACMO2",
                            "OURANOSMRCC4.2.1",
                            "SMHIRCA_BCM",
                            "SMHIRCA_ECHAM5-r3",
                            "SMHIRCA_HadCM3Q3",
                            "UCLM-PROMES",
                            "VMGO-RRCM")
  ##-------------------------------------------------------------------


  ##-------------------------------------------------------------------
  ## ENSEMBLES RCMs QM corrected Schoener
  ensembles.qmschoener.labels <- c("ENSEMBLES_QM_SCHOENER",
                                   "ensembles_qm_schoener")
  ensembles.qmschoener.modelnames <- c('METO-HC_HadRM3Q0_QM_SCHOENER',
                                       'METO-HC_HadRM3Q16_QM_SCHOENER',
                                       'METO-HC_HadRM3Q3_QM_SCHOENER',
                                       'ETHZ-CLM_QM_SCHOENER',
                                       'METNOHIRHAM_HadCM3Q0_QM_SCHOENER',
                                       'METNOHIRHAM_BCM_QM_SCHOENER',
                                       'ICTP-REGCM3_QM_SCHOENER',
                                       'MPI-M-REMO_QM_SCHOENER',
                                       'C4IRCA3_QM_SCHOENER',
                                       'CNRM-RM5.1_QM_SCHOENER',
                                       'CNRM-RM4.5_QM_SCHOENER',
                                       'DMI-HIRHAM5_ARPEGE_QM_SCHOENER',
                                       'DMI-HIRHAM5_ECHAM5_QM_SCHOENER',
                                       'DMI-HIRHAM5_BCM_QM_SCHOENER',
                                       'GKSS-CCLM4.8_QM_SCHOENER',
                                       'KNMI-RACMO2_QM_SCHOENER',
                                       'OURANOSMRCC4.2.1_QM_SCHOENER',
                                       'SMHIRCA_BCM_QM_SCHOENER',
                                       'SMHIRCA_ECHAM5-r3_QM_SCHOENER',
                                       'SMHIRCA_HadCM3Q3_QM_SCHOENER',
                                       'UCLM-PROMES_QM_SCHOENER',
                                       'VMGO-RRCM_QM_SCHOENER',
                                       'AIT-CCLM_QM_SCHOENER',
                                       'WEGC-CCLM_QM_SCHOENER')
  ##-------------------------------------------------------------------

  ## extract all known modelnames from dictionary
  dictionary.names <- names(model.input)
  ## check whether modelnames are unique, give error else
  if (length(unique(dictionary.names)) != length(model.input)){
    msg1 <- "MODELNAMES IN DICTIONARY MUST BE UNIQUE! "
    msg2 <- "THINK ABOUT A CREATIVE MODELNAME ;)"
    stop(msg1, msg2)
  }

  ## models matching dictnames
  model.input[modelnames[modelnames %in% dictionary.names]]

  unknown.modelnames <- NULL
  unknown.modelnames <- modelnames[! modelnames %in% dictionary.names]
  known.modelnames   <- modelnames[  modelnames %in% dictionary.names]

  ## not a single modelname is found in dictionary, set list to length 0, else
  ## read the funky shit!
  match.modelnames <- match(known.modelnames,modelnames)
  model.list <- model.input[modelnames[match.modelnames]]

  ## if there is at least one model not found in InitModelsDictionary then
  ## look it up in labels defined previously in this function
  if (length(unknown.modelnames) > 0) {
    if (any(unknown.modelnames %in% cmip3.sresa2.labels)) {
      model.list <- c(model.list, model.input[cmip3.sresa2.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip3.sresa2.labels]
    }
    if (any(unknown.modelnames %in% cmip3.sresb1.labels)) {
      model.list <- c(model.list, model.input[cmip3.sresb1.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip3.sresb1.labels]
    }
    if (any(unknown.modelnames %in% cmip3.sresa1b.labels)) {
      model.list <- c(model.list, model.input[cmip3.sresa1b.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip3.sresa1b.labels]
    }
    if (any(unknown.modelnames %in% cmip5.rcp26.labels)) {
      model.list <- c(model.list, model.input[cmip5.rcp26.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip5.rcp26.labels]
    }
   if (any(unknown.modelnames %in% cmip5.rcp26.shortperiods.labels)) {
      model.list <- c(model.list, model.input[cmip5.rcp26.shortperiods.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip5.rcp26.shortperiods.labels]
    }
    if (any(unknown.modelnames %in% cmip5.rcp45.labels)) {
      model.list <- c(model.list, model.input[cmip5.rcp45.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip5.rcp45.labels]
    }
   if (any(unknown.modelnames %in% cmip5.rcp45.shortperiods.labels)) {
      model.list <- c(model.list, model.input[cmip5.rcp45.shortperiods.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip5.rcp45.shortperiods.labels]
    }
     if (any(unknown.modelnames %in% cmip5.rcp60.labels)) {
      model.list <- c(model.list, model.input[cmip5.rcp60.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip5.rcp60.labels]
    }
   if (any(unknown.modelnames %in% cmip5.rcp60.shortperiods.labels)) {
      model.list <- c(model.list, model.input[cmip5.rcp60.shortperiods.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip5.rcp60.shortperiods.labels]
    }
    if (any(unknown.modelnames %in% cmip5.rcp85.labels)) {
      model.list <- c(model.list, model.input[cmip5.rcp85.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip5.rcp85.labels]
    }
   if (any(unknown.modelnames %in% cmip5.rcp85.shortperiods.labels)) {
      model.list <- c(model.list, model.input[cmip5.rcp85.shortperiods.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% cmip5.rcp85.shortperiods.labels]
    }
    if (any(unknown.modelnames %in% ensembles.labels))  {
      model.list <- c(model.list, model.input[ensembles.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% ensembles.labels]
    }
    if (any(unknown.modelnames %in% ensembles.qmschoener.labels))  {
      model.list <- c(model.list, model.input[ensembles.qmschoener.modelnames])
      unknown.modelnames <-
        unknown.modelnames[!unknown.modelnames %in% ensembles.qmschoener.labels]
    }
  }

  ##  if there is still at least one model not found in the dictionary
  if (length(unknown.modelnames) > 0) {
    msg1 <- paste("MODEL(S) \"", paste(unknown.modelnames, collapse="\" \""),
                 "\" UNKNOWN\n", sep="")
    msg2 <- "ADD THE MODEL TO DICTIONARY FILE, "
    msg3 <- "OR ADD AS A NEW ENSEMBLE INTO ReadFromModelDictionary."
    stop(msg1, msg2, msg3)
  }
  ## there's at least one NA -name in model.list, thats wrong for sure!
  if (any(is.na(names(model.list))))
    stop("THERE IS AT LEAST ONE MODEL CALLED \"NA\" IN THE MODEL.INPUT.")
  ## empty list...
  if (length(model.list) == 0)
    stop("NO MODELS HAVE BEEN SELECTED")

  return(model.list)
}


###################### handling WUX dataset ############################

MakeModelDataFrame <- function(model.name, institute.name,
                               rcm.name, global.model.name,
                               gcm.run,
                               is.reference.period,
                               time.step,
                               season = "none",
                               is.time.series,
                               subregion = "no", is.corrected = NULL,
                               emission.scenario = "A1B", resolution,
                               data, parameter.name) {
  ## Merges aggregated model data ("data") with all the information provided for
  ## the particular climate model and returns a data.frame. This is a crucial
  ## function where loads can go wrong, so pay EXTRA attention when changing
  ## this function!
  ##
  ## Args:
  ##   model.name: Name the model shall be called
  ##   institute.name: Name of the institute running the RCM
  ##   rcm.name: RCM name
  ##   global.model.name: GCM name
  ##   is.reference.period: Is given period the reference period? (yes/no)
  ##   season: string vector with second statistic timestamps
  ##                            e.g. seasons c("DJF", "MAM", "JJA", "SON")
  ##   subregion: string vector with subregion names
  ##   is.corrected: Do we have error corrected data?
  ##   emission.scenario: string with emmission scenario,
  ##   data: list of spatially aggregated data. It has to have EXACTLY the form:
  ##         list of.... subregions containing named seasonal lists with
  ##         containing aggregated data
  ##
  ## Returns:
  ##   Data.Frame with model information and spatio-temporal aggregated data.
  ##
  ## History:
  ##   2010-10-20 | Original code (thm)
  ##   2010-11-03 | renaming variables and set default is.corrected
  ##                "resolution" added
  ##   2010-11-18 | sorting arguments added (thm)
  ##   2010-12-15 | gcm.run added
  ##   2011-04-29 | restructure from MakeModelInfo to MakeModelDataFrame
  ##   2011-10-20 | NA case added for rcm name
  ##   2014-11-20 | NA case added for gcm.run and resolution

  if (is.null(data))
    stop("NO AGGREGATED DATA SPECIFIED")
  if (is.null(parameter.name))
    stop("PARAMETER NAME NOT SPECIFIED")

  ## setting default
  ##"not error corrected data"
  if (is.null(is.corrected))
    is.corrected <- "no"
  ## no resolution given
  if (is.null(resolution) || !is.null(resolution) && resolution == "")
    resolution <- NA
  ## no gcm.run given
  if (is.null(gcm.run) || !is.null(gcm.run) && gcm.run == "")
    gcm.run <- NA
  ## if no RCMs has been passed
  if (is.null(rcm.name) || !is.null(rcm.name) && rcm.name == "")
    rcm.name <- NA

  ## sort arguments to get sorted factor levels
  emission.scenario <- emission.scenario[order(emission.scenario)]
  subregion <- subregion[order(subregion)]
  season.order <- season[order(season)]

  ## generating full data.frame by expanding all possibilities
  WUX.df <- list(acronym = model.name,
                 institute = institute.name,
                 gcm =   global.model.name,
                 gcm.run = gcm.run,
                 rcm =   rcm.name,
                 em.scn = emission.scenario,
                 subreg = subregion,
                 period = time.step,
                 ref.per = is.reference.period,
                 season = season.order,
                 resolution = resolution,
                 corrected = is.corrected)

  WUX.df <- expand.grid(WUX.df)

  ## aggregated data to dataframe,
  data.df <- reshape::melt(data)
  ## these 3 dims identify the current data: the data itself, the seasons
  ## and the subregions
  if (length(data.df) != 3) {
    msg1 <- "OBVIOUSLY YOU HAVE CHANGED THE AGGREGATION DIMENSION, "
    msg2 <- "WATCH OUT WITH MERGING IN MakeModelDataFrame"
    stop(msg1, msg2)
  }
  names(data.df) <- c(parameter.name, "season", "subreg")

  ## merge dataframe-construct and aggregated parameter data
  WUX.df <- merge(WUX.df, data.df)

  ## if dealing with time series split "season" into year and period
  if (is.time.series) {
    ## extract year
    WUX.df[["year"]] <-
      factor(sapply(strsplit(as.character(WUX.df[["season"]]), " "), "[", 1))

    ## extract period
    WUX.df[["season"]] <-
      factor(sapply(strsplit(as.character(WUX.df[["season"]]), " "), "[", 2))

    ## rearrange data frame that "year" is beside "season"
    season.pos <- which(names(WUX.df) == "season")
    year.pos <-  which(names(WUX.df) == "year")
    WUX.df <- data.frame(WUX.df[1:(season.pos - 1)], year = WUX.df[year.pos],
                         WUX.df[season.pos:(length(names(WUX.df)) - 1)])
  }

  return(WUX.df)
}


MergeParameterDataframes <- function(wux.df.list, parameter.names) {
  ## Merges list of data.frames of different parameters of the same model
  ## to a single data.frame.
  ##
  ## Args:
  ##   wux.df.list: Named list containing data.frames. Each list entry
  ##              corresponds to one climate parameter with its name
  ##              being the CF-convention parameter name.
  ##   parameter.names: Character vector containing CF-convention parameter
  ##                    names to be extracted from wux.df.list.
  ##                    parameter.names have to match names(wux.df.list).
  ##
  ## Returns:
  ##   Single data.frame with model.info columns and one column for each
  ##   parameter.
  ##
  ## History:
  ##   2010-10-27 | Original code (thm)
  ##   2011-11-24 | added missing list entry error message (thm)
  ##   2011-12-21 | added case for models with empty wux.df.list (thm)

  ## if any list entry is NULL
  if (any(sapply(wux.df.list, is.null))){
    msg1 <- "AT LEAST ONE DATA.FRAME FOR ONE PARAMETER IS MISSING. "
    msg2 <- "THIS REALLY SHOULD NOT HAPPPEN."
    stop(msg1, msg2)
  }

  ## number of parameters
  n.of.pars <- length(parameter.names)

  ## sort length of parameters, longest parameter first
  ## --> important for merging
  par.length <- c(1:n.of.pars)
  for (ii in c(1:n.of.pars)) {
    if (length(wux.df.list[[ii]])== 0)
      par.length[ii] <- 0
    else
      par.length[ii] <- length(wux.df.list[[ii]][[1]])
  }
  sort.indices <- sort(par.length, decreasing = TRUE, index.return = TRUE)
  parameter.names <- parameter.names[sort.indices$ix]
  par.length <- par.length[sort.indices$ix]
  ## initializing dataframe with first element of list
  data <- wux.df.list[[parameter.names[1]]]

  ## always merge pairwise
  ii <- 2
  while (ii <= n.of.pars) {
    ## case of empty dataframe (i.e. this parameter doesnt exist for this model)
    if (length(wux.df.list[[parameter.names[ii]]]) == 0){
      data <- cbind(data, NA)
      names(data)[length(data)] <- parameter.names[ii]
    } else {
      data <- merge(data, wux.df.list[[parameter.names[ii]]], all.x = TRUE)
   }
     ii <- ii + 1
   }

  return(data)
}


WuxDataDiff <- function(wux.df, parameters) {
  ## Calculating differences between scenario period and reference period
  ##  data values within 'wux.df' data.frame.
  ##
  ## Args:
  ##   wux.df: WUX dataframe with ref.per values yes and no.
  ##   user.input: input file
  ##
  ## Returns:
  ##   WUX data.frame with cliamte change signals.
  ##
  ## History:
  ##   2010-10-27 | Original code (thm)
  ##   2010-11-24 | Only one reference period taken (thm)
  ##   2011-11-24 | aestethic corrections (thm)
  ##   2014-11-20 | all subset commands changed to '[' operation (thm)

  ## sort data frame by season, then by time, and at
  ## last by modelname (acronym)
  sortby <- c("acronym", "period", "season")
  wux.df <-
    wux.df[order(wux.df$acronym, wux.df$period, wux.df$season),]

  ## initialize diffs dataset, take all columns except the parameter
  ## data colummns
  col.noparam <- which(! names(wux.df) %in% parameters)
  wux.diff.df <- wux.df[wux.df[["ref.per"]] == "no", col.noparam]

  for (param in parameters) {

    ## initializing vector obtaining all parameter-diff values and
    ## percentage changes
    diff.column <- NULL
    perc.column <- NULL

    ## get times of reference and scenarios
    timesteps.ref <- unique(wux.df[wux.df[["ref.per"]] == "yes", ]$period)
    timesteps.scn <- unique(wux.df[wux.df[["ref.per"]] == "no" , ]$period)

    ## always take the refernce period as one single ref period,
    ## regardless of length(timesteps.ref)
    number.of.timesteps.ref <- 1
    number.of.timesteps.scn <- length(timesteps.scn)

    ## so far we can only handle one reference period....
    stopifnot(number.of.timesteps.ref == 1)

    ## we have as many data for in the reference period as in the scenario
    ## period
    x <- number.of.timesteps.ref / number.of.timesteps.scn

    ## we have x- times more scenario data than reference period data, where
    ## x has to be integer.
    ## x == round(x) tests whether we deal with integers or not (the R solution)
    if(x == round(x)){
      ## calculate difference for every scenario timestep and the
      ## corresponding reference period
      for (scn.period in timesteps.scn) {

        scn <- wux.df[wux.df[["ref.per"]] == "no" &
                      wux.df[["period"]] == scn.period
                      , c("acronym", "season", param)]
        ref <- wux.df[wux.df[["ref.per"]] == "yes", c("acronym", "season", param)]

        ## calculating differences / percentual change
        diff <- scn[[param]] - ref[[param]]
        perc <- (scn[[param]] / ref[[param]] - 1) * 100

        ## appending to previous timesteps
        diff.column <- c(diff.column, diff)
        perc.column <- c(perc.column, perc)
      }
    } else {
      ## else x was no integer (having x-times more scenario data than
      ## reference period data)
      stop("WuxDataDiff, SOMETHING WENT SERIOUSLY WRONG")
    }

    ## append diff column to diff dataset
    wux.diff.df <- cbind(wux.diff.df, diff.column)
    diff.cols <- which(names(wux.diff.df) == "diff.column")
    names(wux.diff.df)[diff.cols] <- paste("delta.", param, sep="")
    ## get percentual change for precipitation parameter
    if (param == "precipitation_amount") {
      wux.diff.df <- cbind(wux.diff.df, perc.column)
      perc.cols <- which(names(wux.diff.df) == "perc.column")
      names(wux.diff.df)[perc.cols] <- paste("perc.delta.", param, sep="")
    }
  }
  return(wux.diff.df)
}
