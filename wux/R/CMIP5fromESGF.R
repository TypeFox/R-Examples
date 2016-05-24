CMIP5fromESGF <- function(save.to = NULL,
                          variables = NULL,
                          experiments = NULL,
                          models = NULL) {
  ## Downloads available monthly CMIP5 simulations from ESGF data portal
  ## (http://esgf-data.dkrz.de). This
  ## function therefore creates subdircetories for each climate
  ## simulation in the specified folder, autopmatically recieves the
  ## bash scripts needed for the partiular
  ## simulation-variable-experiment combination and executes the bash
  ## scripts one by one. An external python script is called for this task.
  ## If either data or bash-scripts exist, the download will be skipped.
  ## Use this function with care.
  ##
  ## Args:
  ##   save.to: character. Directory where to download data and bash scripts.
  ##            ATTENTION: subdirectories for each model-experiment combination
  ##            will be created!
  ##   variables: character (vecotor). Short variable names for
  ##              meteorological parameters of interest (eg. "tas" for 2m air
  ##              temperature or "pr" for precipitation amount). See eg.
  ##              the IPCC Standard Output from GCMs
  ##              (http://www-pcmdi.llnl.gov/ipcc/standard_output.html).
  ##   experiments: character (vector). Experiment of the climate simulation
  ##                (e.g. c("historical", "rcp45"), see Taylor (2012) for
  ##                detailed description.
  ##   models:  character (vector). Climate simulations to be downloaded. If no
  ##            models are provided (default), all available simulations will be
  ##            retrieved.  ATTENTION: This is a considerable amount of data, so
  ##            watch out for your diskspace!
  ##
  ## Returns:
  ##   Nothing.
  ##
  ## History:
  ##   2014-11-06 | imported original code (thm)
  ##
  
  if (is.null(experiments))
    stop("please provide experiments (e.g. c(\"historical\", \"rcp45,\")")
  if (is.null(save.to))
    stop("please provide a directoy for data storage (save.to parameter).")
  if (is.null(variables))
    stop("please provide meteorological variables (e.g. \"tas\")")
  
  python.filename <- system.file('exec', 'CMIP5_downloader.py', package='wux')
  
  if (is.null(models))
    models <- ''
  for (model in models) {
    for (variable in variables){
      for (pathway in experiments) {
        if (model == '') {
          command <- paste("python", python.filename,
                           "-p", pathway,
                           "-v", variable,
                           "-d", save.to)
        } else {
          command <- paste("python", python.filename,
                           "-p", pathway,
                           "-v", variable,
                           "-d", save.to,
                           "-m", model)
        }
        ## execute python script
        system(command)
      }
    }
  }
  
}
