
# ----------------------------------------------------------------
# $Author: thm $
# $Date: 2016-01-13 15:46:23 +0100 (Wed, 13 Jan 2016) $
# $Rev: 344 $
# ----------------------------------------------------------------

############## Filename handling ###########

GetFileNames <- function(current.model.input,
                         period.begin, period.end,
                         parameter.name,
                         period, what.timesteps = NULL, ...) {
  ## Finds out what files we need for given time period by reading in the time
  ## vector of all NetCDF files supplied in current.model.input.
  ##
  ## Args:
  ##   current.model.input: List. Model input from InitModelsDictionary for one
  ##                        single climate model (i.e. acronym). Filenames
  ##                        will be extracted from this list structure.
  ##   period.begin: POSIXct time object determining the first timestep which
  ##                 will be read in in future.
  ##   period.end: POSIXct time object determining the last timestep which
  ##               will be read in in future.
  ##   parameter.name: Named character. Short parameter name for NetCDF
  ##                   meteorological parameter. Its name is the CF convention
  ##                   parameter name.
  ##   period: Character. "reference.period" or "scenario.perio" to extract
  ##           information from current.model.input.
  ##   ...: Arguments being passed to ReadNetCdfTimeData (maybe different
  ##        calender or time units).
  ##
  ## Returns:
  ##   Character vector of full filenames for model from current.model.input
  ##   ranging from period.begin to period.end.
  ##
  ## History:
  ##   2011-01-25 | original code (thm)
  ##   2013-10-08 | potential pitfall fixed if file.name and file.path
  ##                have a different list structure (thm)
  ##   2015-06-17 | renamed modelinput directory entries from "reference.period" to "historical" and
  ##                "scenario.period" to "scenario"  (thm)
  ##
  ## TODO(thm) give warning or errors if not enough files
  ##       and check for duplicate timeseries! and maybe sort them...

    ## reference periods are marked as "historical" in modelinput file
    if (period == "reference.period")
        period <- "historical"
    if (period == "scenario.period")
        period <- "scenario"
    
  ## extract parameter we are looking for
  parameter.long <- names(parameter.name)

  ## get appropriate filepath (see manual for distinction between file.path.alt
  ## and file.path.default
  file.path.alt <- current.model.input$file.path.alt
  file.path.default <- current.model.input$file.path.default
  if (is.null(file.path.alt) & is.null(file.path.default)){
    msg1 <- "NEITHER filepath.default NOR filepath.alt GIVEN! "
    msg2 <- "ENTER in InitModelDictionary."
    stop(msg1, msg2)
  }

### first extract the file.path for current parameter and current
### period ("scenario" or "historical") for given model

  ## set filepath to either file.path.alt or file.path.default (whatever has
  ## been definied in user.input)
  if (is.null(file.path.alt)) {
    file.path <- file.path.default
  } else {
    file.path <- file.path.alt
    ## case we have more filepaths (according to parameters or periods)
    is.filepath.vector.or.list <- length(file.path) > 1 | is.list(file.path)
    if (is.filepath.vector.or.list) {
      ## if file.path is list, extract in the exact order parameter-period
      ## else file.path is a vector
      if (is.list(file.path)) {
        ## file.path <- file.path[[parameter.long]][[period]]
        file.path <- file.path[[parameter.long]]
      } else {
        ## are the filepaths different for parameters or periods?
        is.split.by.period <- any(grepl(period, names(file.path)))
        is.split.by.parameter <- any(grepl(parameter.long, names(file.path)))
        ## extract filepath for current parameter/period
        if (!is.split.by.period & !is.split.by.parameter)
          ## case if neither current parameter nor current period where found
          ##  -> use file.path.default instead of file.path.alt
          file.path <- file.path.default
        else if (is.split.by.period)
          ## vector element of current period
          file.path <- file.path[period]
        else if (is.split.by.parameter)
          ## we found the current parameter in the path vector
          file.path <- file.path[parameter.long]
        else
          ## if a file.path.alt is a vector containing both a name attribute
          ## as the current parameter AND the current period, which is nonsense
          stop("YOU SHOULDN'T HAVE GOTTEN HERE... SOMETHING IS WRONG WITH YOUR FIL.PATH SETTING IN InitModelDictionary")
      }
      ## if no data are available (the according file.path
      ## in InitModelDictionary is set NA)
      if ( all(is.na(file.path)) ) {
        return(NA)
      }
    }
  }

  ## file.path on this machine?
  for (path.count in c(1:length(file.path))) {
    if (!file.exists(file.path[[path.count]])) {
      msg1 <- "ERROR WHILE GETTING FILENAMES FROM InitModelDictionary:FILEPATH \""
      msg2 <- "\" DOES NOT EXIST"
      stop(msg1, file.path, msg2)
    }
  }
### now the filenames will be extracted for current parameter and current
### period ("scenario" or "historical") for given model

  ## we try to figure out whether file.names are stored hierarchically in a list
  ## as eg. [[parameter]][[period or [[period]][[parameter]] or [[period]] on
  ## its own. file.name.tmp will be overwritten every time the hierarchy
  ## structure is obviously wrong (file.name.tmp being NULL)
  ## file.name is a safety copy used to overwrite file.name.tmp.
  file.name <- current.model.input$file.name
  file.name.tmp <- file.name
  if (is.list(file.name)) {
    ## if elements of list are lists as well... dig deeper, else take this list
    file.name.tmp <- file.name[[parameter.long]]

    ## hierarchy file.name[[period]]([[parameter.long]])?
    if (is.null(file.name.tmp)) {
      file.name.tmp <- file.name[[period]]
      if (is.null(file.name.tmp))
        stop("INVALID file.name ATTRIBUTE SET")
      ## file.name[[period]][[parameter.long]] exists, so extract that
      if (!is.null(file.name.tmp[[parameter.long]]))
        file.name.tmp <- file.name.tmp[[parameter.long]]
    }

    ## hierarchy file.name[[parameter.long]]([[period]])?
    if (is.list(file.name.tmp)) {
      ## file.name[[parameter.long]][[period]] exists, so extract that
      if (!is.null(file.name.tmp[[period]])) {
        ## file.name.tmp <- file.name.tmp[[period]]
        file.name.tmp <- file.name.tmp
      }
      ## else leave file.name.tmp as it is
    }
  }

  ## we where not successful in finding filenames
  if (is.null(file.name.tmp))
    stop("NO FILENAME IN InitModelsDictionary")

  ## if no data are available (in this case it means that the according
  ## file.name in InitModelDictionary is set NA) this will omit reading in the
  ## data, thus generating NAs in the climate change signals later. This is
  ## needed if a particular parameter does not exist for this particular model
  if ( all(is.na(file.name.tmp)) ) {
    return(NA)
  }

  ## check if files "file.name.tmp" are in the given directory.
  ## here regular expressions are being checked as well
  if (is.list(file.path)) {
     for (period.name in names(file.path)) {
       file.name.tmp[[period.name]] <-
         sapply(file.name.tmp[[period.name]],
                function(x) list.files(file.path[[period.name]], pattern = x),
                USE.NAMES = FALSE)
       ## if no files in "file.path", then a list containing one 0 length variable
       ## has been returned
       if (any(sapply(file.name.tmp[[period.name]], length) == 0))
         stop("NO FILE FOUND")
     }
   } else {
     ## WARNING, file.path is not a list, whereas file.name.tmp IS a list.
     ## we therefore unlist file.name.tmp not to loop in list entries
     ## as e.g. the periods, but instead to loop over the proper filename
     ## entries, for all (refrence and scenario) periods simultanously
     if (is.list(file.name.tmp)) 
       file.name.tmp <- unlist(file.name.tmp)
     ## check for filnames
     file.name.tmp <-
       sapply(file.name.tmp,
              function(x) list.files(file.path, pattern = x),
              USE.NAMES = FALSE)
     ## if no files in "file.path", then a list containing one 0 length variable
     ## has been returned
     if (any(sapply(file.name.tmp, length) == 0))
       stop("NO FILE FOUND, MAYBE YOU HAVE A TYPO IN InitModelsDictionary.")
   }

### check which files from file.name.tmp are within startdate and enddate
### we keep those and return the full filename (including path)
  file.name.tmp <- ScanThroughAllTimeVectors(file.names = file.name.tmp,
                                             file.path = file.path,
                                             parameter = parameter.name,
                                             startdate = period.begin,
                                             enddate = period.end,
                                             what.timesteps = what.timesteps,
                                             act.period = period,
                                             ...)

  if (length(file.name.tmp) == 0)
    stop("NO FILENAMES FOUND, SOMETHINGS WRONG WITH YOUR INPUT OR WUX...")

### excitedly we pass the filenames furter to wux for a glorious future of this
### very model. thou shall be processed.
  return(file.name.tmp)
}


ScanThroughAllTimeVectors <- function(file.names, file.path,
                                      parameter, startdate, enddate,
                                      what.timesteps, act.period, ...) {
  ## Drops filenames not being within startdate and enddate.
  ##
  ## Args:
  ##   file.names: Character vector of filenames being checked.
  ##   file.path:  Character filepath where file.names are stored.
  ##   parameter:  Named character. Short parameter name for NetCDF
  ##               meteorological parameter. Its name is the CF convention
  ##               parameter name.
  ##   startdate: POSIXct time object determining the first timestep which
  ##              will be read in in future.
  ##   enddate:   POSIXct time object determining the last timestep which
  ##              will be read in in future.
  ##   ...: Arguments being passed to ReadNetCdfTimeData (maybe different
  ##        calender or time units).
  ##
  ## Returns:
  ##   Vector of full filenames (including directory) which lie within startdate
  ##   and enddate.
  ##
  ## History:
  ##   2011-01-25 | original code (thm)
  ##   2014-11-21 | platform independent (thm)

  ## get filenames containing "file.name"
  if (is.list(file.path)) {
    ls.with.path <- list()
    for (period.name in names(file.path)) {
      tmp.ls.with.path <- list()
      tmp.path <- paste(file.path[[period.name]], file.names[[period.name]], sep="/")
      tmp.ls.with.path[tmp.path] <- period.name
      ls.with.path <- c(ls.with.path, tmp.ls.with.path)
    }
    rm(tmp.path)
  } else {
    ## ls.filtered <- grep(paste(file.name, collapse="|"), ls, value=TRUE)
    ls.with.path <- file.path(file.path, file.names)
  }

  ## and extract all their time vectors
  timevectors <- ReadNetCdfTimeData(ls.with.path,
                                    what.timesteps=what.timesteps, ...)
  timevectors <- GetTimeCountsOffsets(timevectors, startdate, enddate, check.period=act.period)

  ## get filenames whose time-vector is within start- and enddate
  files.within <- lapply(timevectors, function(x) !is.na(x[["count"]]))
  files.within <- unlist(files.within)
  filenames.in.timeperiod <- names(files.within[files.within==TRUE])

  return(filenames.in.timeperiod)
}

