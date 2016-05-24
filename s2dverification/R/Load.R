Load <- function(var, exp = NULL, obs = NULL, sdates, nmember = NULL, 
                 nmemberobs = NULL, nleadtime = NULL, leadtimemin = 1, 
                 leadtimemax = NULL, storefreq = 'monthly', sampleperiod = 1, 
                 lonmin = 0, lonmax = 360, latmin = -90, latmax = 90, 
                 output = 'areave', method = 'conservative', grid = NULL, 
                 maskmod = vector("list", 15), maskobs = vector("list", 15), 
                 configfile = NULL, varmin = NULL, varmax = NULL, 
                 silent = FALSE, nprocs = NULL, dimnames = NULL, 
                 remapcells = 2) {
  #library(parallel)
  #library(bigmemory)

  # Print a stamp of the call the user issued.
  parameter_names <- ls()
  if (length(parameter_names) < 3 || is.null(var) ||
      is.null(sdates) || (is.null(exp) && is.null(obs))) {
    stop("Error: At least 'var', 'exp'/'obs' and 'sdates' must be provided.")
  }
  load_parameters <- lapply(parameter_names, get, envir = environment())
  names(load_parameters) <- parameter_names
  parameters_to_show <- c('var', 'exp', 'obs', 'sdates', 'grid', 'output', 'storefreq')
  load_parameters <- c(load_parameters[parameters_to_show], load_parameters[-match(parameters_to_show, names(load_parameters))])
  cat(paste("* The load call you issued is:\n*   Load(", 
            paste(strwrap(
              paste(unlist(lapply(names(load_parameters[1:length(parameters_to_show)]), 
              function(x) paste(x, '=', 
                if (x == 'sdates' && length(load_parameters[[x]]) > 4) {
                  paste0("c('", load_parameters[[x]][1], "', '", load_parameters[[x]][2], 
                         "', ..., '", tail(load_parameters[[x]], 1), "')")
                } else {
                  paste(deparse(load_parameters[[x]]), collapse = '')
                }))), 
              collapse = ', '), width = getOption('width') - 9, indent = 0, exdent = 8), collapse = '\n*'),
            ", ...)\n* See the full call in '$load_parameters' after Load() finishes.\n", sep = ''))

  # Run Load() error-aware, so that it always returns something
  errors <- try({

  # Check and sanitize parameters
  # var
  if (is.null(var) || !(is.character(var) && nchar(var) > 0)) {
    stop("Error: parameter 'var' should be a character string of length >= 1.")
  }

  # exp
  exps_to_fetch <- c()
  exp_info_names <- c('name', 'path', 'nc_var_name', 'suffix', 
                      'var_min', 'var_max', 'dimnames')
  if (!is.null(exp) && !(is.character(exp) && all(nchar(exp) > 0)) && !is.list(exp)) {
    stop("Error: parameter 'exp' should be a vector of strings or a list with information of the experimental datasets to load. Check 'exp' in ?Load for details.")
  } else if (!is.null(exp)) {
    if (!is.list(exp)) {
      exp <- lapply(exp, function (x) list(name = x))
    }
    for (i in 1:length(exp)) {
      if (!is.list(exp[[i]])) {
        stop("Error: parameter 'exp' is incorrect. It should be a list of lists.")
      }
      if (!(all(names(exp[[i]]) %in% exp_info_names))) {
        stop("Error: parameter 'exp' is incorrect. There are unrecognized components in the information of some of the experiments. Check 'exp' in ?Load for details.")
      }
      if (!('name' %in% names(exp[[i]]))) {
        exp[[i]][['name']] <- paste0('exp', i)
        if (!('path' %in% names(exp[[i]]))) {
          stop("Error: parameter 'exp' is incorrect. A 'path' should be provided for each experimental dataset if no 'name' is provided. See 'exp' in ?Load for details.")
        }
      } else if (!('path' %in% names(exp[[i]]))) {
        exps_to_fetch <- c(exps_to_fetch, i)
      }
      if ('path' %in% names(exp[[i]])) {
        if (!('nc_var_name' %in% names(exp[[i]]))) {
          exp[[i]][['nc_var_name']] <- '$VAR_NAME$'
        }
        if (!('suffix' %in% names(exp[[i]]))) {
          exp[[i]][['suffix']] <- ''
        }
        if (!('var_min' %in% names(exp[[i]]))) {
          exp[[i]][['var_min']] <- ''
        }
        if (!('var_max' %in% names(exp[[i]]))) {
          exp[[i]][['var_max']] <- ''
        }
      }
    }
    if ((length(exps_to_fetch) > 0) && (length(exps_to_fetch) < length(exp))) {
      cat("! Warning: 'path' was provided for some experimental datasets in 'exp'. Any \n*   information in the configuration file related to these will be ignored.\n")
    }
  }

  # obs
  obs_to_fetch <- c()
  obs_info_names <- c('name', 'path', 'nc_var_name', 'suffix', 
                      'var_min', 'var_max')
  if (!is.null(obs) && !(is.character(obs) && all(nchar(obs) > 0)) && !is.list(obs)) {
    stop("Error: parameter 'obs' should be a vector of strings or a list with information of the observational datasets to load. Check 'obs' in ?Load for details.")
  } else if (!is.null(obs)) {
    if (!is.list(obs)) {
      obs <- lapply(obs, function (x) list(name = x))
    }
    for (i in 1:length(obs)) {
      if (!is.list(obs[[i]])) {
        stop("Error: parameter 'obs' is incorrect. It should be a list of lists.")
      }
      if (!(all(names(obs[[i]]) %in% obs_info_names))) {
        stop("Error: parameter 'obs' is incorrect. There are unrecognized components in the information of some of the observations. Check 'obs' in ?Load for details.")
      }
      if (!('name' %in% names(obs[[i]]))) {
        obs[[i]][['name']] <- paste0('obs', i)
        if (!('path' %in% names(obs[[i]]))) {
          stop("Error: parameter 'obs' is incorrect. A 'path' should be provided for each observational dataset if no 'name' is provided. See 'obs' in ?Load for details.")
        }
      } else if (!('path' %in% names(obs[[i]]))) {
        obs_to_fetch <- c(obs_to_fetch, i)
      }
      if ('path' %in% names(obs[[i]])) {
        if (!('nc_var_name' %in% names(obs[[i]]))) {
          obs[[i]][['nc_var_name']] <- '$VAR_NAME$'
        }
        if (!('suffix' %in% names(obs[[i]]))) {
          obs[[i]][['suffix']] <- ''
        }
        if (!('var_min' %in% names(obs[[i]]))) {
          obs[[i]][['var_min']] <- ''
        }
        if (!('var_max' %in% names(obs[[i]]))) {
          obs[[i]][['var_max']] <- ''
        }
      }
    }
    if (length(c(obs_to_fetch, exps_to_fetch) > 1) && (length(obs_to_fetch) < length(obs))) {
      cat("! Warning: 'path' was provided for some observational datasets in 'obs'. Any \n*   information in the configuration file related to these will be ignored.\n")
    }
  }

  # sdates
  if (is.null(sdates)) {
    stop("Error: parameter 'sdates' must be provided.")
  }
  if (!is.character(sdates) || !all(nchar(sdates) == 8) || any(is.na(strtoi(sdates)))) {
    stop("Error: parameter 'sdates' is incorrect. All starting dates should be a character string in the format 'YYYYMMDD'.")
  }

  # nmember
  if (!is.null(nmember) && !is.null(exp)) {
    if (!is.numeric(nmember)) {
      stop("Error: parameter 'nmember' is incorrect. It should be numeric.")
    }
    if (length(nmember) == 1) {
      cat(paste("! Warning: 'nmember' should specify the number of members of each experimental dataset. Forcing to", nmember, "for all experiments.\n"))
      nmember <- rep(nmember, length(exp))
    }
    if (length(nmember) != length(exp)) {
      stop("Error: 'nmember' must contain as many values as 'exp'.")
    } else if (any(is.na(nmember))) {
      nmember[which(is.na(nmember))] <- max(nmember, na.rm = TRUE)
    }
  }

  # nmemberobs
  if (!is.null(nmemberobs) && !is.null(obs)) {
    if (!is.numeric(nmemberobs)) {
      stop("Error: parameter 'nmemberobs' is incorrect. It should be numeric.")
    }
    if (length(nmemberobs) == 1) {
      cat(paste("! Warning: 'nmemberobs' should specify the number of members of each observational dataset. Forcing to", nmemberobs, "for all observations.\n"))
      nmemberobs <- rep(nmemberobs, length(obs))
    }
    if (length(nmemberobs) != length(obs)) {
      stop("Error: 'nmemberobs' must contain as many values as 'obs'.")
    } else if (any(is.na(nmemberobs))) {
      nmemberobs[which(is.na(nmemberobs))] <- max(nmemberobs, na.rm = TRUE)
    }
  }

  # nleadtime
  if (!is.null(nleadtime) && !is.numeric(nleadtime)) {
    stop("Error: parameter 'nleadtime' is wrong. It should be numeric.")
  }

  # leadtimemin
  if (is.null(leadtimemin) || !is.numeric(leadtimemin)) {
    stop("Error: parameter 'leadtimemin' is wrong. It should be numeric.")
  }

  # leadtimemax
  if (!is.null(leadtimemax) && !is.numeric(leadtimemax)) {
    stop("Error: parameter 'leadtimemax' is wrong. It should be numeric.")
  }

  # storefreq       
  if (!is.character(storefreq) || !(storefreq %in% c('monthly', 'daily'))) {
    stop("Error: parameter 'storefreq' is wrong, can take value 'daily' or 'monthly'.")
  }

  # sampleperiod
  if (is.null(sampleperiod) || !is.numeric(sampleperiod)) {
    stop("Error: parameter 'sampleperiod' is wrong. It should be numeric.")
  }

  # lonmin
  if (is.null(lonmin) || !is.numeric(lonmin)) {
    stop("Error: parameter 'lonmin' is wrong. It should be numeric.")
  }
  if (lonmin < -360 || lonmin > 360) {
    stop("Error: parameter 'lonmin' must be in the range [-360, 360]")
  }
  if (lonmin < 0) {
    lonmin <- lonmin + 360
  }

  # lonmax
  if (is.null(lonmax) || !is.numeric(lonmax)) {
    stop("Error: parameter 'lonmax' is wrong. It should be numeric.")
  }
  if (lonmax < -360 || lonmax > 360) {
    stop("Error: parameter 'lonmax' must be in the range [-360, 360]")
  }
  if (lonmax < 0) {
    lonmax <- lonmax + 360
  }

  # latmin
  if (is.null(latmin) || !is.numeric(latmin)) {
    stop("Error: parameter 'latmin' is wrong. It should be numeric.")
  }
  if (latmin > 90 || latmin < -90) {
    stop("Error: 'latmin' must be in the interval [-90, 90].")
  }
  
  # latmax
  if (is.null(latmax) || !is.numeric(latmax)) {
    stop("Error: parameter 'latmax' is wrong. It should be numeric.")
  }
  if (latmax > 90 || latmax < -90) {
    stop("Error: 'latmax' must be in the interval [-90, 90].")
  }

  # output
  if (is.null(output) || !(output %in% c('lonlat', 'lon', 'lat', 'areave'))) {
    stop("Error: 'output' can only take values 'lonlat', 'lon', 'lat' or 'areave'.")
  }

  # method
  if (is.null(method) || !(method %in% c('bilinear', 'bicubic', 'conservative', 'distance-weighted'))) {
    stop("Error: parameter 'method' is wrong, can take value 'bilinear', 'bicubic', 'conservative' or 'distance-weighted'.")
  }
  remap <- switch(method, 'bilinear' = 'remapbil', 'bicubic' = 'remapbic', 
                  'conservative' = 'remapcon', 'distance-weighted' = 'remapdis')

  # grid
  if (!is.null(grid)) {
    if (is.character(grid)) {
      supported_grids <- list('r[0-9]{1,}x[0-9]{1,}', 't[0-9]{1,}grid')
      grid_matches <- unlist(lapply(lapply(supported_grids, regexpr, grid), .IsFullMatch, grid))
      if (sum(grid_matches) < 1) {
        stop("The specified grid in the parameter 'grid' is incorrect. Must be one of r<NX>x<NY> or t<RES>grid.")
      }
    } else {
      stop("Error: parameter 'grid' should be a character string, if specified.")
    }
  }

  # maskmod
  if (!is.list(maskmod)) {
    stop("Error: parameter 'maskmod' must be a list.")
  }
  if (length(maskmod) < length(exp)) {
    stop("Error: 'maskmod' must contain a numeric mask or NULL for each experiment in 'exp'.")
  }
  for (i in 1:length(maskmod)) {
    if (is.list(maskmod[[i]])) {
      if ((length(maskmod[[i]]) > 2) || !all(names(maskmod[[i]]) %in% c('path', 'nc_var_name'))) {
        stop("Error: all masks in 'maskmod' must be a numeric matrix, or a list with the components 'path' and optionally 'nc_var_name', or NULL.")
      }
    } else if (!(is.numeric(maskmod[[i]]) || is.null(maskmod[[i]]))) {
      stop("Error: all masks in 'maskmod' must be a numeric matrix, or a list with the components 'path' and optionally 'nc_var_name', or NULL.")
    }
  }

  # maskobs
  if (!is.list(maskobs)) {
    stop("Error: parameter 'maskobs' must be a list.")
  }
  if (length(maskobs) < length(obs)) {
    stop("Error: 'maskobs' must contain a numeric mask or NULL for each obseriment in 'obs'.")
  }
  for (i in 1:length(maskobs)) {
    if (is.list(maskobs[[i]])) {
      if ((length(maskobs[[i]]) > 2) || !all(names(maskobs[[i]]) %in% c('path', 'nc_var_name'))) {
        stop("Error: all masks in 'maskobs' must be a numeric matrix, or a list with the components 'path' and optionally 'nc_var_name', or NULL.")
      }
    } else if (!(is.numeric(maskobs[[i]]) || is.null(maskobs[[i]]))) {
      stop("Error: all masks in 'maskobs' must be a numeric matrix, or a list with the components 'path' and optionally 'nc_var_name', or NULL.")
    }
  }

  ## Force the observational masks to be the same as the experimental when
  ## possible.
  if ((output != 'areave' || !is.null(grid)) && length(exp) > 0) {
    if (!all(unlist(lapply(maskobs, is.null)))) {
      cat("! Warning: 'maskobs' will be ignored. 'maskmod[[1]]' will be applied to observations instead.\n")
    }
    maskobs <- lapply(maskobs, function(x) x <- maskmod[[1]])
  }

  # configfile
  if (is.null(configfile)) {
    configfile <- system.file("config", "BSC.conf", package = "s2dverification")
  } else if (!is.character(configfile) || !(nchar(configfile) > 0)) {
    stop("Error: parameter 'configfile' must be a character string with the path to an s2dverification configuration file, if specified.")
  }

  # varmin
  if (!is.null(varmin) && !is.numeric(varmin)) {
    stop("Error: parameter 'varmin' must be numeric, if specified.")
  }

  # varmax
  if (!is.null(varmax) && !is.numeric(varmax)) {
    stop("Error: parameter 'varmax' must be numeric, if specified.")
  }

  # silent
  if (!is.logical(silent)) {
    stop("Error: parameter 'silent' must be TRUE or FALSE.")
  }

  # nprocs
  if (!is.null(nprocs) && (!is.numeric(nprocs) || nprocs < 1)) {
    stop("Error: parameter 'nprocs' must be a positive integer, if specified.")
  }

  # dimnames
  if (!is.null(dimnames) && (!is.list(dimnames))) {
    stop("Error: parameter 'dimnames' must be a list, if specified.")
  }
  if (!all(names(dimnames) %in% c('member', 'lat', 'lon'))) {
    stop("Error: parameter 'dimnames' is wrong. There are unrecognized component names. See 'dimnames' in ?Load for details.")
  }

  # remapcells  
  if (!is.numeric(remapcells) || remapcells < 0) {
    stop("Error: 'remapcells' must be an integer >= 0.")
  }

  # If not all data has been provided in 'exp' and 'obs', configuration file is read.
  if (length(exps_to_fetch) > 0 || length(obs_to_fetch) > 0) {
    cat("* Some 'path's not explicitly provided in 'exp' and 'obs', so will now proceed to open the configuration file.\n")
    data_info <- ConfigFileOpen(configfile, silent, TRUE)

    # Check that the var, exp and obs parameters are right and keep the entries
    # that match for each dataset.
    # Afterwards, the matching entries are applied sequentially (as specified 
    # in ?ConfigFileOpen) and the replace_values are applied to the result.
    # Finally a path pattern for each dataset is provided.
    matches <- ConfigApplyMatchingEntries(data_info, var, sapply(exp[exps_to_fetch], '[[', 'name'), 
                 sapply(obs[obs_to_fetch], '[[', 'name'), show_entries = FALSE, show_result = FALSE)
    # 'replace_values' is a named list that associates a variable name to an 
    # associated value. Initially it is filled with variables and values parsed
    # from the configuration file, but we can add or modify some values during
    # the execution to choose for example which start date we want to load.
    # When '.ConfigReplaceVariablesInString' is called, all the variable accesses 
    # ($VARIABLE_NAME$) that appear in the string given as parameter are 
    # replaced by the associated value in 'replace_values'.
    replace_values <- data_info$definitions
    if (!is.null(exp) && length(exps_to_fetch) > 0) {
      counter <- 1
      exp[exps_to_fetch] <- lapply(matches$exp_info, 
        function (x) {
          x[names(exp[[exps_to_fetch[counter]]])] <- exp[[exps_to_fetch[counter]]]
          x[['path']] <- paste0(x[['main_path']], x[['file_path']])
          counter <<- counter + 1
          x
        })
    }
    if (!is.null(obs) && length(obs_to_fetch) > 0) {
      counter <- 1
      obs[obs_to_fetch] <- lapply(matches$obs_info, 
        function (x) {
          x[names(obs[[obs_to_fetch[counter]]])] <- obs[[obs_to_fetch[counter]]]
          x[['path']] <- paste0(x[['main_path']], x[['file_path']])
          counter <<- counter + 1
          x
        })
    }
    if (!silent) {
      cat("* All pairs (var, exp) and (var, obs) have matching entries.\n")
    }
  } else {
    replace_values <- list(DEFAULT_NC_VAR_NAME = '$VAR_NAME$',
                           DEFAULT_VAR_MIN = '',
                           DEFAULT_VAR_MAX = '',
                           DEFAULT_SUFFIX = '',
                           DEFAULT_DIM_NAME_LONGITUDES = 'longitude',
                           DEFAULT_DIM_NAME_LATITUDES = 'latitude',
                           DEFAULT_DIM_NAME_MEMBERS = 'ensemble')
  }
  # We take the dimnames that haven't been explicitly specified from the 
  # configuration file.
  # If the configuration file wasn't opened, we take the default values from
  # the dictionary 'replace_values'.
  dimnames <- list(lon = ifelse(is.null(dimnames[["lon"]]), 
                                replace_values[["DEFAULT_DIM_NAME_LONGITUDES"]],
                                dimnames[['lon']]),
                   lat = ifelse(is.null(dimnames[["lat"]]), 
                                replace_values[["DEFAULT_DIM_NAME_LATITUDES"]],
                                dimnames[['lat']]),
                   member = ifelse(is.null(dimnames[["member"]]), 
                                   replace_values[["DEFAULT_DIM_NAME_MEMBERS"]],
                                   dimnames[['member']]))
  if (!is.null(exp)) {
    exp <- lapply(exp, function (x) {
      if (!('dimnames' %in% names(x))) {
        x[['dimnames']] <- dimnames
        x
      } else {
        dimnames2 <- dimnames
        dimnames2[names(x[['dimnames']])] <- x[['dimnames']]
        x[['dimnames']] <- dimnames2
        x
      }
    })
  }
  if (!is.null(obs)) {
    obs <- lapply(obs, function (x) {
      if (!('dimnames' %in% names(x))) {
        x[['dimnames']] <- dimnames
        x
      } else {
        dimnames2 <- dimnames
        dimnames2[names(x[['dimnames']])] <- x[['dimnames']]
        x[['dimnames']] <- dimnames2
        x 
      }
    })
  }
  single_dataset <- (length(obs) + length(exp) == 1)

  ## We add some predefined values in the dictionary.
  replace_values[["VAR_NAME"]] <- var
  replace_values[["STORE_FREQ"]] <- storefreq

  # Initialize some variables that will take various values along the
  # execution
  latitudes <- longitudes <- NULL
  leadtimes <- NULL
  var_exp <- var_obs <- NULL
  units <- var_long_name <- NULL
  is_2d_var <- FALSE

  # Start defining the dimensions of the output matrices
  nmod <- length(exp)
  nobs <- length(obs)
  nsdates <- length(sdates)

  # We will iterate over all the experiments, start dates and members and will open
  # the file pointed by the data in the configuration file.
  # If a file is found, we will open it and read its metadata to work out the 
  # remaining dimensions: members, leadtimes, longitudes and latitudes.
  #
  # At each iteration we will build a 'work piece' that will contain information
  # on the data we want to load from a file. For each file we will have one
  # work piece. These work pieces will be packages of information to be sent to
  # the various parallel processes. Each process will need this information to
  # access and manipulate the data according to the output type and other 
  # parameters. 
  if (!silent) {
    cat("* Fetching first experimental files to work out 'var_exp' size...\n")
  }

  dataset_type <- 'exp'
  dim_exp <- NULL
  filename <- file_found <- tmp <- nltime <- NULL
  dims2define <- TRUE
  is_file_per_member_exp <- rep(nmod, FALSE)
  exp_work_pieces <- list()
  jmod <- 1
  while (jmod <= nmod) {
    tags_to_find <- c('MEMBER_NUMBER')
    position_of_tags <- na.omit(match(tags_to_find, names(replace_values)))
    if (length(position_of_tags) > 0) {
      quasi_final_path <- .ConfigReplaceVariablesInString(exp[[jmod]][['path']], 
                            replace_values[-position_of_tags], TRUE)
    } else {
      quasi_final_path <- .ConfigReplaceVariablesInString(exp[[jmod]][['path']], 
                            replace_values, TRUE)
    }
    is_file_per_member_exp[jmod] <- grepl('$MEMBER_NUMBER$', 
                                          quasi_final_path, fixed = TRUE)
    replace_values[["EXP_NAME"]] <- exp[[jmod]][['name']]
    replace_values[["NC_VAR_NAME"]] <- exp[[jmod]][['nc_var_name']]
    namevar <- .ConfigReplaceVariablesInString(exp[[jmod]][['nc_var_name']], replace_values)
    replace_values[["SUFFIX"]] <- exp[[jmod]][['suffix']]
    if (is.null(varmin)) {
      mod_var_min <- as.numeric(.ConfigReplaceVariablesInString(exp[[jmod]][['var_min']], replace_values))
    } else {
      mod_var_min <- varmin
    }
    if (is.null(varmax)) {
      mod_var_max <- as.numeric(.ConfigReplaceVariablesInString(exp[[jmod]][['var_max']], replace_values))
    } else {
      mod_var_max <- varmax
    }
    jsdate <- 1
    while (jsdate <= nsdates) {
      replace_values[["START_DATE"]] <- sdates[jsdate]
      replace_values[["YEAR"]] <- substr(sdates[jsdate], 1, 4)
      replace_values[["MONTH"]] <- substr(sdates[jsdate], 5, 6)
      replace_values[["DAY"]] <- substr(sdates[jsdate], 7, 8)
      # If the dimensions of the output matrices are still to define, we try to read
      # the metadata of the data file that corresponds to the current iteration
      if (dims2define) {
        if (is_file_per_member_exp[jmod]) {
          replace_values[["MEMBER_NUMBER"]] <- '*'
        }
        # We must build a work piece that will be sent to the .LoadDataFile function
        # in 'explore_dims' mode. We will obtain, if success, the dimensions of the
        # data in the file.
        work_piece <- list(dataset_type = dataset_type,
                           filename = .ConfigReplaceVariablesInString(exp[[jmod]][['path']], replace_values),
                           namevar = namevar, grid = grid, remap = remap, remapcells = remapcells,
                           is_file_per_member = is_file_per_member_exp[jmod],
                           is_file_per_dataset = FALSE,
                           lon_limits = c(lonmin, lonmax),
                           lat_limits = c(latmin, latmax), dimnames = exp[[jmod]][['dimnames']],
                           single_dataset = single_dataset)
        found_data <- .LoadDataFile(work_piece, explore_dims = TRUE, silent = silent)
        found_dims <- found_data$dims
        var_long_name <- found_data$var_long_name
        units <- found_data$units
        if (!is.null(found_dims)) {
          is_2d_var <- found_data$is_2d_var
          if (!is_2d_var && (output != 'areave')) {
            cat(paste("! Warning: '", output, "' output format not allowed when loading global mean variables. Forcing to 'areave'.\n",
                sep = ''))
            output <- 'areave'
          }
          if (output != 'areave' && is.null(grid)) {
            grid <- found_data$grid
          }
          if (is.null(nmember)) {
            if (is.null(found_dims[['member']])) {
              cat("! Warning: loading data from a server but 'nmember' not specified. Loading only one member.\n")
              nmember <- rep(1, nmod)
            } else {
              nmember <- rep(found_dims[['member']], nmod)
            }
          }
          if (is.null(nleadtime)) {
            nleadtime <- found_dims[['time']]
          }
          if (is.null(leadtimemax)) {
            leadtimemax <- nleadtime
          } else if (leadtimemax > nleadtime) {
            stop("Error: 'leadtimemax' argument is greater than the number of loaded leadtimes. Put first the experiment with the greatest number of leadtimes or adjust properly the parameters 'nleadtime' and 'leadtimemax'.")
          }
          leadtimes <- seq(leadtimemin, leadtimemax, sampleperiod)
          latitudes <- found_dims[['lat']]
          longitudes <- found_dims[['lon']]
          
          if (output == 'lon' || output == 'lonlat') {
            dim_exp[['lon']] <- length(longitudes)
          }
          if (output == 'lat' || output == 'lonlat') {
            dim_exp[['lat']] <- length(latitudes)
          }
          dim_exp[['time']] <- length(leadtimes)
          dim_exp[['member']] <- max(nmember)
          dim_exp[['sdate']] <- nsdates
          dim_exp[['dataset']] <- nmod
          dims2define <- FALSE
        }
      }
      # We keep on iterating through members to build all the work pieces.
      if (is_file_per_member_exp[jmod]) {
        jmember <- 1
        while (jmember <= nmember[jmod]) {
          replace_values[["MEMBER_NUMBER"]] <- sprintf(paste("%.", (nmember[jmod] %/% 10) + 1, "i", sep = ''), jmember - 1)
          work_piece <- list(filename = .ConfigReplaceVariablesInString(exp[[jmod]][['path']], replace_values),
                             namevar = namevar, indices = c(1, jmember, jsdate, jmod), 
                             nmember = nmember[jmod], leadtimes = leadtimes, mask = maskmod[[jmod]],
                             is_file_per_dataset = FALSE, dimnames = exp[[jmod]][['dimnames']],
                             var_limits = c(mod_var_min, mod_var_max), remapcells = remapcells)
          exp_work_pieces <- c(exp_work_pieces, list(work_piece))
          jmember <- jmember + 1
        }
      } else {
        work_piece <- list(filename = .ConfigReplaceVariablesInString(exp[[jmod]][['path']], replace_values),
                           namevar = namevar, indices = c(1, 1, jsdate, jmod), 
                           nmember = nmember[jmod], leadtimes = leadtimes, mask = maskmod[[jmod]],
                           is_file_per_dataset = FALSE, dimnames = exp[[jmod]][['dimnames']],
                           var_limits = c(mod_var_min, mod_var_max), remapcells = remapcells)
        exp_work_pieces <- c(exp_work_pieces, list(work_piece))
      }
      jsdate <- jsdate + 1
    }
    jmod <- jmod + 1
  }
  if (dims2define && length(exp) > 0) {
    cat("! Warning: no data found in file system for any experimental dataset.\n")
  }

  dims <- dim_exp[na.omit(match(c('dataset', 'member', 'sdate', 'time', 'lat', 'lon'), names(dim_exp)))]
  if (is.null(dims[['member']]) || any(is.na(unlist(dims))) || any(unlist(dims) == 0)) {
    dims <- 0
    dim_exp <- NULL
  }
  if (!silent) {
    message <- "* Success. Detected dimensions of experimental data: "
    cat(paste0(message, paste(unlist(dims), collapse = ', '), '\n'))
    cat("* Fetching first observational files to work out 'var_obs' size...\n")
  }

  # If there are no experiments to load we need to choose a number of time steps
  # to load from observational datasets. We load from the first start date to 
  # the current date.
  if (is.null(exp) || dims == 0) {
    if (is.null(leadtimemax)) {
      cat("! Warning: loading observations only and no 'leadtimemax' specified. Data will be loaded from each starting date to current time.\n")
      diff <- Sys.time() - as.POSIXct(paste(substr(sdates[1], 1, 4), '-',
              substr(sdates[1], 5, 6), '-', substr(sdates[1], 7, 8), sep=''))
      if (storefreq == 'monthly') { 
        leadtimemax <- as.integer(diff/30) 
      } else {
        leadtimemax <- as.integer(diff)
      }
    }
    if (is.null(nleadtime)) {
      nleadtime <- leadtimemax
    }
    leadtimes <- seq(leadtimemin, leadtimemax, sampleperiod)
  }
  
  # Now we start iterating over observations. We try to find the output matrix
  # dimensions and we build anyway the work pieces corresponding to the observational
  # data that time-corresponds the experimental data or the time-steps until the
  # current date if no experimental datasets were specified.
  dataset_type <- 'obs'
  dim_obs <- NULL
  dims2define <- TRUE
  lat_indices <- lon_indices <- NULL
  obs_work_pieces <- list()
  is_file_per_dataset_obs <- rep(FALSE, nobs)
  is_file_per_member_obs <- rep(FALSE, nobs)
  jobs <- 1
  while (jobs <= nobs) {
    tags_to_find <- c('MONTH', 'DAY', 'YEAR', 'MEMBER_NUMBER')
    position_of_tags <- na.omit(match(tags_to_find, names(replace_values)))
    if (length(position_of_tags) > 0) {
      quasi_final_path <- .ConfigReplaceVariablesInString(obs[[jobs]][['path']], 
                            replace_values[-position_of_tags], TRUE)
    } else {
      quasi_final_path <- .ConfigReplaceVariablesInString(obs[[jobs]][['path']], 
                            replace_values, TRUE)
    }
    is_file_per_dataset_obs[jobs] <- !any(sapply(c("$MONTH$", "$DAY$", "$YEAR$"), 
                                          grepl, quasi_final_path, fixed = TRUE))
    is_file_per_member_obs[jobs] <- grepl("$MEMBER_NUMBER$", quasi_final_path, fixed = TRUE)
    replace_values[["OBS_NAME"]] <- obs[[jobs]][['name']]
    replace_values[["NC_VAR_NAME"]] <- obs[[jobs]][['nc_var_name']]
    namevar <- .ConfigReplaceVariablesInString(obs[[jobs]][['nc_var_name']], replace_values)
    replace_values[["SUFFIX"]] <- obs[[jobs]][['suffix']]
    if (is.null(varmin)) {
      obs_var_min <- as.numeric(.ConfigReplaceVariablesInString(obs[[jobs]][['var_min']], replace_values))
    } else {
      obs_var_min <- varmin
    }
    if (is.null(varmax)) {
      obs_var_max <- as.numeric(.ConfigReplaceVariablesInString(obs[[jobs]][['var_max']], replace_values))
    } else {
      obs_var_max <- varmax
    }
    # This file format (file per whole dataset) is only supported in observations.
    # However a file per whole dataset experiment could be seen as a file per
    # member/ensemble experiment with a single start date, so still loadable.
    # Nonetheless file per whole dataset observational files do not need to contain
    # a year and month in the filename, the time correspondance relies on the 
    # month and years associated to each timestep inside the NetCDF file.
    # So file per whole dataset experiments need to have a start date in the filename.
    if (is_file_per_dataset_obs[jobs]) {
      ## TODO: Open file-per-dataset-files only once.
      if (dims2define) {
        work_piece <- list(dataset_type = dataset_type,
                           filename = .ConfigReplaceVariablesInString(obs[[jobs]][['path']], replace_values),
                           namevar = namevar, grid = grid, remap = remap, remapcells = remapcells,
                           is_file_per_member = is_file_per_member_obs[jobs],
                           is_file_per_dataset = is_file_per_dataset_obs[jobs],
                           lon_limits = c(lonmin, lonmax), 
                           lat_limits = c(latmin, latmax), dimnames = obs[[jobs]][['dimnames']],
                           single_dataset = single_dataset)
        found_data <- .LoadDataFile(work_piece, explore_dims = TRUE, silent = silent)
        found_dims <- found_data$dims
        var_long_name <- found_data$var_long_name
        units <- found_data$units
        if (!is.null(found_dims)) {
          is_2d_var <- found_data$is_2d_var
          if (!is_2d_var && (output != 'areave')) {
            cat(paste("! Warning: '", output, "' output format not allowed when loading global mean variables. Forcing to 'areave'.\n",
                sep = ''))
            output <- 'areave'
          }
          if (output != 'areave' && is.null(grid)) {
            grid <- found_data$grid
          }
          if (is.null(nmemberobs)) {
            if (is.null(found_dims[['member']])) {
              cat("! Warning: loading observational data from a server but 'nmemberobs' not specified. Loading only one member.\n")
              nmemberobs <- rep(1, nobs)
            } else {
              nmemberobs <- rep(found_dims[['member']], nobs)
            }
          }
          if (is.null(dim_exp)) {
            longitudes <- found_dims[['lon']]
            latitudes <- found_dims[['lat']]
          }
          
          if (output == 'lon' || output == 'lonlat') {
            dim_obs[['lon']] <- length(longitudes)
          }
          if (output == 'lat' || output == 'lonlat') {
            dim_obs[['lat']] <- length(latitudes)
          }
          dim_obs[['time']] <- length(leadtimes)
          dim_obs[['member']] <- max(nmemberobs)
          dim_obs[['sdate']] <- nsdates
          dim_obs[['dataset']] <- nobs
          dims2define <- FALSE
        }
      }
      work_piece <- list(filename = .ConfigReplaceVariablesInString(obs[[jobs]][['path']], replace_values),
                         namevar = namevar, indices = c(1, 1, 1, jobs), 
                         nmember = nmemberobs[jobs], 
                         mask = maskobs[[jobs]], leadtimes = leadtimes, 
                         is_file_per_dataset = is_file_per_dataset_obs[jobs], 
                         startdates = sdates, dimnames = obs[[jobs]][['dimnames']],
                         var_limits = c(obs_var_min, obs_var_max), remapcells = remapcells)
      obs_work_pieces <- c(obs_work_pieces, list(work_piece))
    } else {
      jsdate <- 1
      while (jsdate <= nsdates) {
        replace_values[["START_DATE"]] <- sdates[jsdate]        
        sdate <- sdates[jsdate]

        if (storefreq == 'daily') {
          day <- substr(sdate, 7, 8)
          if (day == '') {
            day <- '01'
          }
          day <- as.integer(day)
          startdate <- as.POSIXct(paste(substr(sdate, 1, 4), '-',
                       substr(sdate, 5, 6), '-', day, ' 12:00:00', sep = '')) + 
                       (leadtimemin - 1) * 86400
          year <- as.integer(substr(startdate, 1, 4))
          month <- as.integer(substr(startdate, 6, 7))
        } else {
          month <- (as.integer(substr(sdate, 5, 6)) + leadtimemin - 2) %% 12 + 1
          year <- as.integer(substr(sdate, 1, 4)) + (as.integer(substr(sdate, 
                  5, 6)) + leadtimemin - 2) %/% 12
        }
        jleadtime <- 1
        while (jleadtime <= length(leadtimes)) {
          replace_values[["YEAR"]] <- paste(year, '', sep = '')
          replace_values[["MONTH"]] <- sprintf("%2.2i", month)
          if (storefreq == 'daily') {
            replace_values[["DAY"]] <- sprintf("%2.2i", day)
            days_in_month <- ifelse(LeapYear(year), 29, 28)
            days_in_month <- switch(paste(month, '', sep = ''), '1' = 31, 
                                    '3' = 31, '4' = 30, '5' = 31, '6' = 30, 
                                    '7' = 31, '8' = 31, '9' = 30, '10' = 31, 
                                    '11' = 30, '12' = 31, days_in_month)
            ## This condition must be fulfilled to put all the month time steps
            ## in the dimension of length nleadtimes. Otherwise it must be cut:
            #(length(leadtimes) - 1)*sampleperiod + 1 - (jleadtime - 1)*sampleperiod >= days_in_month - day + 1
            obs_file_indices <- seq(day, min(days_in_month, (length(leadtimes) - jleadtime) * sampleperiod + day), sampleperiod)
          } else {
            obs_file_indices <- 1
          }
          if (dims2define) {
            if (is_file_per_member_obs[jobs]) {
              replace_values[["MEMBER_NUMBER"]] <- '*'
            }
            work_piece <- list(dataset_type = dataset_type,
                               filename = .ConfigReplaceVariablesInString(obs[[jobs]][['path']], replace_values),
                               namevar = namevar, grid = grid, remap = remap, remapcells = remapcells,
                               is_file_per_member = is_file_per_member_obs[jobs],
                               is_file_per_dataset = is_file_per_dataset_obs[jobs],
                               lon_limits = c(lonmin, lonmax),
                               lat_limits = c(latmin, latmax), 
                               dimnames = obs[[jobs]][['dimnames']], single_dataset = single_dataset)
            found_data <- .LoadDataFile(work_piece, explore_dims = TRUE, silent = silent)
            found_dims <- found_data$dims
            var_long_name <- found_data$var_long_name
            units <- found_data$units
            if (!is.null(found_dims)) {
              is_2d_var <- found_data$is_2d_var
              if (!is_2d_var && (output != 'areave')) {
                cat(paste("! Warning: '", output, "' output format not allowed when loading global mean variables. Forcing to 'areave'\n.",
                    sep = ''))
                output <- 'areave'
              }
              if (output != 'areave' && is.null(grid)) {
                grid <- found_data$grid
              }
              if (is.null(nmemberobs)) {
                if (is.null(found_dims[['member']])) {
                  cat("! Warning: loading observational data from a server but 'nmemberobs' not specified. Loading only one member.\n")
                  nmemberobs <- rep(1, nobs)
                } else {
                  nmemberobs <- rep(found_dims[['member']], nobs)
                }
              }
              if (is.null(dim_exp)) {
                longitudes <- found_dims[['lon']]
                latitudes <- found_dims[['lat']]
              }
              
              if (output == 'lon' || output == 'lonlat') {
                dim_obs[['lon']] <- length(longitudes)
              }
              if (output == 'lat' || output == 'lonlat') {
                dim_obs[['lat']] <- length(latitudes)
              }
              dim_obs[['time']] <- length(leadtimes)
              dim_obs[['member']] <- max(nmemberobs)
              dim_obs[['sdate']] <- nsdates
              dim_obs[['dataset']] <- nobs
              dims2define <- FALSE
            }
          }
          if (is_file_per_member_obs[jobs]) {
            jmember <- 1
            while (jmember <= nmemberobs[jobs]) {
              replace_values[["MEMBER_NUMBER"]] <- sprintf(paste("%.", (nmemberobs[jobs] %/% 10) + 1, "i", sep = ''), jmember - 1)
              work_piece <- list(filename = .ConfigReplaceVariablesInString(obs[[jobs]][['path']], replace_values),
                                 namevar = namevar, indices = c(jleadtime, jmember, jsdate, jobs), 
                                 nmember = nmemberobs[jobs], leadtimes = obs_file_indices, 
                                 mask = maskobs[[jobs]], dimnames = obs[[jobs]][['dimnames']],
                                 is_file_per_dataset = is_file_per_dataset_obs[jobs], 
                                 var_limits = c(obs_var_min, obs_var_max), remapcells = remapcells)
              obs_work_pieces <- c(obs_work_pieces, list(work_piece))
              jmember <- jmember + 1
            }
          } else {
            work_piece <- list(filename = .ConfigReplaceVariablesInString(obs[[jobs]][['path']], replace_values),
                               namevar = namevar, indices = c(jleadtime, 1, jsdate, jobs), 
                               nmember = nmemberobs[jobs], leadtimes = obs_file_indices, 
                               mask = maskobs[[jobs]], dimnames = obs[[jobs]][['dimnames']],
                               is_file_per_dataset = is_file_per_dataset_obs[jobs], 
                               var_limits = c(obs_var_min, obs_var_max), remapcells)
            obs_work_pieces <- c(obs_work_pieces, list(work_piece))
          }
          
          if (storefreq == 'daily') {
            startdate <- startdate + 86400 * sampleperiod * length(obs_file_indices)
            year <- as.integer(substr(startdate, 1, 4))
            month <- as.integer(substr(startdate, 6, 7))
            day <- as.integer(substr(startdate, 9, 10))
          } else {
            month <- month + sampleperiod
            year <- year + (month - 1) %/% 12
            month <- (month - 1) %% 12 + 1
          }
          jleadtime <- jleadtime + length(obs_file_indices)
        }
        
        jsdate <- jsdate + 1
      }
    }
    jobs <- jobs + 1
  }
  if (dims2define && length(obs) > 0) {
    cat("! Warning: no data found in file system for any observational dataset.\n")
  }
  dims <- dim_obs[na.omit(match(c('dataset', 'member', 'sdate', 'time', 'lat', 'lon'), names(dim_obs)))]
  if (is.null(dims[['member']]) || any(is.na(unlist(dims))) || any(unlist(dims) == 0)) {
    dims <- 0
    dim_obs <- NULL
  }
  if (!silent) {
    message <- "* Success. Detected dimensions of observational data: "
    cat(paste0(message, paste(unlist(dims), collapse = ', '), '\n'))
  }

  if (!(is.null(dim_obs) && is.null(dim_exp))) {

  # We build two matrices in shared memory for the parallel processes to
  # store their results
  # These matrices will contain data arranged with the following
  # dimension order, to maintain data spacial locality during the 
  # parallel fetch:
  #   longitudes, latitudes, leadtimes, members, startdates, nmod/nobs
  # So [1, 1, 1, 1, 1, 1] will be next to [2, 1, 1, 1, 1, 1] in memory
  pointer_var_exp <- pointer_var_obs <- NULL
  if (!is.null(dim_exp) && (length(unlist(dim_exp)) == length(dim_exp)) && 
      !any(is.na(unlist(dim_exp))) && !any(unlist(dim_exp) == 0)) {
    var_exp <- big.matrix(nrow = prod(unlist(dim_exp)), ncol = 1)
    pointer_var_exp <- describe(var_exp)
  }
  if (!is.null(dim_obs) && (length(unlist(dim_obs)) == length(dim_obs)) && 
      !any(is.na(unlist(dim_obs))) && !any(unlist(dim_obs) == 0)) {
    var_obs <- big.matrix(nrow = prod(unlist(dim_obs)), ncol = 1)
    pointer_var_obs <- describe(var_obs)
  }
  if (is.null(nprocs)) {
    nprocs <- detectCores()
  }
  # We calculate the % of total progress that each work piece represents so 
  # that progress bar can be updated properly
  exp_work_piece_percent <- prod(dim_exp) / (prod(dim_obs) + prod(dim_exp))
  obs_work_piece_percent <- prod(dim_obs) / (prod(dim_obs) + prod(dim_exp))
  # Add some important extra fields in the work pieces before sending
  exp_work_pieces <- lapply(exp_work_pieces, function (x) c(x, list(dataset_type = 'exp', dims = dim_exp, out_pointer = pointer_var_exp)))###, progress_amount = exp_work_piece_progress)))
  obs_work_pieces <- lapply(obs_work_pieces, function (x) c(x, list(dataset_type = 'obs', dims = dim_obs, out_pointer = pointer_var_obs)))###, progress_amount = obs_work_piece_progress)))
  work_pieces <- c(exp_work_pieces, obs_work_pieces)
  # Calculate the progress %s that will be displayed and assign them to the 
  # appropriate work pieces
  if (length(work_pieces)/nprocs >= 2 && !silent) {
    if (length(work_pieces)/nprocs < 10) {
      amount <- 100/ceiling(length(work_pieces)/nprocs)
      reps <- ceiling(length(work_pieces)/nprocs)
    } else {
      amount <- 10
      reps <- 10
    }
    progress_steps <- rep(amount, reps)
    if (length(exp_work_pieces) == 0) {
      selected_exp_pieces <- c()
    } else if (length(exp_work_pieces) < floor(reps*exp_work_piece_percent) + 1) {
      selected_exp_pieces <- length(exp_work_pieces)
      progress_steps <- c(sum(head(progress_steps, 
                              floor(reps*exp_work_piece_percent))),
                          tail(progress_steps,
                               ceiling(reps*obs_work_piece_percent)))
    } else {
      selected_exp_pieces <- round(seq(1, length(exp_work_pieces), 
                                       length.out = floor(reps*exp_work_piece_percent) + 1))[-1]
    }
    if (length(obs_work_pieces) == 0) {
      selected_obs_pieces <- c()
    } else if (length(obs_work_pieces) < ceiling(reps*obs_work_piece_percent) + 1) {
      selected_obs_pieces <- length(obs_work_pieces)
      progress_steps <- c(head(progress_steps, 
                               floor(reps*exp_work_piece_percent)),
                          sum(tail(progress_steps,
                               ceiling(reps*obs_work_piece_percent))))
    } else {
      selected_obs_pieces <- round(seq(1, length(obs_work_pieces), 
                                       length.out = ceiling(reps*obs_work_piece_percent) + 1))[-1]
    }
    selected_pieces <- c(selected_exp_pieces, selected_obs_pieces + length(exp_work_pieces))
    progress_steps <- paste0(' + ', round(progress_steps, 2), '%')
    progress_message <- '* Progress: 0%'
  } else {
    progress_message <- ''
    selected_pieces <- NULL
  }
  piece_counter <- 1
  step_counter <- 1
  work_pieces <- lapply(work_pieces, 
                 function (x) {
                   wp <- c(x, list(is_2d_var = is_2d_var, grid = grid, remap = remap,  
                                   lon_limits = c(lonmin, lonmax), 
                                   lat_limits = c(latmin, latmax), 
                                   output = output, remapcells = remapcells,
                                   single_dataset = single_dataset))
                   if (piece_counter %in% selected_pieces) {
                     wp <- c(wp, list(progress_amount = progress_steps[step_counter]))
                     step_counter <<- step_counter + 1
                   }
                   piece_counter <<- piece_counter + 1
                   wp
                 })
  if (!silent) {
    cat(paste("* Will now proceed to read and process ", length(work_pieces), " data files:\n", sep = ''))
    if (length(work_pieces) < 30) {
      lapply(work_pieces, function (x) cat(paste("*   ", x[['filename']], '\n', sep = '')))
    } else {
      cat(paste("*   The list of files is long. You can check it after Load() finishes in the output '$source_files'.\n"))
    }
    if (length(dim_obs) == 0) {
      bytes_obs <- 0
      obs_dim_sizes <- '0'
    } else {
      bytes_obs <- prod(c(dim_obs, 8))
      obs_dim_sizes <- paste(na.omit(as.vector(dim_obs[c('dataset', 'member', 'sdate', 'time', 'lat', 'lon')])), collapse = ' x ')
    }
    if (length(dim_exp) == 0) {
      bytes_exp <- 0
      exp_dim_sizes <- '0'
    } else {
      bytes_exp <- prod(c(dim_exp, 8))
      exp_dim_sizes <- paste(na.omit(as.vector(dim_exp[c('dataset', 'member', 'sdate', 'time', 'lat', 'lon')])), collapse = ' x ')
    }
    cat(paste("* Total size of requested data: ", bytes_obs + bytes_exp, "bytes.\n"))
    cat(paste("*   - Experimental data:  (", exp_dim_sizes, ") x 8 bytes =", bytes_exp, "bytes.\n"))
    cat(paste("*   - Observational data: (", obs_dim_sizes, ") x 8 bytes =", bytes_obs, "bytes.\n"))
    cat(paste("* If size of requested data is close to or above the free shared RAM memory, R will crash.\n"))
  }
  # Build the cluster of processes that will do the work and dispatch work pieces.
  # The function .LoadDataFile is applied to each work package. This function will
  # open the data file, regrid if needed, trim (select time steps, longitudes, 
  # latitudes, members), apply the mask, compute and apply the weights if needed,
  # disable extreme values and store in the shared memory matrices. 
  if (nprocs == 1) {
    found_files <- lapply(work_pieces, .LoadDataFile, silent = silent)
  } else {
    cluster <- makeCluster(nprocs, outfile = "")
    # Open connections to keep track of progress
    ###range_progress_ports <- c(49000, 49999)
    ###progress_ports <- as.list(sample(range_progress_ports[2] - range_progress_ports[1], nprocs) + range_progress_ports[1])

    # Open from master side
    ###connection_set_up_job <- mcparallel({
    ###  progress_connections <- vector('list', length(progress_ports))
    ###  for (connection in 1:length(progress_ports)) {
    ###    attempts <- 0
    ###    max_attempts <- 3
    ###    while (is.null(progress_connections[[connection]]) && attempts < max_attempts) {
    ###      Sys.sleep(2)
    ###      suppressWarnings({
    ###        progress_connections[[connection]] <- try({
    ###          socketConnection(port = progress_ports[[connection]], open = 'w+b')
    ###        }, silent = TRUE)
    ###      })
    ###      if (!('sockconn' %in% class(progress_connections[[connection]]))) {
    ###        progress_connections[[connection]] <- NULL
    ###      }
    ###      attempts <- attempts + 1
    ###    }
    ###  }

      # And start polling the sockets and update progress bar
    ###  if (!any( lapply is.null!!! is.null(progress_connections))) {
    ###    progress <- 0.0
    ###    pb <- txtProgressBar(0, 1, style = 3)
    ###    stop_polling <- FALSE
    ###    attempts <- 0
    ###    max_attempts <- 3
    ###    while (progress < 0.999 && !stop_polling) {
    ###      Sys.sleep(3)
    ###      progress_obtained <- lapply(progress_connections, function(x) as.numeric(readBin(x, 'double')))
    ###      total_progress_obtained <- sum(unlist(progress_obtained))
    ###      if (total_progress_obtained > 0) {
    ###        progress <- progress + total_progress_obtained
    ###        setTxtProgressBar(pb, progress)
    ###        attempts <- 0
    ###      } else {
    ###        attempts <- attempts + 1
    ###        if (attempts >= max_attempts) {
    ###          stop_polling <- TRUE
    ###        }
    ###      }
    ###    }
    ###  }
    ###})

    # Open from the workers side
    ###open_connections <- clusterApply(cluster, progress_ports, 
    ###  function (x) {
    ###    progress_connection <<- NULL
    ###    progress_connection <<- try({
    ###      socketConnection(server = TRUE, port = x, open = 'w+b')
    ###    })
    ###    if ('sockconn' %in% class(progress_connection)) {
    ###      TRUE
    ###    } else {
    ###      progress_connection <<- NULL
    ###      FALSE
    ###    }
    ###  })

    ###if (!all(unlist(open_connections))) {
    ###  if (!silent) {
    ###    cat(paste("! Warning: failed to open connections in ports", process_track_ports[1], "to", process_track_ports[2], "to keep track of progress. Progress bar will not be displayed\n"))
    ###  }
    ###}

    if (!silent) {
      cat(paste("* Loading... This may take several minutes...\n", sep = ''))
      cat(progress_message)
    }
    # Send the heavy work to the workers
    work_errors <- try({
      found_files <- clusterApplyLB(cluster, work_pieces, .LoadDataFile, silent = silent)
    })
    stopCluster(cluster)
  }
  if (!silent) {
    if (progress_message != '') {
      cat("\n")
    }
    if (any(unlist(lapply(found_files, is.null)))) {
      if (sum(unlist(lapply(found_files, is.null))) < 30) {
        cat("! WARNING: The following files were not found in the file system. Filling with NA values instead.\n")
        lapply(work_pieces[which(unlist(lapply(found_files, is.null)))], function (x) cat(paste("*   ", x[['filename']], '\n', sep = '')))
      } else {
        cat("! WARNING: Some files were not found in the file system. The list is long. You can check it in the output '$not_found_files'. Filling with NA values instead.\n")
      }
    }
  }
  source_files <- unlist(found_files[which(!unlist(lapply(found_files, is.null)))])
  not_found_files <- unlist(lapply(work_pieces[which(unlist(lapply(found_files, is.null)))], '[[', 'filename'))

  } else {
    error_message <- "Error: No found files for any dataset. Check carefully the file patterns and correct either the pattern or the provided parameters:\n"
    if (!is.null(exp)) {
      lapply(exp, function (x) error_message <<- paste0(error_message, paste0(x[['path']], '\n')))
    }
    if (!is.null(obs)) {
      lapply(obs, function (x) error_message <<- paste0(error_message, paste0(x[['path']], '\n')))
    }
    stop(error_message)
  }

  })

  if (class(errors) == 'try-error') {
    invisible(list(load_parameters = load_parameters))
  } else {
    variable <- list()
    variable[['varName']] <- var
    variable[['level']] <- NULL
    attr(variable, 'is_standard') <- FALSE
    attr(variable, 'units') <- units
    attr(variable, 'longname') <- var_long_name
    attr(variable, 'daily_agg_cellfun') <- 'none'
    attr(variable, 'monthly_agg_cellfun') <- 'none'
    attr(variable, 'verification_time') <- 'none'
    
    if (is.null(var_exp)) {
      mod_data <- NULL
    } else { 
      dim_reorder <- length(dim_exp):1
      dim_reorder[2:3] <- dim_reorder[3:2]
      old_dims <- dim_exp
      dim_exp <- dim_exp[dim_reorder]
      mod_data <- aperm(array(bigmemory::as.matrix(var_exp), dim = old_dims), dim_reorder)
      attr(mod_data, 'dimensions') <- names(dim_exp)
    }

    if (is.null(var_obs)) {
      obs_data <- NULL
    } else {                     
      dim_reorder <- length(dim_obs):1
      dim_reorder[2:3] <- dim_reorder[3:2]
      old_dims <- dim_obs
      dim_obs <- dim_obs[dim_reorder]
      obs_data <- aperm(array(bigmemory::as.matrix(var_obs), dim = old_dims), dim_reorder)
      attr(obs_data, 'dimensions') <- names(dim_obs)
    }

    if (is.null(latitudes)) {
      lat <- 0
      attr(lat, 'cdo_grid_name') <- 'none'
    } else {
      lat <- latitudes
      attr(lat, 'cdo_grid_name') <- if (is.null(grid)) 'none' else grid
    }
    attr(lat, 'projection') <- 'none'

    if (is.null(longitudes)) {
      lon <- 0
      attr(lon, 'cdo_grid_name') <- 'none' 
    } else {
      lon <- longitudes
      attr(lon, 'cdo_grid_name') <- if (is.null(grid)) 'none' else grid
    }
    attr(lon, 'projection') <- 'none'

    dates <- list()
    dates[['start']] <- NULL
    dates[['end']] <- NULL

    models <- NULL
    if (length(exp) > 0 && !is.null(dim_exp)) {
      models <- list()
      for (jmod in 1:length(exp)) {
        models[[exp[[jmod]][['name']]]] <- list(
          members = paste0('Member_', 1:nmember[jmod]),
          source = if ((nchar(exp[[jmod]][['path']]) - 
                        nchar(gsub("/", "", exp[[jmod]][['path']])) > 2) &&
                       (length(sdates) > 1 && !is_file_per_member_exp[jmod])) {
                     parts <- strsplit(exp[[jmod]][['path']], '/')[[1]]
                     paste(parts[-length(parts)], sep = '', collapse = '/')
                   } else {
                     exp[[jmod]][['path']]
                   })
      }
    }

    observations <- NULL
    if (length(obs) > 0 && !is.null(dim_obs)) {
      observations <- list()
      for (jobs in 1:length(obs)) {
        observations[[obs[[jobs]][['name']]]] <- list(
          members = paste0('Member_', 1:nmemberobs[jobs]),
          source = if ((nchar(obs[[jobs]][['path']]) - 
                        nchar(gsub("/", "", obs[[jobs]][['path']])) > 2) &&
                       !is_file_per_dataset_obs[jobs]) {
                     parts <- strsplit(obs[[jobs]][['path']], '/')[[1]]
                     paste(parts[-length(parts)], sep = '', collapse = '/')
                   } else {
                     obs[[jobs]][['path']]
                   })
      }
    }

    # Before ending, the data is arranged in the common format, with the following
    # dimension order:
    #  nmod/nobs, members, startdates, leadtimes, latitudes, longitudes
    invisible(list(mod = mod_data,
                   obs = obs_data,
                   lon = lon,
                   lat = lat,
                   Variable = variable,
                   Datasets = list(exp = models, obs = observations),
                   Dates = dates,
                   InitializationDates = lapply(sdates, 
                     function (x) {
                       sink('/dev/null')
                       date <- print(as.POSIXct(as.Date(x, format = '%Y%m%d')))
                       sink()
                       date
                     }),
                   when = Sys.time(),
                   source_files = source_files,
                   not_found_files = not_found_files,
                   load_parameters = load_parameters))
  }
}
