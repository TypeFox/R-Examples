## Function to tell if a regexpr() match is a complete match to a specified name
.IsFullMatch <- function(x, name) {
  ifelse(x > 0 && attributes(x)$match.length == nchar(name), TRUE, FALSE)
}

.ConfigReplaceVariablesInString <- function(string, replace_values, allow_undefined_key_vars = FALSE) {
  # This function replaces all the occurrences of a variable in a string by 
  # their corresponding string stored in the replace_values.
  if (length(strsplit(string, "\\$")[[1]]) > 1) {
    parts <- strsplit(string, "\\$")[[1]]
    output <- ""
    i <- 0
    for (part in parts) {
      if (i %% 2 == 0) {
        output <- paste(output, part, sep = "")
      } else {
        if (part %in% names(replace_values)) {
          output <- paste(output, .ConfigReplaceVariablesInString(replace_values[[part]], replace_values, allow_undefined_key_vars), sep = "")
        } else if (allow_undefined_key_vars) {
          output <- paste0(output, "$", part, "$")
        } else {
          stop(paste('Error: The variable $', part, '$ was not defined in the configuration file.', sep = ''))
        }
      }
      i <- i + 1
    }
    output
  } else {
    string
  }
}

.LoadDataFile <- function(work_piece, explore_dims = FALSE, silent = FALSE) {
  # The purpose, working modes, inputs and outputs of this function are
  # explained in ?LoadDataFile
  #suppressPackageStartupMessages({library(ncdf4)})
  #suppressPackageStartupMessages({library(bigmemory)})
  #suppressPackageStartupMessages({library(plyr)})
  # Auxiliar function to convert array indices to lineal indices
  arrayIndex2VectorIndex <- function(indices, dims) {
    if (length(indices) > length(dims)) {
      stop("Error: indices do not match dimensions in arrayIndex2VectorIndex.")
    }
    position <- 1
    dims <- rev(dims)
    indices <- rev(indices)
    for (i in 1:length(indices)) {
      position <- position + (indices[i] - 1) * prod(dims[-c(1:i)])
    }
    position
  }
  .t2nlatlon <- function(t) {
    ## As seen in cdo's griddes.c: ntr2nlat()
    nlats <- (t * 3 + 1) / 2
    if ((nlats > 0) && (nlats - trunc(nlats) >= 0.5)) {
      nlats <- ceiling(nlats)
    } else {
      nlats <- round(nlats)
    }
    if (nlats %% 2 > 0) {
      nlats <- nlats + 1
    }
    ## As seen in cdo's griddes.c: compNlon(), and as specified in ECMWF
    nlons <- 2 * nlats
    keep_going <- TRUE
    while (keep_going) {
      n <- nlons
      if (n %% 8 == 0) n <- trunc(n / 8)
      while (n %% 6 == 0) n <- trunc(n / 6)
      while (n %% 5 == 0) n <- trunc(n / 5)
      while (n %% 4 == 0) n <- trunc(n / 4)
      while (n %% 3 == 0) n <- trunc(n / 3)
      if (n %% 2 == 0) n <- trunc(n / 2)
      if (n <= 8) {
        keep_going <- FALSE
      } else {
        nlons <- nlons + 2
        if (nlons > 9999) {
          stop("Error: pick another gaussian grid truncation. It doesn't fulfill the standards to apply FFT.")
        }
      }
    }
    c(nlats, nlons)
  }
  
  .nlat2t <- function(nlats) {
    trunc((nlats * 2 - 1) / 3)
  }

  found_file <- NULL
  dims <- NULL
  grid_name <- units <- var_long_name <- is_2d_var <- NULL

  filename <- work_piece[['filename']]
  namevar <- work_piece[['namevar']]
  output <- work_piece[['output']]
  # The names of all data files in the directory of the repository that match 
  # the pattern are obtained.
  if (length(grep("^http", filename)) > 0) {
    is_url <- TRUE
    files <- filename
    ## TODO: Check that the user is not using shell globbing exps.
  } else {
    is_url <- FALSE
    files <- Sys.glob(filename)
  }

  # If we don't find any, we leave the flag 'found_file' with a NULL value.
  if (length(files) > 0) {
    # The first file that matches the pattern is chosen and read.
    filename <- files[length(files)]
    filein <- filename
    found_file <- filename
    mask <- work_piece[['mask']]

    if (!silent) {
      if (explore_dims) {
        cat(paste("* Exploring dimensions...", filename, '\n'))
      }
      ##} else {
      ##  cat(paste("* Reading & processing data...", filename, '\n'))
      ##}
    }

    # We will fill in 'expected_dims' with the names of the expected dimensions of
    # the data array we'll retrieve from the file.
    expected_dims <- NULL
    remap_needed <- FALSE
    # But first we open the file and work out whether the requested variable is 2d
    fnc <- nc_open(filein)
    if (!(namevar %in% names(fnc$var))) {
      stop(paste("Error: The variable", namevar, "is not defined in the file", filename))
    }
    var_long_name <- fnc$var[[namevar]]$longname
    units <- fnc$var[[namevar]]$units
    if (is.null(work_piece[['is_2d_var']])) {
      is_2d_var <- all(c(work_piece[['dimnames']][['lon']], 
                         work_piece[['dimnames']][['lat']]) %in%
                       unlist(lapply(fnc$var[[namevar]][['dim']], 
                                     '[[', 'name')))
    } else {
      is_2d_var <- work_piece[['is_2d_var']]
    }
    if ((is_2d_var || work_piece[['is_file_per_dataset']]) && (Sys.which("cdo")[[1]] == "")) {
      stop("Error: CDO libraries not available")
    }
    # If the variable to load is 2-d, we need to determine whether:
    #  - interpolation is needed
    #  - subsetting is requested
    if (is_2d_var) {
      ## We read the longitudes and latitudes from the file.
      lon <- ncvar_get(fnc, work_piece[['dimnames']][['lon']])
      lat <- ncvar_get(fnc, work_piece[['dimnames']][['lat']])
      # If a common grid is requested or we are exploring the file dimensions
      # we need to read the grid type and size of the file to finally work out the 
      # CDO grid name.
      if (!is.null(work_piece[['grid']]) || explore_dims) {
        # Here we read the grid type and its number of longitudes and latitudes
        file_info <- system(paste('cdo -s griddes', filein, '2> /dev/null'), intern = TRUE)
        grids_positions <- grep('# gridID', file_info)
        grids_first_lines <- grids_positions + 2
        grids_last_lines <- c((grids_positions - 2)[-1], length(file_info))
        grids_info <- as.list(1:length(grids_positions))
        grids_info <- lapply(grids_info, function (x) file_info[grids_first_lines[x]:grids_last_lines[x]])
        grids_info <- lapply(grids_info, function (x) gsub("  *", " ", x))
        grids_info <- lapply(grids_info, function (x) gsub("^ | $", "", x))
        grids_info <- lapply(grids_info, function (x) unlist(strsplit(x, " | = ")))
        grids_types <- unlist(lapply(grids_info, function (x) x[grep('gridtype', x) + 1]))
        grids_matches <- unlist(lapply(grids_info, function (x) {
          nlons <- if (length(grep('xsize', x)) > 0) {
                     as.integer(x[grep('xsize', x) + 1])
                   } else {
                     NA
                   }
          nlats <- if (length(grep('ysize', x)) > 0) {
                    as.integer(x[grep('ysize', x) + 1])
                  } else {
                    NA
                  }
          if (identical(nlons, length(lon)) && 
              identical(nlats, length(lat))) {
            TRUE
          } else {
            FALSE
          }
        }))
        grids_matches <- grids_matches[which(grids_types %in% c('gaussian', 'lonlat'))]
        grids_info <- grids_info[which(grids_types %in% c('gaussian', 'lonlat'))]
        grids_types <- grids_types[which(grids_types %in% c('gaussian', 'lonlat'))]
        if (length(grids_matches) == 0) {
          stop("Error: Only 'gaussian' and 'lonlat' grids supported. See e.g: cdo sinfo ", filename)
        }
        if (sum(grids_matches) > 1) {
          if ((all(grids_types[which(grids_matches)] == 'gaussian') || 
               all(grids_types[which(grids_matches)] == 'lonlat')) && 
               all(unlist(lapply(grids_info[which(grids_matches)], identical, 
                                 grids_info[which(grids_matches)][[1]])))) {
            grid_type <- grids_types[which(grids_matches)][1]
          } else {
            stop("Error: Load() can't disambiguate: More than one lonlat/gaussian grids with the same size as the requested variable defined in ", filename)
          }
        } else {
          grid_type <- grids_types[which(grids_matches)]
        }
        grid_lons <- length(lon)
        grid_lats <- length(lat)
        # Convert to CDO grid name as seen in cdo's griddes.c: nlat2ntr()
        if (grid_type == 'lonlat') {
          grid_name <- paste0('r', grid_lons, 'x', grid_lats)
        } else {
          grid_name <- paste0('t', .nlat2t(grid_lats), 'grid')
        }
      }
      # If a common grid is requested, we will also calculate its size which we will use
      # later on.
      if (!is.null(work_piece[['grid']])) {
        # Now we calculate the common grid type and its lons and lats
        if (length(grep('^t\\d{1,+}grid$', work_piece[['grid']])) > 0) {
          common_grid_type <- 'gaussian'
          common_grid_res <- as.integer(strsplit(work_piece[['grid']], '[^0-9]{1,+}')[[1]][2])
          nlonlat <- .t2nlatlon(common_grid_res)
          common_grid_lats <- nlonlat[1]
          common_grid_lons <- nlonlat[2]
        } else if (length(grep('^r\\d{1,+}x\\d{1,+}$', work_piece[['grid']])) > 0) {
          common_grid_type <- 'lonlat'
          common_grid_lons <- as.integer(strsplit(work_piece[['grid']], '[^0-9]{1,+}')[[1]][2])
          common_grid_lats <- as.integer(strsplit(work_piece[['grid']], '[^0-9]{1,+}')[[1]][3])
        } else {
          stop("Error: Only supported grid types in parameter 'grid' are t<RES>grid and r<NX>x<NY>")
        }
      } else {
        ## If no 'grid' is specified, there is no common grid.
        ## But these variables are filled in for consistency in the code.
        common_grid_lons <- length(lon)
        common_grid_lats <- length(lat)
      }
      first_common_grid_lon <- 0
      last_common_grid_lon <- 360 - 360/common_grid_lons
      ## This is not true for gaussian grids or for some regular grids, but 
      ## is a safe estimation
      first_common_grid_lat <- -90
      last_common_grid_lat <- 90
      # And finally determine whether interpolation is needed or not
      remove_shift <- FALSE
      if (!is.null(work_piece[['grid']])) {
        if ((grid_lons != common_grid_lons) || 
            (grid_lats != common_grid_lats) || 
            (grid_type != common_grid_type) ||
            ((lon[1] != first_common_grid_lon) 
             && !work_piece[['single_dataset']])) {
          if (grid_lons == common_grid_lons && grid_lats == common_grid_lats &&
              grid_type == common_grid_type && lon[1] != first_common_grid_lon &&
              !work_piece[['single_dataset']]) {
            remove_shift <- TRUE
          }
          remap_needed <- TRUE
          common_grid_name <- work_piece[['grid']]
        }
      } else if ((lon[1] != first_common_grid_lon) && explore_dims && 
                 !work_piece[['single_dataset']]) {
        remap_needed <- TRUE
        common_grid_name <- grid_name
        remove_shift <- TRUE
      }
      if (remove_shift && !explore_dims) {
        if (!is.null(work_piece[['progress_amount']])) {
          cat("\n")
        }
        cat(paste0("! Warning: The dataset with index ", 
            tail(work_piece[['indices']], 1), " in '", 
            work_piece[['dataset_type']], "' doesn't start at longitude 0 and will be re-interpolated in order to align its longitudes with the standard CDO grids definable with the names 't<RES>grid' or 'r<NX>x<NY>', which are by definition starting at the longitude 0.\n"))
        if (!is.null(mask)) {
          cat(paste0("! Warning: A mask was provided for the dataset with index ",    
              tail(work_piece[['indices']], 1), " in '",
              work_piece[['dataset_type']], "'. This dataset has been re-interpolated to align its longitudes to start at 0. You must re-interpolate the corresponding mask to align its longitudes to start at 0 as well, if you haven't done so yet. Running cdo remapcon,", common_grid_name, " original_mask_file.nc new_mask_file.nc will fix it.\n"))
        }
      }
      # Now calculate if the user requests for a lonlat subset or for the 
      # entire field
      lonmin <- work_piece[['lon_limits']][1]
      lonmax <- work_piece[['lon_limits']][2]
      latmin <- work_piece[['lat_limits']][1]
      latmax <- work_piece[['lat_limits']][2]
      lonlat_subsetting_requested <- FALSE
      if (lonmin <= lonmax) {
        if ((lonmin > first_common_grid_lon) || (lonmax < last_common_grid_lon)) {
          lonlat_subsetting_requested <- TRUE
        }
      } else {
        if ((lonmin - lonmax) > 360/common_grid_lons) {
          lonlat_subsetting_requested <- TRUE
        } else {
          gap_width <- floor(lonmin / (360/common_grid_lons)) - 
                       floor(lonmax / (360/common_grid_lons))
          if (gap_width > 0) { 
            if (!(gap_width == 1 && (lonmin %% (360/common_grid_lons) == 0) && 
                  (lonmax %% (360/common_grid_lons) == 0))) {
              lonlat_subsetting_requested <- TRUE
            }
          }
        }
      }
      if ((latmin > first_common_grid_lat) || (latmax < last_common_grid_lat)) {
        lonlat_subsetting_requested <- TRUE
      }

      # When remap is needed but no subsetting, the file is copied locally
      # so that cdo works faster, and then interpolated.
      # Otherwise the file is kept as is and the subset will have to be 
      # interpolated still.
      if (!lonlat_subsetting_requested && remap_needed) {
        nc_close(fnc)
        filecopy <- tempfile(pattern = "load", fileext = ".nc")
        file.copy(filein, filecopy)
        filein <- tempfile(pattern = "loadRegridded", fileext = ".nc")
        system(paste0("cdo -s ", work_piece[['remap']], ",", 
                      common_grid_name, 
                      " -selname,", namevar, " ", filecopy, " ", filein, 
                      " 2>/dev/null", sep = ""))
        file.remove(filecopy)
        work_piece[['dimnames']][['lon']] <- 'lon'
        work_piece[['dimnames']][['lat']] <- 'lat'
        fnc <- nc_open(filein)
        lon <- ncvar_get(fnc, work_piece[['dimnames']][['lon']])
        lat <- ncvar_get(fnc, work_piece[['dimnames']][['lat']])
      }

      # Read and check also the mask
      if (!is.null(mask)) {
        ###mask_file <- tempfile(pattern = 'loadMask', fileext = '.nc')
        if (is.list(mask)) {
          if (!file.exists(mask[['path']])) {
            stop(paste("Error: Couldn't find the mask file", mask[['path']]))
          }
          mask_file <- mask[['path']]
          ###file.copy(work_piece[['mask']][['path']], mask_file)
          fnc_mask <- nc_open(mask_file)
          vars_in_mask <- sapply(fnc_mask$var, '[[', 'name')
          if ('nc_var_name' %in% names(mask)) {
            if (!(mask[['nc_var_name']] %in% 
                  vars_in_mask)) {
              stop(paste("Error: couldn't find variable", mask[['nc_var_name']], 
                         "in the mask file", mask[['path']]))
            }
          } else {
            if (length(vars_in_mask) != 1) {
              stop(paste("Error: one and only one non-coordinate variable should be defined in the mask file", 
                   mask[['path']], "if the component 'nc_var_name' is not specified. Currently found: ", 
                   paste(vars_in_mask, collapse = ', '), "."))
            } else {
              mask[['nc_var_name']] <- vars_in_mask
            }
          }
          if (sum(fnc_mask$var[[mask[['nc_var_name']]]]$size > 1) != 2) {
            stop(paste0("Error: the variable '", 
                 mask[['nc_var_name']], 
                 "' must be defined only over the dimensions '", 
                 work_piece[['dimnames']][['lon']], "' and '", 
                 work_piece[['dimnames']][['lat']], 
                 "' in the mask file ", 
                 mask[['path']]))
          }
          mask <- ncvar_get(fnc_mask, mask[['nc_var_name']], collapse_degen = TRUE)
          nc_close(fnc_mask)
        ###  mask_lon <- ncvar_get(fnc_mask, work_piece[['dimnames']][['lon']])
        ###  mask_lat <- ncvar_get(fnc_mask, work_piece[['dimnames']][['lat']])
        ###} else {
        ###  dim_longitudes <- ncdim_def(work_piece[['dimnames']][['lon']], "degrees_east", lon)
        ###  dim_latitudes <- ncdim_def(work_piece[['dimnames']][['lat']], "degrees_north", lat)
        ###  ncdf_var <- ncvar_def('LSM', "", list(dim_longitudes, dim_latitudes), NA, 'double')
        ###  fnc_mask <- nc_create(mask_file, list(ncdf_var))
        ###  ncvar_put(fnc_mask, ncdf_var, work_piece[['mask']])
        ###  nc_close(fnc_mask)
        ###  fnc_mask <- nc_open(mask_file)
        ###  work_piece[['mask']] <- list(path = mask_file, nc_var_name = 'LSM')
        ###  mask_lon <- lon
        ###  mask_lat <- lat
        ###}
      ###}
        ### Now ready to check that the mask is right
        ##if (!(lonlat_subsetting_requested && remap_needed)) {
        ###  if ((dim(mask)[2] != length(lon)) || (dim(mask)[1] != length(lat))) {
        ###    stop(paste("Error: the mask of the dataset with index ", tail(work_piece[['indices']], 1), " in '", work_piece[['dataset_type']], "' is wrong. It must be on the common grid if the selected output type is 'lonlat', 'lon' or 'lat', or 'areave' and 'grid' has been specified. It must be on the grid of the corresponding dataset if the selected output type is 'areave' and no 'grid' has been specified. For more information check ?Load and see help on parameters 'grid', 'maskmod' and 'maskobs'.", sep = ""))
        ###  }
          ###if (!(identical(mask_lon, lon) && identical(mask_lat, lat))) {
          ###  stop(paste0("Error: the longitudes and latitudes in the masks must be identical to the ones in the corresponding data files if output = 'areave' or, if the selected output is 'lon', 'lat' or 'lonlat', the longitudes in the mask file must start by 0 and the latitudes must be ordered from highest to lowest. See\n  ", 
          ###     work_piece[['mask']][['path']], " and ", filein))
          ###}
        }
      }

      lon_indices <- 1:length(lon)
      if (!(lonlat_subsetting_requested && remap_needed)) {
        lon[which(lon < 0)] <- lon[which(lon < 0)] + 360
      }
      if (lonmax >= lonmin) {
        lon_indices <- lon_indices[which(((lon %% 360) >= lonmin) & ((lon %% 360) <= lonmax))]
      } else if (!remap_needed) {
        lon_indices <- lon_indices[which(((lon %% 360) <= lonmax) | ((lon %% 360) >= lonmin))]
      }
      lat_indices <- which(lat >= latmin & lat <= latmax)
      ## In most of the cases the latitudes are ordered from -90 to 90. 
      ## We will reorder them to be in the order from 90 to -90, so mostly 
      ## always the latitudes are reordered.
      ## TODO: This could be avoided in future.
      if (lat[1] < lat[length(lat)]) {
        lat_indices <- lat_indices[length(lat_indices):1]
      }
      if (!is.null(mask) && !(lonlat_subsetting_requested && remap_needed)) {
        if ((dim(mask)[1] != length(lon)) || (dim(mask)[2] != length(lat))) {
          stop(paste("Error: the mask of the dataset with index ", tail(work_piece[['indices']], 1), " in '", work_piece[['dataset_type']], "' is wrong. It must be on the common grid if the selected output type is 'lonlat', 'lon' or 'lat', or 'areave' and 'grid' has been specified. It must be on the grid of the corresponding dataset if the selected output type is 'areave' and no 'grid' has been specified. For more information check ?Load and see help on parameters 'grid', 'maskmod' and 'maskobs'.", sep = ""))
        }
        mask <- mask[lon_indices, lat_indices]
      }
      ## If the user requests subsetting, we must extend the lon and lat limits if possible
      ## so that the interpolation after is done properly
      maximum_extra_points <- work_piece[['remapcells']]
      if (lonlat_subsetting_requested && remap_needed) {
        if ((maximum_extra_points > (head(lon_indices, 1) - 1)) ||
            (maximum_extra_points > (length(lon) - tail(lon_indices, 1)))) {
          ## if the requested number of points goes beyond the left or right
          ## sides of the map, we need to take the entire map so that the 
          ## interpolation works properly
          lon_indices <- 1:length(lon)
        } else {
          extra_points <- min(maximum_extra_points, head(lon_indices, 1) - 1)
          if (extra_points > 0) {
            lon_indices <- c((head(lon_indices, 1) - extra_points):(head(lon_indices, 1) - 1), lon_indices)
          }
          extra_points <- min(maximum_extra_points, length(lon) - tail(lon_indices, 1))
          if (extra_points > 0) {
            lon_indices <- c(lon_indices, (tail(lon_indices, 1) + 1):(tail(lon_indices, 1) + extra_points))
          }
        }
        min_lat_ind <- min(lat_indices)
        max_lat_ind <- max(lat_indices)
        extra_points <- min(maximum_extra_points, min_lat_ind - 1)
        if (extra_points > 0) {
          if (lat[1] < tail(lat, 1)) {
            lat_indices <- c(lat_indices, (min_lat_ind - 1):(min_lat_ind - extra_points))
          } else {
            lat_indices <- c((min_lat_ind - extra_points):(min_lat_ind - 1), lat_indices)
          }
        }
        extra_points <- min(maximum_extra_points, length(lat) - max_lat_ind)
        if (extra_points > 0) {
          if (lat[1] < tail(lat, 1)) {
            lat_indices <- c((max_lat_ind + extra_points):(max_lat_ind + 1), lat_indices)
          } else {
            lat_indices <- c(lat_indices, (max_lat_ind + 1):(max_lat_ind + extra_points))
          }
        }
      }
      lon <- lon[lon_indices]
      lat <- lat[lat_indices]
      expected_dims <- c(work_piece[['dimnames']][['lon']],
                         work_piece[['dimnames']][['lat']])
    } else {
      lon <- 0
      lat <- 0
    }
    # We keep on filling the expected dimensions
    var_dimnames <- unlist(lapply(fnc$var[[namevar]][['dim']], '[[', 'name'))
    nmemb <- nltime <- NULL
    ## Sometimes CDO renames 'members' dimension to 'lev'
    old_members_dimname <- NULL
    if (('lev' %in% var_dimnames) && !(work_piece[['dimnames']][['member']] %in% var_dimnames)) {
      old_members_dimname <- work_piece[['dimnames']][['member']]
      work_piece[['dimnames']][['member']] <- 'lev'
    }
    if (work_piece[['dimnames']][['member']] %in% var_dimnames) {
      nmemb <- fnc$var[[namevar]][['dim']][[match(work_piece[['dimnames']][['member']], var_dimnames)]]$len
      expected_dims <- c(expected_dims, work_piece[['dimnames']][['member']])
    } else {
      nmemb <- 1
    }
    if (length(expected_dims) > 0) {
      dim_matches <- match(expected_dims, var_dimnames)
      if (any(is.na(dim_matches))) {
        if (!is.null(old_members_dimname)) {
          expected_dims[which(expected_dims == 'lev')] <- old_members_dimname
        }
        stop(paste("Error: the expected dimension(s)", 
                   paste(expected_dims[which(is.na(dim_matches))], collapse = ', '), 
                   "were not found in", filename))
      }
      time_dimname <- var_dimnames[-dim_matches]
    } else {
      time_dimname <- var_dimnames
    }
    if (length(time_dimname) > 0) {
      if (length(time_dimname) == 1) {
        nltime <- fnc$var[[namevar]][['dim']][[match(time_dimname, var_dimnames)]]$len
        expected_dims <- c(expected_dims, time_dimname)
        dim_matches <- match(expected_dims, var_dimnames)
      } else {
        if (!is.null(old_members_dimname)) {
          expected_dims[which(expected_dims == 'lev')] <- old_members_dimname
        }
        stop(paste("Error: the variable", namevar, 
                   "is defined over more dimensions than the expected (", 
                   paste(c(expected_dims, 'time'), collapse = ', '), 
                   "). It could also be that the members dimension in 'dimnames' or in the configuration file is incorrect. If not, it could also be that the members dimension is named incorrectly. In that case, either rename the dimension in the file or adjust Load() to recognize this name with the parameter 'dimnames'. See file", filename))
      }
    } else {
      nltime <- 1
    }

    # Now we must retrieve the data from the file, but only the asked indices.
    # So we build up the indices to retrieve.
    # Longitudes or latitudes have been retrieved already.
    if (explore_dims) {
      # If we're exploring the file we only want one time step from one member, 
      # to regrid it and work out the number of longitudes and latitudes.
      # We don't need more.
      members <- 1
      ltimes_list <- list(c(1))
    } else {
      # The data is arranged in the array 'tmp' with the dimensions in a 
      # common order:
      #   1) Longitudes 
      #   2) Latitudes
      #   3) Members (even if is not a file per member experiment)
      #   4) Lead-times
      if (work_piece[['is_file_per_dataset']]) {
        time_indices <- 1:nltime
        mons <- strsplit(system(paste('cdo showmon ', filein, 
                         ' 2>/dev/null'), intern = TRUE), split = ' ')
        years <- strsplit(system(paste('cdo showyear ', filein, 
                          ' 2>/dev/null'), intern = TRUE), split = ' ')
        mons <- as.integer(mons[[1]][which(mons[[1]] != "")])
        years <- as.integer(years[[1]][which(years[[1]] != "")])
        time_indices <- ts(time_indices, start = c(years[1], mons[1]), 
                           end = c(years[length(years)], mons[length(mons)]),
                           frequency = 12)
        ltimes_list <- list()
        for (sdate in work_piece[['startdates']]) {
          selected_time_indices <- window(time_indices, start = c(as.integer(
                                   substr(sdate, 1, 4)), as.integer(substr(sdate, 5, 6))), 
                                   end = c(3000, 12), frequency = 12, extend = TRUE)
          selected_time_indices <- selected_time_indices[work_piece[['leadtimes']]]
          ltimes_list <- c(ltimes_list, list(selected_time_indices))
        }
      } else {
        ltimes <- work_piece[['leadtimes']]
        #if (work_piece[['dataset_type']] == 'exp') {
          ltimes_list <- list(ltimes[which(ltimes <= nltime)])
        #}
      }
      ## TODO: Put, when reading matrices, this kind of warnings
      #  if (nmember < nmemb) {
      #    cat("Warning:
      members <- 1:work_piece[['nmember']]
      members <- members[which(members <= nmemb)]
    }

    # Now, for each list of leadtimes to load (usually only one list with all leadtimes), 
    # we'll join the indices and retrieve data
    found_disordered_dims <- FALSE
    for (ltimes in ltimes_list) {
      if (is_2d_var) {
        start <- c(min(lon_indices), min(lat_indices))
        end <- c(max(lon_indices), max(lat_indices))
        if (lonlat_subsetting_requested && remap_needed) {
          subset_indices <- list(min(lon_indices):max(lon_indices) - min(lon_indices) + 1,
                                 lat_indices - min(lat_indices) + 1)
          dim_longitudes <- ncdim_def(work_piece[['dimnames']][['lon']], "degrees_east", lon)
          dim_latitudes <- ncdim_def(work_piece[['dimnames']][['lat']], "degrees_north", lat)
          ncdf_dims <- list(dim_longitudes, dim_latitudes)
        } else {
          subset_indices <- list(lon_indices - min(lon_indices) + 1,
                                 lat_indices - min(lat_indices) + 1)
          ncdf_dims <- list()
        }
        final_dims <- c(length(subset_indices[[1]]), length(subset_indices[[2]]), 1, 1)
      } else {
        start <- end <- c()
        subset_indices <- list()
        ncdf_dims <- list()
        final_dims <- c(1, 1, 1, 1)
      }
      
      if (work_piece[['dimnames']][['member']] %in% expected_dims) {
        start <- c(start, head(members, 1))
        end <- c(end, tail(members, 1))
        subset_indices <- c(subset_indices, list(members - head(members, 1) + 1))
        dim_members <- ncdim_def(work_piece[['dimnames']][['member']], "", members)
        ncdf_dims <- c(ncdf_dims, list(dim_members))
        final_dims[3] <- length(members)
      }
      if (time_dimname %in% expected_dims) {
        if (any(!is.na(ltimes))) {
          start <- c(start, head(ltimes[which(!is.na(ltimes))], 1))
          end <- c(end, tail(ltimes[which(!is.na(ltimes))], 1))
          subset_indices <- c(subset_indices, list(ltimes - head(ltimes[which(!is.na(ltimes))], 1) + 1))
        } else {
          start <- c(start, NA)
          end <- c(end, NA)
          subset_indices <- c(subset_indices, list(ltimes))
        }
        dim_time <- ncdim_def(time_dimname, "", 1:length(ltimes), unlim = TRUE)
        ncdf_dims <- c(ncdf_dims, list(dim_time))
        final_dims[4] <- length(ltimes)
      }
      count <- end - start + 1
      start <- start[dim_matches]
      count <- count[dim_matches]
      subset_indices <- subset_indices[dim_matches]
      # Now that we have the indices to retrieve, we retrieve the data
      if (prod(final_dims) > 0) {
        tmp <- take(ncvar_get(fnc, namevar, start, count, 
                    collapse_degen = FALSE), 
                    1:length(subset_indices), subset_indices)
        # The data is regridded if it corresponds to an atmospheric variable. When
        # the chosen output type is 'areave' the data is not regridded to not 
        # waste computing time unless the user specified a common grid.
        if (is_2d_var) {
          ###if (!is.null(work_piece[['mask']]) && !(lonlat_subsetting_requested && remap_needed)) {
          ###  mask <- take(ncvar_get(fnc_mask, work_piece[['mask']][['nc_var_name']], 
          ###               start[dim_matches[1:2]], count[dim_matches[1:2]],
          ###               collapse_degen = FALSE), 1:2, subset_indices[dim_matches[1:2]])
          ###}
          if (lonlat_subsetting_requested && remap_needed) {
            filein <- tempfile(pattern = "loadRegridded", fileext = ".nc")
            filein2 <- tempfile(pattern = "loadRegridded2", fileext = ".nc")
            ncdf_var <- ncvar_def(namevar, "", ncdf_dims[dim_matches], 
                                  fnc$var[[namevar]]$missval, 
                                  prec = if (fnc$var[[namevar]]$prec == 'int') {
                                           'integer'
                                         } else {
                                           fnc$var[[namevar]]$prec
                                         })
            nc_close(fnc)
            fnc <- nc_create(filein2, list(ncdf_var))
            ncvar_put(fnc, ncdf_var, tmp)
            nc_close(fnc)
            system(paste0("cdo -s -sellonlatbox,", if (lonmin > lonmax) {
                                                     "0,360,"
                                                   } else {
                                                     paste0(lonmin, ",", lonmax, ",")
                                                   }, latmin, ",", latmax,
                   " -", work_piece[['remap']], ",", common_grid_name, 
                   " ", filein2, " ", filein, " 2>/dev/null", sep = ""))
            file.remove(filein2)
            fnc <- nc_open(filein)
            lon <- ncvar_get(fnc, 'lon')
            lat <- ncvar_get(fnc, 'lat')
            ## We read the longitudes and latitudes from the file.
            ## In principle cdo should put in order the longitudes
            ## and slice them properly unless data is across greenwich
            lon[which(lon < 0)] <- lon[which(lon < 0)] + 360
            lon_indices <- 1:length(lon)
            if (lonmax < lonmin) {
              lon_indices <- lon_indices[which((lon <= lonmax) | (lon >= lonmin))]
            }
            lat_indices <- 1:length(lat)
            ## In principle cdo should put in order the latitudes
            if (lat[1] < lat[length(lat)]) {
              lat_indices <- length(lat):1
            }
            final_dims[c(1, 2)] <- c(length(lon_indices), length(lat_indices))
            subset_indices[[dim_matches[1]]] <- lon_indices
            subset_indices[[dim_matches[2]]] <- lat_indices

            tmp <- take(ncvar_get(fnc, namevar, collapse_degen = FALSE), 
                        1:length(subset_indices), subset_indices)

            if (!is.null(mask)) {
              ## We create a very simple 2d netcdf file that is then interpolated to the common
              ## grid to know what are the lons and lats of our slice of data
              mask_file <- tempfile(pattern = 'loadMask', fileext = '.nc')
              mask_file_remap <- tempfile(pattern = 'loadMask', fileext = '.nc')
              dim_longitudes <- ncdim_def(work_piece[['dimnames']][['lon']], "degrees_east", c(0, 360))
              dim_latitudes <- ncdim_def(work_piece[['dimnames']][['lat']], "degrees_north", c(-90, 90))
              ncdf_var <- ncvar_def('LSM', "", list(dim_longitudes, dim_latitudes), NA, 'double')
              fnc_mask <- nc_create(mask_file, list(ncdf_var))
              ncvar_put(fnc_mask, ncdf_var, array(rep(0, 4), dim = c(2, 2)))
              nc_close(fnc_mask)
              system(paste0("cdo -s ", work_piece[['remap']], ",", common_grid_name, 
                     " ", mask_file, " ", mask_file_remap, " 2>/dev/null", sep = ""))
              fnc_mask <- nc_open(mask_file_remap)
              mask_lons <- ncvar_get(fnc_mask, 'lon')
              mask_lats <- ncvar_get(fnc_mask, 'lat')
              nc_close(fnc_mask)
              file.remove(mask_file, mask_file_remap)
              if ((dim(mask)[1] != common_grid_lons) || (dim(mask)[2] != common_grid_lats)) {
                stop(paste("Error: the mask of the dataset with index ", tail(work_piece[['indices']], 1), " in '", work_piece[['dataset_type']], "' is wrong. It must be on the common grid if the selected output type is 'lonlat', 'lon' or 'lat', or 'areave' and 'grid' has been specified. It must be on the grid of the corresponding dataset if the selected output type is 'areave' and no 'grid' has been specified. For more information check ?Load and see help on parameters 'grid', 'maskmod' and 'maskobs'.", sep = ""))
              }
              mask_lons[which(mask_lons < 0)] <- mask_lons[which(mask_lons < 0)] + 360
              if (lonmax >= lonmin) {
                mask_lon_indices <- which((mask_lons >= lonmin) & (mask_lons <= lonmax))
              } else {
                mask_lon_indices <- which((mask_lons >= lonmin) | (mask_lons <= lonmax))
              }
              mask_lat_indices <- which((mask_lats >= latmin) & (mask_lats <= latmax))
              if (lat[1] < lat[length(lat)]) {
                mask_lat_indices <- mask_lat_indices[length(mask_lat_indices):1]
              }
              mask <- mask[mask_lon_indices, mask_lat_indices]
            }
            lon <- lon[lon_indices]
            lat <- lat[lat_indices]
            ###  nc_close(fnc_mask)
            ###  system(paste0("cdo -s -sellonlatbox,", if (lonmin > lonmax) {
            ###                                           "0,360,"
            ###                                         } else {
            ###                                           paste0(lonmin, ",", lonmax, ",")
            ###                                         }, latmin, ",", latmax,
            ###         " -", work_piece[['remap']], ",", common_grid_name, 
            ###This is wrong: same files  
            ###         " ", mask_file, " ", mask_file, " 2>/dev/null", sep = ""))
            ###  fnc_mask <- nc_open(mask_file) 
            ###  mask <- take(ncvar_get(fnc_mask, work_piece[['mask']][['nc_var_name']],
            ###               collapse_degen = FALSE), 1:2, subset_indices[dim_matches[1:2]])
            ###}
          }
        }
        if (!all(dim_matches == sort(dim_matches))) {
          if (!found_disordered_dims && rev(work_piece[['indices']])[2] == 1 && rev(work_piece[['indices']])[3] == 1) {
            found_disordered_dims <- TRUE
            cat(paste0("! Warning: the dimensions for the variable ", namevar, " in the files of the experiment with index ", tail(work_piece[['indices']], 1), " are not in the optimal order for loading with Load(). The optimal order would be '", paste(expected_dims, collapse = ', '), "'. One of the files of the dataset is stored in ", filename))
          }
          tmp <- aperm(tmp, dim_matches)
        }
        dim(tmp) <- final_dims
        # If we are exploring the file we don't need to process and arrange
        # the retrieved data. We only need to keep the dimension sizes.
        if (explore_dims) {
          if (work_piece[['is_file_per_member']]) {
            ## TODO: When the exp_full_path contains asterisks and is file_per_member
            ##       members from different datasets may be accounted.
            ##       Also if one file member is missing the accounting will be wrong.
            ##       Should parse the file name and extract number of members.
            if (is_url) {
              nmemb <- NULL
            } else {
              nmemb <- length(files)
            }
          }
          dims <- list(member = nmemb, time = nltime, lon = lon, lat = lat)
        } else {
        # If we are not exploring, then we have to process the retrieved data
          if (is_2d_var) {
            tmp <- apply(tmp, c(3, 4), function(x) {
              # Disable of large values.
              if (!is.na(work_piece[['var_limits']][2])) {
                x[which(x > work_piece[['var_limits']][2])] <- NA
              }
              if (!is.na(work_piece[['var_limits']][1])) {
                x[which(x < work_piece[['var_limits']][1])] <- NA
              }
              if (!is.null(mask)) {
                x[which(mask < 0.5)] <- NA
              }
  
              if (output == 'areave' || output == 'lon') {
                weights <- InsertDim(cos(lat * pi / 180), 1, length(lon))
                weights[which(is.na(x))] <- NA
                if (output == 'areave') {
                  weights <- weights / mean(weights, na.rm = TRUE)
                  mean(x * weights, na.rm = TRUE) 
                } else {
                  weights <- weights / InsertDim(Mean1Dim(weights, 2, narm = TRUE), 2, length(lat))
                  Mean1Dim(x * weights, 2, narm = TRUE)
                }
              } else if (output == 'lat') {
                Mean1Dim(x, 1, narm = TRUE)
              } else if (output == 'lonlat') {
                signif(x, 5)
              }
            })
            if (output == 'areave') {
              dim(tmp) <- c(1, 1, final_dims[3:4])
            } else if (output == 'lon') {
              dim(tmp) <- c(final_dims[1], 1, final_dims[3:4])
            } else if (output == 'lat') {
              dim(tmp) <- c(1, final_dims[c(2, 3, 4)])
            } else if (output == 'lonlat') {
              dim(tmp) <- final_dims
            }
          }
          var_data <- attach.big.matrix(work_piece[['out_pointer']])
          if (work_piece[['dims']][['member']] > 1 && nmemb > 1 && 
              work_piece[['dims']][['time']] > 1 && 
              nltime < work_piece[['dims']][['time']]) {
            work_piece[['indices']][2] <- work_piece[['indices']][2] - 1
            for (jmemb in members) {
              work_piece[['indices']][2] <- work_piece[['indices']][2] + 1
              out_position <- arrayIndex2VectorIndex(work_piece[['indices']], work_piece[['dims']])
              out_indices <- out_position:(out_position + length(tmp[, , jmemb, ]) - 1)
              var_data[out_indices] <- as.vector(tmp[, , jmemb, ])
            }
            work_piece[['indices']][2] <- work_piece[['indices']][2] - tail(members, 1) + 1
          } else {
            out_position <- arrayIndex2VectorIndex(work_piece[['indices']], work_piece[['dims']])
            out_indices <- out_position:(out_position + length(tmp) - 1)
            a <- aperm(tmp, c(1, 2, 4, 3))
            as.vector(a)
            var_data[out_indices] <- as.vector(aperm(tmp, c(1, 2, 4, 3)))
          }
          work_piece[['indices']][3] <- work_piece[['indices']][3] + 1
        }
      }
    }
    nc_close(fnc)    
    if (is_2d_var && remap_needed) {
      file.remove(filein)
      ###if (!is.null(mask) && lonlat_subsetting_requested) {
      ###  file.remove(mask_file)
      ###}
    }
    
  }
  if (explore_dims) {
    list(dims = dims, is_2d_var = is_2d_var, grid = grid_name, units = units, 
         var_long_name = var_long_name)
  } else {
    ###if (!silent && !is.null(progress_connection) && !is.null(work_piece[['progress_amount']])) {
    ###  foobar <- writeBin(work_piece[['progress_amount']], progress_connection)
    ###}
    if (!silent && !is.null(work_piece[['progress_amount']])) {
      cat(paste0(work_piece[['progress_amount']]))
    }
    found_file
  }
}

.LoadSampleData <- function(var, exp = NULL, obs = NULL, sdates, 
                            nmember = NULL, nmemberobs = NULL, 
                            nleadtime = NULL, leadtimemin = 1, 
                            leadtimemax = NULL, storefreq = 'monthly', 
                            sampleperiod = 1, lonmin = 0, lonmax = 360, 
                            latmin = -90, latmax = 90, output = 'areave', 
                            method = 'conservative', grid = NULL, 
                            maskmod = vector("list", 15), 
                            maskobs = vector("list", 15), 
                            configfile = NULL, suffixexp = NULL, 
                            suffixobs = NULL, varmin = NULL, varmax = NULL, 
                            silent = FALSE, nprocs = NULL) {
  ## This function loads and selects sample data stored in sampleMap and 
  ## sampleTimeSeries and is used in the examples instead of Load() so as
  ## to avoid nco and cdo system calls and computation time in the stage 
  ## of running examples in the CHECK process on CRAN.
  selected_start_dates <- match(sdates, c('19851101', '19901101', '19951101', 
                                          '20001101', '20051101'))
  start_dates_position <- 3
  lead_times_position <- 4

  if (output == 'lonlat') {
    sampleData <- s2dverification::sampleMap
    if (is.null(leadtimemax)) {
      leadtimemax <- dim(sampleData$mod)[lead_times_position]
    }
    selected_lead_times <- leadtimemin:leadtimemax

    dataOut <- sampleData
    dataOut$mod <- sampleData$mod[, , selected_start_dates, selected_lead_times, , ]
    dataOut$obs <- sampleData$obs[, , selected_start_dates, selected_lead_times, , ]
  }
  else if (output == 'areave') {
    sampleData <- s2dverification::sampleTimeSeries
    if (is.null(leadtimemax)) {
      leadtimemax <- dim(sampleData$mod)[lead_times_position]
    }
    selected_lead_times <- leadtimemin:leadtimemax

    dataOut <- sampleData
    dataOut$mod <- sampleData$mod[, , selected_start_dates, selected_lead_times]
    dataOut$obs <- sampleData$obs[, , selected_start_dates, selected_lead_times]
  }

  dims_out <- dim(sampleData$mod)
  dims_out[start_dates_position] <- length(selected_start_dates)
  dims_out[lead_times_position] <- length(selected_lead_times)
  dim(dataOut$mod) <- dims_out

  dims_out <- dim(sampleData$obs)
  dims_out[start_dates_position] <- length(selected_start_dates)
  dims_out[lead_times_position] <- length(selected_lead_times)
  dim(dataOut$obs) <- dims_out

  invisible(list(mod = dataOut$mod, obs = dataOut$obs, 
                 lat = dataOut$lat, lon = dataOut$lon))
}

.ConfigGetDatasetInfo <- function(matching_entries, table_name) {
  # This function obtains the information of a dataset and variable pair,
  # applying all the entries that match in the configuration file.
  if (table_name == 'experiments') {
    id <- 'EXP'
  } else {
    id <- 'OBS'
  }
  defaults <- c(paste0('$DEFAULT_', id, '_MAIN_PATH$'), paste0('$DEFAULT_', id, '_FILE_PATH$'), '$DEFAULT_NC_VAR_NAME$', '$DEFAULT_SUFFIX$', '$DEFAULT_VAR_MIN$', '$DEFAULT_VAR_MAX$')
  info <- NULL

  for (entry in matching_entries) {
    if (is.null(info)) {
      info <- entry[-1:-2]
      info[which(info == '*')] <- defaults[which(info == '*')]
    } else {
      info[which(entry[-1:-2] != '*')] <- entry[-1:-2][which(entry[-1:-2] != '*')]
    }
  }

  info <- as.list(info)
  names(info) <- c('main_path', 'file_path', 'nc_var_name', 'suffix', 'var_min', 'var_max')  
  info
}
