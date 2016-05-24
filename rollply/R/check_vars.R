# 
# 
# 
# Checks variables in input data.frame and/or grid.
# 

check_coords <- function(coords, grid) { 
  
  # Check if all coordinates involved are numeric
  all_numeric <- all(sapply(coords, is.numeric))
  if (!is.null(grid)) { 
    all_numeric <- all_numeric & all(sapply(coords, is.numeric))
  }
  
  if ( ! all_numeric ) { 
    stop("Can not work with non-numeric coordinates.")
  }
  
  # Check for the presence of NA's
  na_in_coords <- apply(coords, 2, function(X) any(is.na(X)))
  if (! is.null(grid) ) {
    na_in_grid   <- apply(grid, 2, function(X) any(is.na(X)))
  } else { 
    na_in_grid <- FALSE
  }
  
  # We do not remove NAs as this would create a copy of 
  # a potentially big dataset.
  if (na_in_coords || na_in_grid) { 
    stop(paste0('NA in coordinates are not supported. ',
                'Try removing them before calling rollply.'))
  }
  
  if (!is.null(grid)) { 
    if ( ! all(colnames(coords) %in% colnames(grid)) ) { 
      stop(paste0("Coordinates columns not found in grid, please check ",
                  "the provided grid."))
    } 
  }
  
}

check_args <- function(.rollvars, 
                       .data, 
                       grid, 
                       grid.type, 
                       wdw.size) { 
  
  # Checks for proper grid type
  if ( ! exists(paste0('build_grid_', grid.type)) ) {
    stop('Unknown grid type. See ?build_grid for a list of supported grid types')
  }
  
  if (  ! is_mat_or_df(.data) ) {
    stop("I do not know what to do with argument .data of class ", class(.data))
  }
  
  # If a grid object was passed
  if ( ! is.null(grid) ) { 
    
    # Rollply accepts matrices as grid, too
    if ( ! is_mat_or_df(grid) ) { 
      stop("I do not know what to do with argument grid of class ", class(grid))
    }
    
    # Need colnames !
    if ( is.null(colnames(grid)) ) { 
      stop("Grid object should have column names")
    }
    
  }
  
  if ( ! is.numeric(wdw.size) || length(wdw.size) > 1 ) { 
    stop('Invalid or missing window size (argument wdw.size), ', 
         'please check arguments')
  }
}

is_mat_or_df <- function(obj) { 
  inherits(obj,'matrix') | inherits(obj, 'data.frame')
}
