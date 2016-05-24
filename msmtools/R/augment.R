#' A fast and general method for building augmented data.
#'
#' A fast and general method for reshaping standard longitudinal data, where each observation contains
#' the exact starting and ending time of the process, into a new structure called 'augmented'
#' which is suitable for multi-state analyses using the \code{\link[msm]{msm}} package.
#'
#' @param data A \code{data.table} object where each row represents an observation.
#' @param data_key A keying variable which \code{augment} uses to define a key for \code{data}.
#' This represents the subject ID.
#' @param n_events An integer variable indicating the progressive (monotonic) event number
#' of a given ID. If missing, \code{augment} fastly creates a variable named \code{"n_events"}.
#' \code{augment} always checks whether \code{n_events} is monotonic increasing within the
#' provided \code{data_key} and stops the execution in case the check fails (see 'Details').
#' @param pattern Either an integer, a factor or a characer with 2 or 3 unique values which
#' provides the ID status at the end of the study. \code{pattern} has a predefined structure.
#' When 2 values are detected, they must be in the format: 0 = "alive", 1 = "dead". When 3 values
#' are detected, then the format must be: 0 = "alive", 1 = "dead during a transition",
#' 2 = "dead after a transition has ended".
#' @param state A list of 3 possible states which a subject can reach. \code{state} has a predefined
#' structure as follows: IN, OUT, DEAD. If \code{augment} does not detect exactly 3 states, it
#' stops with error.
#' @param t_start The starting time of an observation. It can be passed as date, integer, or numeric
#' format.
#' @param t_end The ending time of an observation. It can be passed as date, integer, or numeric
#' format.
#' @param t_cens The censoring time of the study
#' @param t_death The exact death time of a subject ID. If \code{t_death} is missing,
#' \code{t_cens} is assumed to contain both censoring and death time.
#' @param t_augmented The new time variable of the process in the augmented format.
#' If \code{t_augmented} is missing, then the default name 'augmented' is assumed.
#' The variable is added to \code{data}. If \code{augment} detects a date or a difftime
#' format in \code{t_start}, then \code{t_augmented} is cast to integer or to numeric, respectively,
#' into a new variable with the suffix '_int' or '_num' added to \code{t_augmented}
#' This is done because \code{\link[msm]{msm}} can't deal correclty if time variable is a date
#' or a difftime. Both variables are positioned before \code{t_start}.
#' @param more_status A variable which marks further transitions beside the default given by
#' \code{state}. \code{more_status} can be a factor or a character. In particular, if the
#' corresponding observation is a standard admission (i.e. no other information available), then
#' \code{more_status} must be set to 'df' which stands for 'Default' (see 'Examples' and ?hosp).
#' If missing, \code{augment} ignores it.
#' @param verbose If \code{FALSE}, all information produced by \code{print}, \code{cat} and
#' \code{message} are suppressed. All is done internally so that no global
#' options are changed. \code{verbose} can be set to \code{FALSE} on all common OS
#' (see also \code{\link[base]{sink}} and \code{\link[base]{options}}). Default is \code{TRUE}.
#' @details In order to get the data processed, a monotonic increasing process needs to be ensured.
#' \code{augment} checks this both in case \code{n_events} is missing or not. The data are
#' fastly ordered through \code{\link[data.table]{setkey}} function and using \code{data_key} as
#' the primary key and \code{t_start} as secondary key. Then it checks \code{n_events} and if it fails,
#' it returns the subjects gived by \code{data_key} where issues occurred before giving an
#' error and stopping. If \code{n_events} is not passed, then the ordering procedure remains the
#' same, but the progression number is created internally with the name \code{n_events}.
#' @return A restructured long format dataset of class \code{"data.table"} where each row
#' represents a specific transition for a given subject. Moreover, \code{augment} adds some
#' important variables:\cr
#' -----\cr
#' \emph{augmented}: the new timing variable for the process when looking at transitions. If
#' \code{t_augmented} is missing, then \code{augment} creates \emph{augmented} by default.
#' \emph{augmented}. The function looks directly to \code{t_start} and \code{t_end} to build
#' it and thus it inherits their class.
#' In particular, if \code{t_start} is a date format, then \code{augment} computes a new variable
#' cast as integer and names it \emph{augmented_int}. If \code{t_start} is a difftime format,
#' then \code{augment} computes a new variable cast as a numeric and names it \emph{augmented_num};\cr
#' \emph{status}: a status flag which looks at \code{state}. \code{augment} automatically checks
#' whether argument \code{pattern} has 2 or 3 unique values and computes the correct structure of
#' a given subject. The variable is cast as character;\cr
#' \emph{status_num}: the corresponding integer version of \emph{status};\cr
#' \emph{n_status}: a mix of \emph{status} and \code{n_events} cast as character.
#' \emph{n_status} comes into play when a model on the progression of the process is intended.\cr
#' -----\cr
#'
#' If \code{more_status} is passed, then \code{augment} computes some more variables. They mimic the
#' meaning of \emph{status}, \emph{status_num}, and \emph{n_status} but they account for the more
#' complex structure defined. They are: \code{status_exp}, \code{status_exp_num}, and
#' \code{n_status_exp}.
#' @examples
#' # 1.
#' # loading data
#' data( hosp )
#'
#' # augmenting hosp
#' hosp_augmented = augment( data = hosp, data_key = subj, n_events = adm_number, pattern = label_3,
#'                           t_start = dateIN, t_end = dateOUT, t_cens = dateCENS )
#'
#' # 2.
#' # augmenting hosp by passing more information regarding transition with arg. more_status
#' hosp_augmented_more = augment( data = hosp, data_key = subj, n_events = adm_number,
#'                                pattern = label_3, t_start = dateIN, t_end = dateOUT,
#'                                t_cens = dateCENS, more_status = rehab_it )
#'
#' \dontrun{
#' augmented = augment( data = hosp, data_key = subj, n_events = dateIN,
#'                      pattern = label_3, t_start = dateIN, t_end = dateOUT, t_cens = dateCENS ) }
#' @references Jackson, C.H. (2011). Multi-State Models for Panel Data:
#' The \emph{msm} Package for R. Journal of Statistical Software, 38(8), 1-29.
#' URL \url{http://www.jstatsoft.org/v38/i08/}.
#' @author Francesco Grossetti \email{francesco.grossetti@@polimi.it}.
#' @export
augment = function( data, data_key, n_events, pattern, state = list ( 'IN', 'OUT', 'DEAD' ),
                    t_start, t_end, t_cens, t_death, t_augmented = 'augmented',
                    more_status, verbose = TRUE ) {

  tic = proc.time()
  status         = NULL
  status_num     = NULL
  n_status       = NULL
  status_exp     = NULL
  status_exp_num = NULL
  n_status_exp   = NULL
  .              = NULL
  oldw = getOption( "warn" )
  if ( verbose == TRUE ) {
    options( warn = 1 )
  }
  if ( verbose == FALSE ) {
    options( warn = 1 )
    if ( .Platform$OS.type == 'windows' ) {
      sink( file = "NUL" )
    } else {
      sink( file = "/dev/null" )
    }
    cat( '---\n' )
  }
  if ( missing( data ) ) {
    stop( 'a dataset of class data.table must be provided' )
  }
  if ( !inherits( data, "data.table" ) ) {
    stop( "a dataset of class data.table must be provided" )
  }
  if ( missing( data_key ) ) {
    stop( 'a variable of keying must be provided' )
  }
  if ( missing( pattern ) ) {
    stop( "a pattern must be provided" )
  }
  if ( !inherits( state, "list" ) || length( state ) != 3 ) {
    stop( "state pattern must be a list of 3 elements" )
  }
  if ( missing( t_start ) || missing( t_end ) ) {
    stop( 'a starting and an ending event times must be provided' )
  }
  if ( missing( t_cens ) ) {
    stop( 'a censoring time must be provided' )
  }
  if ( missing( t_death ) ) {
    warning( 'no t_death has been passed. Assuming that ', substitute( t_cens ),
             ' contains both censoring and death time' )
  }
  if ( verbose == TRUE ) {
    cat( '---\n' )
    cat( '# # # # setting everything up # # # #\n' )
    cat( '---\n' )
  }
  pattern = as.character( substitute( list( pattern ) )[ -1L ] )
  t_start = as.character( substitute( list( t_start ) )[ -1L ] )
  t_end = as.character( substitute( list( t_end ) )[ -1L ] )
  t_cens = as.character( substitute( list( t_cens ) )[ -1L ] )

  if ( eval( substitute( class( data$t_start ) ) ) != eval( substitute( class( data$t_end ) ) ) ) {
    stop( 'the starting and the ending event times must be of the same class' )
  } else if ( eval( substitute( class( data$t_start ) ) ) != eval( substitute( class( data$t_cens ) ) ) ) {
    stop( 'the starting and the censoring event times must be of the same class' )
  }
  setkey( data, NULL )
  if ( !missing( n_events ) ) {
    cols = as.character( substitute( list( data_key, n_events ) )[ -1L ] )
    if ( !length( cols ) )
      cols = colnames( data )
    if ( !inherits( eval( substitute( data$n_events ) ), "integer" ) ) {
      stop( 'n_events must be an integer' )
    }
    if ( verbose == TRUE ) {
      message( 'checking monotonicity of ', cols[[ 2 ]] )
    }
    ev = data[ , .( ev = all( get( cols[[ 2 ]] ) == cummax( get( cols[[ 2 ]] ) ) ) ),
               by = eval( cols[[ 1 ]] ) ]
    setkeyv( data, c( cols[[ 1 ]], t_start ) )
    if ( all( ev$ev ) == FALSE ) {
      message( substitute( n_events ), ' is not monotonic increasing within ',
               substitute( data_key ) )
      message( 'the corresponding subjects are:' )
      message( paste( ev[ ev == FALSE ][ , get( cols[[ 1 ]] ) ], collapse = '; ' ) )
      stop( 'Please, fix the issues and relaunch augment()' )
    } else {
      cat( 'Ok, ', cols[[ 2 ]], ' is monotonic\n', sep = '' )
      cat( '---\n' )
    }
    setkeyv( data, cols )
  } else {
    cols = as.character( substitute( list( data_key ) )[ -1L ] )
    if ( !length( cols ) )
      cols = colnames( data )
    setkeyv( data, c( cols, t_start ) )
    data[ , n_events := seq( .N ), by = eval( cols ) ]
    cols = c( cols, names( data )[ dim( data )[ 2 ] ] )
    if ( verbose == TRUE ) {
      message( 'checking monotonicity of ', cols[[ 2 ]] )
    }
    ev = data[ , .( ev = all( get( cols[[ 2 ]] ) == cummax( get( cols[[ 2 ]] ) ) ) ),
               by = eval( cols[[ 1 ]] ) ]
    setkeyv( data, c( cols[[ 1 ]], t_start ) )
    if ( all( ev$ev ) == FALSE ) {
      message( substitute( n_events ), ' is not monotonic increasing within ',
               substitute( data_key ) )
      message( 'the corresponding subjects are:' )
      message( paste( ev[ ev == FALSE ][ , get( cols[[ 1 ]] ) ], collapse = '; ' ) )
      stop( 'Please, fix the issues and relaunch augment()' )
    } else {
      cat( 'Ok, n_events is monotonic\n' )
      cat( '---\n' )
    }
    setkeyv( data, cols )
  }
  checks = c( cols, pattern, t_start, t_end )
  if ( !missing( t_death ) ) {
    t_death = as.character( substitute( list( t_death ) )[ -1L ] )
    if ( eval( substitute( class( data$t_cens ) ) ) != eval( substitute( class( data$t_death ) ) ) ) {
      stop( 'the censoring and the death event times must be of the same class' )
    }
  }
  if ( verbose == TRUE ) {
    message( 'checking for any missing data in function arguments' )
  }
  test = apply( data[ , checks, with = FALSE ], 2, function( x ) any( sum( is.na( x ) ) > 0 ) )
  if ( any ( test ) ) {
    cat( '---\n' )
    if ( verbose == TRUE ) {
      message( 'detected missing values in the following variables:' )
    }
    invisible( sapply( names( test[ test == TRUE ] ), function( x ) cat( x, '\n' ) ) )
    stop( 'Please, fix the issues and relaunch augment()' )
  } else {
    cat( 'Ok, no missing data detected\n')
    cat( '---\n' )
  }

  l = lapply( 1:2, function( x ) data )
  expand = rbindlist( l )
  setkeyv( expand, cols )
  values = sort( eval( substitute( unique( data$pattern ) ) ) )

  if ( verbose == TRUE ) {
    message( 'checking ', substitute( pattern ), ' and defining patterns' )
  }
  if ( length( values ) < 2 ) {
    stop( 'unit identification label must be an integer, a factor or a character
          with at least 2 elements' )
  } else if ( length( values ) == 2 ) {
    cat( 'detected only 2 values\n' )
    cat( '---\n' )
    if ( inherits( eval( substitute( unique( data$pattern ) ) ), 'integer' ) ||
         inherits( eval( substitute( unique( data$pattern ) ) ), 'numeric' ) ) {
      match1 = data[ data[ get( pattern ) == 0, .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
      if ( missing( t_death ) ) {
        match3 = data[ data[ get( pattern ) == 1,
                             .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ][ get( t_end ) != get( t_cens ) ]
      } else {
        match3 = data[ data[ get( pattern ) == 1,
                             .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ][ get( t_end ) != get( t_death ) ]
      }
    } else if ( inherits( eval( substitute( unique( data$pattern ) ) ), 'factor' ) ) {
      match1 = data[ data[ as.integer( get( pattern ) ) - 1 == 0,
                           .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
      if ( missing( t_death ) ) {
        match3 = data[ data[ as.integer( get( pattern ) ) - 1 == 1,
                             .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ][ get( t_end ) != get( t_cens ) ]
      } else {
        match3 = data[ data[ as.integer( get( pattern ) ) - 1 == 1,
                             .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ][ get( t_end ) != get( t_death ) ]
      }
    } else if ( inherits( eval( substitute( unique( data$pattern ) ) ), 'character' ) ) {
      match1 = data[ data[ get( pattern ) == values[ 1 ], .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
      if ( missing( t_death ) ) {
        match3 = data[ data[ get( pattern ) == values[ 2 ],
                             .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ][ get( t_end ) != get( t_cens ) ]
      } else {
        match3 = data[ data[ get( pattern ) == values[ 2 ],
                             .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ][ get( t_end ) != get( t_death ) ]
      }
    }
  } else if ( length( values ) == 3 ) {
    cat( 'Ok, detected 3 values\n' )
    cat( '---\n' )
    if ( inherits( eval( substitute( unique( data$pattern ) ) ), 'integer' ) ||
         inherits( eval( substitute( unique( data$pattern ) ) ), 'numeric' ) ) {
      match1 = data[ data[ get( pattern ) == 0, .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
      match3 = data[ data[ get( pattern ) == 2, .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
    } else if ( inherits( eval( substitute( unique( data$pattern ) ) ), 'factor' ) ) {
      match1 = data[ data[ as.integer( get( pattern ) ) - 1 == 0,
                           .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
      match3 = data[ data[ as.integer( get( pattern ) ) - 1 == 2,
                           .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
    } else if ( inherits( eval( substitute( unique( data$pattern ) ) ), 'character' ) ) {
      match1 = data[ data[ get( pattern ) == values[ 1 ], .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
      match3 = data[ data[ get( pattern ) == values[ 3 ], .I[ .N ], by = eval( cols[[ 1 ]] ) ]$V1 ]
    }
  }
  if ( verbose == TRUE ) {
    message( 'augmenting data' )
  }
  l = list( expand, match1, match3 )
  final = rbindlist( l )
  setkeyv( final, cols )
  cat( 'data have been augmented!\n' )
  cat( '---\n' )

  if ( length( values ) == 2 ) {
    if ( missing( t_death ) ) {
      counter_events = final[ , .( N = .N,
                                   j = unique( get( pattern ) ),
                                   t_end = max( get( t_end ) ),
                                   t_cens = max( get( t_cens ) ) ), by = eval( cols[[ 1 ]] ) ]
      s = dim( counter_events )[ 1 ]
      status_flag_temp = vector( mode = 'list', dim( counter_events )[ 1 ] )
    } else {
      counter_events = final[ , .( N = .N,
                                   j = unique( get( pattern ) ),
                                   t_end = max( get( t_end ) ),
                                   t_death = max( get( t_death ) ) ), by = eval( cols[[ 1 ]] ) ]
      s = dim( counter_events )[ 1 ]
      status_flag_temp = vector( mode = 'list', dim( counter_events )[ 1 ] )
    }
  } else if ( length( values ) == 3 ) {
    counter_events = data[ , .( N = .N,
                                j = unique( get( pattern ) ) ), by = eval( cols[[ 1 ]] ) ]
    s = dim( counter_events )[ 1 ]
    status_flag_temp = vector( mode = 'list', dim( counter_events )[ 1 ] )
  }
  if ( verbose == TRUE ) {
    message( 'adding status flag' )
  }
  for ( i in seq_along( counter_events$N ) ) {
    if ( s > 10000 ) {
      if ( i %% 5000 == 0 ) {
        cat( '* * * iteration', i, 'of', s, '\n' )
      }
    } else {
      if ( i %% 1000 == 0 ) {
        cat( '* * * iteration', i, 'of', s, '\n' )
      }
    }
    if ( length( values ) == 2 ) {
      if ( missing( t_death ) ) {
        if ( counter_events$j[ i ] == values[ 1 ] ) {
          status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ),
                                            ( floor( counter_events$N[ i ] / 2 ) ) ),
                                       state [[ 2 ]] )
        } else if ( counter_events$j[ i ] == values[ 2 ] &
                    counter_events$t_end[ i ] == counter_events$t_cens[ i ] ) {
          status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ),
                                            ( floor( counter_events$N[ i ] / 2 ) - 1 ) ),
                                       state[[ 1 ]], state[[ 3 ]] )
        } else if ( counter_events$j[ i ] == values[ 2 ] &
                    counter_events$t_end[ i ] != counter_events$t_cens[ i ] ) {
          status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ),
                                            ( floor( counter_events$N[ i ] / 2 ) ) ), state[[ 3 ]] )
        }
      } else {
        if ( counter_events$j[ i ] == values[ 1 ] ) {
          status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ),
                                            ( floor( counter_events$N[ i ] / 2 ) ) ),
                                       state [[ 2 ]] )
        } else if ( counter_events$j[ i ] == values[ 2 ] &
                    counter_events$t_end[ i ] == counter_events$t_death[ i ] ) {
          status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ),
                                            ( floor( counter_events$N[ i ] / 2 ) - 1 ) ),
                                       state[[ 1 ]], state[[ 3 ]] )
        } else if ( counter_events$j[ i ] == values[ 2 ] &
                    counter_events$t_end[ i ] != counter_events$t_death[ i ] ) {
          status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ),
                                            ( floor( counter_events$N[ i ] / 2 ) ) ), state[[ 3 ]] )
        }
      }
    } else if ( length( values ) == 3 ) {
      if ( counter_events$j[ i ] == values[ 1 ] ) {
        status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ),
                                          counter_events$N[ i ] ), state [[ 2 ]] )
      } else if ( counter_events$j[ i ] == values[ 2 ] ) {
        status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ),
                                          ( counter_events$N[ i ] - 1 ) ), state[[ 1 ]],
                                     state[[ 3 ]] )
      } else if ( counter_events$j[ i ] == values[ 3 ] ) {
        status_flag_temp[[ i ]] = c( rep( c( state[[ 1 ]], state[[ 2 ]] ), counter_events$N[ i ] ),
                                     state[[ 3 ]] )
      }
    }
  }
  status_flag = unlist( status_flag_temp, recursive = FALSE )
  final[ , status := status_flag ]
  if ( sum( is.na( final$status ) ) == 0 ) {
    cat( 'status flag has been added successfully \n' )
    cat( '---\n' )
  } else {
    stop( 'status flag has not been build correctly' )
  }
  if ( verbose == TRUE ) {
    message( 'adding numeric status flag' )
  }
  k = uniqueN( final$status )
  lev = unique( final$status )
  for ( i in 1:k ) {
    final[ status == lev[ i ], status_num := i ]
  }
  if ( i == k ) {
    cat( 'numeric status has been added successfully \n' )
    cat( '---\n' )
  } else {
    stop( 'numeric status status has not been build correctly' )
  }
  if ( verbose == TRUE ) {
    message( 'adding sequential status flag' )
  }
  if ( missing( n_events ) ) {
    final[ , n_status := ifelse( status != state[[ 3 ]],
                                 paste( n_events, ' ', status, sep = '' ),
                                 state[[ 3 ]] ) ]
  } else {
    final[ , n_status := ifelse( status != state[[ 3 ]],
                                 paste( eval( substitute( n_events ) ), ' ', status, sep = '' ),
                                 state[[ 3 ]] ) ]
  }
  if ( sum( is.na( final$n_status ) ) == 0 ) {
    cat( 'sequential status flag has been added successfully \n' )
    cat( '---\n' )
  } else {
    stop( 'sequential status flag has not been build correctly' )
  }
  if ( verbose == TRUE ) {
    message( 'adding variable ', substitute( t_augmented ), ' as new time variable' )
  }
  final[ status == state[[ 1 ]], substitute( t_augmented ) := get( t_start ) ]
  final[ status == state[[ 2 ]], substitute( t_augmented ) := get( t_end ) ]
  if ( missing( t_death ) ) {
    final[ status == state[[ 3 ]], substitute( t_augmented ) := get( t_cens ) ]
  } else {
    final[ status == state[[ 3 ]], substitute( t_augmented ) := get( t_death ) ]
  }
  if ( inherits( eval( substitute( data$t_start ) ), 'Date' ) ) {
    final[ , paste( substitute( t_augmented ), '_int', sep = '' ) := as.integer( get( t_augmented ) ) ]
    id_col = which( names( data ) == substitute( t_start ) )
    setcolorder( final, c( 1:( id_col - 1 ), ( dim( final )[ 2 ] - 1 ), dim( final )[ 2 ],
                           id_col:( dim( final )[ 2 ] - 2 ) ) )
    cat( 'variables \"', substitute( t_augmented ), '\" and \"',
         paste( substitute( t_augmented ), '_int', sep = '' ),
         '\" successfully added and repositioned\n', sep = '' )
    cat( '---\n' )
  } else if ( inherits( eval( substitute( data$t_start ) ), 'difftime' ) ) {
    final[ , paste( substitute( t_augmented ), '_num', sep = '' ) := as.numeric( get( t_augmented ) ) ]
    id_col = which( names( data ) == substitute( t_start ) )
    setcolorder( final, c( 1:( id_col - 1 ), ( dim( final )[ 2 ] - 1 ), dim( final )[ 2 ],
                           id_col:( dim( final )[ 2 ] - 2 ) ) )
    cat( 'variables \"', substitute( t_augmented ), '\" and \"',
         paste( substitute( t_augmented ), '_num', sep = '' ),
         '\" successfully added and repositioned\n', sep = '' )
    cat( '---\n' )
  } else if ( inherits( eval( substitute( data$t_start ) ), 'integer' ) ||
              inherits( eval( substitute( data$t_start ) ), 'numeric' ) ) {

    id_col = which( names( data ) == substitute( t_start ) )
    setcolorder( final, c( 1:( id_col - 1 ), dim( final )[ 2 ],
                           id_col:( dim( final )[ 2 ] - 1 ) ) )
    cat( 'variable \"', substitute( t_augmented ),
         '\" successfully added and repositioned\n', sep = '' )
    cat( '---\n' )
  }

  if ( !missing( more_status ) ) {
    more_status = as.character( substitute( list( more_status )  )[ -1L ] )
    test = apply( data[ , more_status, with = FALSE ], 2,
                  function( x ) any( sum( is.na( x ) ) > 0 ) )
    if ( any ( test ) ) {
      cat( '---\n' )
      message( 'detected missing values in variable:' )
      invisible( sapply( names( test[ test == TRUE ] ), function( x ) cat( x, '\n' ) ) )
      stop( 'Please, fix the issues and relaunch augment()' )
    }

    values = eval( substitute( unique( data$more_status ) ) )
    if ( verbose == TRUE ) {
      message( 'adding expanded status flag' )
    }
    final[ status == state[[ 3 ]], status_exp := state[[ 3 ]] ]
    for ( i in seq_along( values ) ) {
      final[ status != state[[ 3 ]] & get( more_status ) == values[ i ],
             status_exp := paste( values[ i ], '_', status, sep = '' ) ]
    }
    if ( sum( is.na( final$status_exp ) ) == 0 ) {
      cat( 'expanded status flag has been added successfully \n' )
      cat( '---\n' )
    } else {
      stop( 'expanded status flag has not been build correctly' )
    }
    if ( verbose == TRUE ) {
      message( 'adding numeric expanded status flag' )
    }
    k = uniqueN( final$status_exp )
    lev = unique( final$status_exp )
    for ( i in 1:k ) {
      final[ status_exp == lev[ i ], status_exp_num := i ]
    }
    if ( i == k ) {
      cat( 'expanded numeric status has been added successfully \n' )
      cat( '---\n' )
    } else {
      stop( 'expanded numeric status has not been build correctly' )
    }
    if ( verbose == TRUE ) {
      message( 'adding sequential expanded status flag' )
    }
    final[ , n_status_exp := ifelse( status_exp != state[[ 3 ]],
                                     paste( eval( substitute( n_events ) ), ' ',
                                            status_exp, sep = '' ), state[[ 3 ]] ) ]
    if ( sum( is.na( final$n_status_exp ) ) == 0 ) {
      cat( 'expanded sequential status flag has been added successfully \n' )
    } else {
      stop( 'expanded sequential status flag has not been build correctly' )
    }
  }

  toc = proc.time()
  time = toc - tic
  cat( '---------------------------\n' )
  cat( 'augment() took:', time[ 3 ], 'sec. \n', sep = ' ' )
  cat( '---------------------------\n' )
  if ( verbose == FALSE ) {
    sink()
  }
  options( warn = oldw )
  return( final )
}

