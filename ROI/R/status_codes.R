## ROI: status_codes.R
## overview of solver status codes and their canonicalization

## adds a new status code to db, default roi_code is 1L, i.e. a failure
add_status_code_to_db <- function(solver, code, symbol, message, roi_code = 1L){
  status_db$set_entry(solver = solver,
                      code = code,
                      symbol = symbol,
                      message = message,
                      roi_code = roi_code)
}

get_status_message_from_db <- function(solver, code){
  status_db[[solver, code]]
}

delete_status_code_from_db <- function(solver, code){
  status_db$delete_entry(solver = solver,
                         code = code)
}

available_in_status_codes_db <- function( )
  unique( status_db$get_field_entries("solver") )

