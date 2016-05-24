#' Common string transformations 
#' 
#' Perform common transformations on strings
#' 
#' @param string character vector
#' 
#' These function take a single character vector argument and return a character
#' vector that has had one or more functions applied to it. They are the 
#' building blocks for building up case transformations.
#' 
#' @return character vector
#' 
#' @examples 
#' 
#' # TRANSFORMATIONS
#'   str_ucfirst( "abc def" )              # Abc Def
#'   str_expand_capwords( "AbcDef")        # Abc Def
#'   
#' # DELETION
#'   str_delete_whitespace( "ABC 123" )    # ABC123
#'   str_delete_separators( "A_B-C.123" )  # ABC123
#'   str_delete_nonword( "ABC & 123" )     # ABC123
#'   
#' @rdname string_transformations
#' @include make_str_replace.R patterns.R
#' @export



#' @rdname string_transformations
#' @export
  str_delete_whitespace <- make_str_delete( pattern=pattern_whitespace )


#' @rdname string_transformations
#' @export
  str_delete_separators <- make_str_delete( pattern=pattern_separators )


#' @rdname string_transformations
#' @export
  str_delete_nonword    <- make_str_delete( pattern=pattern_nonword )


#' @rdname string_transformations
#' @export
  str_expand_capwords <- 
    make_str_replace( pattern = '([a-z])([A-Z])', replacement = '\\1 \\2' )


#' @rdname string_transformations
#' @export
  str_delete_leading_nonword <- 
    make_str_delete( pattern = '^\\W' )


#' @rdname string_transformations
#' @export 
  str_delete_space <- make_str_delete( pattern = fixed(' ') ) 
