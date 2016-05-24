#' @name hashid_defaults
#'
#' @title Default Values for hashid settings
#' @description Default alphabet, separators, and 
#'    ratio of character separators and guards for hashid
#'   
#' @source 
#' http://www.hashids.org
#' 

#' @rdname hashid_defaults
#' @export
DEFAULT_ALPHABET =  "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ1234567890"

#' @rdname hashid_defaults
#' @export
DEFAULT_SEPS = "cfhistuCFHISTU"

#' @rdname hashid_defaults
#' @export
RATIO_SEPARATORS = 3.5

#' @rdname hashid_defaults
#' @export
RATIO_GUARDS = 12