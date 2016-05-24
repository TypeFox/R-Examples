#' @include imports.R

#' @rdname is_uk_car_licence
#' @export
assert_all_are_uk_car_licences <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all UK car licence plates.", get_name_in_parent(x))
  assert_engine(
    is_uk_car_licence, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )        
}

#' @rdname is_uk_car_licence
#' @export
assert_any_are_uk_car_licences <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not UK car licence plates.", get_name_in_parent(x))
  assert_engine(
    is_uk_car_licence, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )        
}

#' @rdname is_uk_car_licence
#' @export
assert_all_are_uk_car_licenses <- assert_all_are_uk_car_licences

#' @rdname is_uk_car_licence
#' @export
assert_any_are_uk_car_licenses <- assert_any_are_uk_car_licences

#' @rdname is_uk_national_insurance_number
#' @export
assert_all_are_uk_national_insurance_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf(
    "%s are not all UK national insurance numbers.", 
    get_name_in_parent(x)
  )
  assert_engine(
    is_uk_national_insurance_number, 
    x, 
    msg = msg,  
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_uk_national_insurance_number
#' @export
assert_any_are_uk_national_insurance_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf(
    "%s are all not UK national insurance numbers.", 
    get_name_in_parent(x)
  )
  assert_engine(
    is_uk_national_insurance_number, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )        
}

#' @rdname is_uk_postcode
#' @export
assert_all_are_uk_postcodes <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all UK postcodes.", get_name_in_parent(x))
  assert_engine(
    is_uk_postcode, 
    x, 
    msg = msg, 
    na_ignore = na_ignore,
    severity = severity
  )          
}

#' @rdname is_uk_postcode
#' @export
assert_any_are_uk_postcodes <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not UK postcodes.", get_name_in_parent(x))
  assert_engine(
    is_uk_postcode, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )     
}

#' @rdname is_uk_telephone_number
#' @export
assert_all_are_uk_telephone_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are not all UK telephone numbers.", get_name_in_parent(x))
  assert_engine(
    is_uk_telephone_number, 
    x, 
    msg = msg,  
    na_ignore = na_ignore,
    severity = severity
  )  
}

#' @rdname is_uk_telephone_number
#' @export
assert_any_are_uk_telephone_numbers <- function(x, na_ignore = FALSE, 
  severity = getOption("assertive.severity", "stop"))
{                                                     
  msg <- gettextf("%s are all not UK telephone numbers.", get_name_in_parent(x))
  assert_engine(
    is_uk_telephone_number, 
    x, 
    msg = msg, 
    what = "any", 
    na_ignore = na_ignore,
    severity = severity
  )  
}
