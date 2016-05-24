# From assertive.types is-type-base.R, assert-is-type-base.R, is-type.R, assert-is-a-type.R

#' Does x belong to these classes?
#' 
#' See \code{\link[assertive.types]{assert_is_all_of}}.
#' @name assert_is_all_of
#' @aliases assert_is_any_of
#' @export assert_is_all_of
#' @export assert_is_any_of
NULL

#' Is the input an array or matrix?
#' 
#' See \code{\link[assertive.types]{is_array}}.
#' @name is_array
#' @aliases is_matrix assert_is_array assert_is_matrix
#' @export is_array
#' @export is_matrix
#' @export assert_is_array
#' @export assert_is_matrix
NULL

#' Is the input of type character?
#' 
#' See \code{\link[assertive.types]{is_character}}.
#' @name is_character
#' @aliases is_a_string assert_is_character assert_is_a_string
#' @export is_character
#' @export is_a_string
#' @export assert_is_character
#' @export assert_is_a_string
NULL

#' Is the input complex?
#' 
#' See \code{\link[assertive.types]{is_complex}}.
#' @name is_complex
#' @aliases is_complex is_a_complex assert_is_complex assert_is_a_complex
#' @export is_complex
#' @export is_a_complex
#' @export assert_is_complex
#' @export assert_is_a_complex
NULL

#' Is the input is a data.frame?
#' 
#' See \code{\link[assertive.types]{is_data.frame}}.
#' @name is_data.frame
#' @aliases assert_is_data.frame
#' @export is_data.frame
#' @export assert_is_data.frame
NULL

#' Is the input an environment?
#' 
#' See \code{\link[assertive.types]{is_environment}}.
#' @name is_environment
#' @aliases assert_is_environment
#' @export is_environment
#' @export assert_is_environment
NULL

#' Is the input a factor?
#' 
#' See \code{\link[assertive.types]{is_factor}}.
#' @name is_factor
#' @aliases is_ordered assert_is_factor assert_is_ordered
#' @export is_factor
#' @export is_ordered
#' @export assert_is_factor
#' @export assert_is_ordered
NULL

#' Is the input a function?
#' 
#' See \code{\link[assertive.types]{is_function}}.
#' @name is_function
#' @aliases is_primitive is_stepfun assert_is_function assert_is_primitive assert_is_stepfun
#' @export is_function
#' @export is_primitive
#' @export is_stepfun
#' @export assert_is_function
#' @export assert_is_primitive
#' @export assert_is_stepfun
NULL

#' Does the object inherit from some class?
#' 
#' See \code{\link[assertive.types]{is_inherited_from}}.
#' @name is_inherited_from
#' @aliases assert_is_inherited_from
#' @export is_inherited_from
#' @export assert_is_inherited_from
NULL

#' Is the input an integer?
#' 
#' See \code{\link[assertive.types]{is_integer}}.
#' @name is_integer
#' @aliases is_an_integer assert_is_integer assert_is_an_integer
#' @export is_integer
#' @export is_an_integer
#' @export assert_is_integer
#' @export assert_is_an_integer
NULL

#' Is the input a language object?
#' 
#' See \code{\link[assertive.types]{is_language}}.
#' @name is_language
#' @aliases is_call is_expression is_name is_symbol assert_is_language assert_is_call assert_is_expression assert_is_name assert_is_symbol
#' @export is_language
#' @export is_call
#' @export is_expression
#' @export is_name
#' @export is_symbol
#' @export assert_is_language
#' @export assert_is_call
#' @export assert_is_expression
#' @export assert_is_name
#' @export assert_is_symbol
NULL

#' Is the input a list?
#' 
#' See \code{\link[assertive.types]{is_list}}.
#' @name is_list
#' @aliases assert_is_list
#' @export is_list
#' @export assert_is_list
NULL

#' Is the input logical?
#' 
#' See \code{\link[assertive.types]{is_logical}}.
#' @name is_logical
#' @aliases is_a_bool assert_is_logical assert_is_a_bool
#' @export is_logical
#' @export is_a_bool
#' @export assert_is_logical
#' @export assert_is_a_bool
NULL

#' Is the input numeric?
#' 
#' See \code{\link[assertive.types]{is_numeric}}.
#' @name is_numeric
#' @aliases is_double is_a_number is_a_double assert_is_numeric assert_is_double assert_is_a_number assert_is_a_double
#' @export is_numeric
#' @export is_double
#' @export is_a_number
#' @export is_a_double
#' @export assert_is_numeric
#' @export assert_is_double
#' @export assert_is_a_number
#' @export assert_is_a_double
NULL

#' Is the input a QR decomposition of a matrix?
#' 
#' See \code{\link[assertive.types]{is_qr}}.
#' @name is_qr
#' @aliases assert_is_qr
#' @export is_qr
#' @export assert_is_qr
NULL

#' Is the input raw?
#' 
#' See \code{\link[assertive.types]{is_raw}}.
#' @name is_raw
#' @aliases is_a_raw assert_is_raw assert_is_a_raw
#' @export is_raw
#' @export is_a_raw
#' @export assert_is_raw
#' @export assert_is_a_raw
NULL

#' Is the input an S4 object?
#' 
#' See \code{\link[assertive.types]{is_s4}}.
#' @name is_s4
#' @aliases is_s4 is_S4 assert_is_s4 assert_is_S4
#' @export is_s4
#' @export is_S4
#' @export assert_is_s4
#' @export assert_is_S4
NULL

#' Is the input a table?
#' 
#' See \code{\link[assertive.types]{is_table}}.
#' @name is_table
#' @aliases is_table assert_is_table
#' @export is_table
#' @export assert_is_table
NULL

# From assertive.types is-condition.R assert-is-condition.R

#' Is the input a condition?
#' 
#' See \code{\link[assertive.types]{is_try_error}}.
#' @name is_try_error
#' @aliases is_simple_error is_error is_simple_warning is_warning is_simple_message is_message is_condition assert_is_try_error assert_is_simple_error assert_is_error assert_is_simple_warning assert_is_warning assert_is_simple_message assert_is_message assert_is_condition
#' @export is_try_error
#' @export is_simple_error
#' @export is_error
#' @export is_simple_warning
#' @export is_warning
#' @export is_simple_message
#' @export is_message
#' @export is_condition
#' @export assert_is_try_error
#' @export assert_is_simple_error
#' @export assert_is_error
#' @export assert_is_simple_warning
#' @export assert_is_warning
#' @export assert_is_simple_message
#' @export assert_is_message
#' @export assert_is_condition
NULL

# From assertive.types is-date.R  assert-is-date.R

#' Is the input a date?
#' 
#' See \code{\link[assertive.types]{is_date}}.
#' @name is_date
#' @aliases is_posixct is_posixlt assert_is_date assert_is_posixct assert_is_posixlt
#' @export is_date
#' @export is_posixct
#' @export is_posixlt
#' @export assert_is_date
#' @export assert_is_posixct
#' @export assert_is_posixlt
NULL

# From assertive.types is-formula.R, assert-is-formula.R

#' Is the input a formula?
#' 
#' See \code{\link[assertive.types]{is_formula}}.
#' @name is_formula
#' @aliases is_one_sided_formula is_two_sided_formula assert_is_formula assert_is_one_sided_formula assert_is_two_sided_formula
#' @export is_formula
#' @export is_one_sided_formula
#' @export is_two_sided_formula
#' @export assert_is_formula
#' @export assert_is_one_sided_formula
#' @export assert_is_two_sided_formula
NULL

# From assertive.types is-type-data.table.R, assert-is-type-data.table.R

#' Is the input a data.table?
#' 
#' See \code{\link[assertive.types]{is_data.table}}.
#' @name is_data.table
#' @aliases assert_is_data.table
#' @export is_data.table
#' @export assert_is_data.table
NULL

# From assertive.types is-type-dplyr.R, assert-is-type-dplyr.R

#' Is the input a tbl?
#' 
#' See \code{\link[assertive.types]{is_tbl}}.
#' @name is_tbl
#' @aliases is_tbl_cube is_tbl_df is_tbl_dt assert_is_tbl assert_is_tbl_cube assert_is_tbl_df assert_is_tbl_dt
#' @export is_tbl
#' @export is_tbl_cube
#' @export is_tbl_df
#' @export is_tbl_dt
#' @export assert_is_tbl
#' @export assert_is_tbl_cube
#' @export assert_is_tbl_df
#' @export assert_is_tbl_dt
NULL

# From assertive.types is-type-grDevices.R, assert-is-type-grDevices.R

#' Is the input a raster?
#' 
#' See \code{\link[assertive.types]{is_raster}}.
#' @name is_raster
#' @aliases assert_is_raster
#' @export is_raster
#' @export assert_is_raster
NULL

# From assertive.types is-type-methods.R, assert-is-type-methods.R

#' Is the input the name of a (formally defined) class?
#' 
#' See \code{\link[assertive.types]{is_class}}.
#' @name is_class
#' @aliases assert_all_are_classes assert_any_are_classes
#' @export is_class
#' @export assert_all_are_classes
#' @export assert_any_are_classes
NULL

# From assertive.types is-type-stats.R, assert-is-type-stats.R

#' Is the input a (dendrogram) leaf?
#' 
#' See \code{\link[assertive.types]{is_leaf}}.
#' @name is_leaf
#' @aliases is_leaf assert_is_leaf
#' @export is_leaf
#' @export assert_is_leaf
NULL

#' Is the input a time series?
#' 
#' See \code{\link[assertive.types]{is_ts}}.
#' @name is_ts
#' @aliases is_mts is_tskernel assert_is_ts assert_is_mts assert_is_tskernel
#' @export is_ts
#' @export is_mts
#' @export is_tskernel
#' @export assert_is_ts
#' @export assert_is_mts
#' @export assert_is_tskernel
NULL

# From assertive.types is-type-utils.R, assert-is-type-utils.R

#' Is the input relistable?
#' 
#' See \code{\link[assertive.types]{is_relistable}}.
#' @name is_relistable
#' @aliases assert_is_relistable
#' @export is_relistable
#' @export assert_is_relistable
NULL
