# From assertive.base are-same-size.R assert-are-same-size.R

#' Are the inputs the same length/dimension?
#' 
#' See \code{\link[assertive.properties]{are_same_length}}.
#' @name are_same_length
#' @aliases have_same_dims are_same_length_legacy assert_are_same_length assert_have_same_dims assert_all_are_same_length_legacy assert_any_are_same_length_legacy assert_all_are_same_length assert_any_are_same_length
#' @export are_same_length
#' @export are_same_length_legacy
#' @export have_same_dims
#' @export assert_are_same_length
#' @export assert_all_are_same_length_legacy
#' @export assert_any_are_same_length_legacy
#' @export assert_all_are_same_length
#' @export assert_any_are_same_length
#' @export assert_have_same_dims
NULL 

# From assertive.properties is-monotonic.R, assert-is-monotonic.R

#' Is the vector monotonically increasing or decreasing?
#' 
#' See \code{\link[assertive.properties]{is_monotonic_increasing}}.
#' @name is_monotonic_increasing
#' @aliases is_monotonic is_monotonic_decreasing assert_is_monotonic_increasing assert_is_monotonic_decreasing
#' @export is_monotonic_increasing
#' @export is_monotonic_decreasing
#' @export assert_is_monotonic_increasing
#' @export assert_is_monotonic_decreasing
NULL

# From assertive.properties has-attributes.R, assert-has-attributes.R

#' Does the input have any attributes?
#' 
#' See \code{\link[assertive.properties]{has_any_attributes}}.
#' @name has_any_attributes
#' @aliases has_no_attributes
#' @export has_any_attributes
#' @export has_no_attributes
NULL

#' Does the input have the specified attributes?
#' 
#' See \code{\link[assertive.properties]{has_attributes}}.
#' @name has_attributes
#' @aliases assert_has_all_attributes assert_has_any_attributes
#' @export has_attributes
#' @export assert_has_all_attributes
#' @export assert_has_any_attributes
NULL

# From assertive.properties has-dims.R, assert-has-dims.R

#' Does the input have rows/columns?
#' 
#' See \code{\link[assertive.properties]{has_cols}}.
#' @name has_cols
#' @aliases has_rows assert_has_cols assert_has_rows
#' @export has_cols
#' @export has_rows
#' @export assert_has_cols
#' @export assert_has_rows
NULL

#' Does the input have dimensions?
#' 
#' See \code{\link[assertive.properties]{has_dims}}.
#' @name has_dims
#' @aliases assert_has_dims
#' @export has_dims
#' @export assert_has_dims
NULL

# From assertive.properties has-dupes.R, assert-has-dupes.R

#' Does the input have duplicates?
#' 
#' See \code{\link[assertive.properties]{has_duplicates}}.
#' @name has_duplicates
#' @aliases has_no_duplicates assert_has_duplicates assert_has_no_duplicates
#' @export has_duplicates
#' @export has_no_duplicates
#' @export assert_has_duplicates
#' @export assert_has_no_duplicates
NULL

# From assertive.properties has-names.R, assert-has-names.R

#' Does the input have names?
#' 
#' See \code{\link[assertive.properties]{has_names}}.
#' @name has_names
#' @aliases has_rownames has_colnames has_dimnames assert_has_names assert_has_rownames assert_has_colnames assert_has_dimnames 
#' @export has_names
#' @export has_rownames
#' @export has_colnames
#' @export has_dimnames
#' @export assert_has_names
#' @export assert_has_rownames
#' @export assert_has_colnames
#' @export assert_has_dimnames
NULL

# From assertive.properties has-slot.R, assert-has-slot.R

#' Does the S4 input have a slot?
#' 
#' See \code{\link[assertive.properties]{has_slot}}.
#' @name has_slot
#' @aliases assert_has_slot
#' @export has_slot
#' @export assert_has_slot
NULL

# From assertive.properties is-atomic-recursive-vector.R, assert-is-atomic-recursive-vector.R

#' Is the input atomic/recursive/vector?
#' 
#' See \code{\link[assertive.properties]{is_atomic}}.
#' @name is_atomic
#' @aliases is_recursive is_vector is_nested is_non_nested assert_is_atomic assert_is_recursive assert_is_vector assert_is_nested assert_is_non_nested
#' @export is_atomic
#' @export is_recursive
#' @export is_vector
#' @export is_nested
#' @export is_non_nested
#' @export assert_is_atomic
#' @export assert_is_recursive
#' @export assert_is_vector
#' @export assert_is_nested
#' @export assert_is_non_nested
NULL

# From assertive.properties is-empty-scalar.R, assert-is-empty-scalar.R

#' Is the input empty/scalar?
#' 
#' See \code{\link[assertive.properties]{is_empty}}.
#' @name is_empty
#' @aliases is_non_empty is_scalar is_non_scalar has_elements is_of_dimension is_of_length assert_is_empty assert_is_non_empty assert_is_scalar assert_is_non_scalar assert_has_elements assert_is_of_dimension assert_is_of_length
#' @export is_empty
#' @export is_non_empty
#' @export is_scalar
#' @export is_non_scalar
#' @export has_elements
#' @export is_of_dimension
#' @export is_of_length
#' @export assert_is_empty
#' @export assert_is_non_empty
#' @export assert_is_scalar
#' @export assert_is_non_scalar
#' @export assert_has_elements
#' @export assert_is_of_dimension
#' @export assert_is_of_length
NULL

# From assertive.properties is-null.R, assert-is-null.R

#' Checks to see if the input is (not) null.
#' 
#' See \code{\link[assertive.properties]{is_null}}.
#' @name is_null
#' @aliases is_not_null assert_is_null assert_is_not_null
#' @export is_null
#' @export is_not_null
#' @export assert_is_null
#' @export assert_is_not_null
NULL

# From assertive.properties is-unsorted.R, assert-is-unsorted.R

#' Is the input unsorted?
#' 
#' See \code{\link[assertive.properties]{is_unsorted}}.
#' @name is_unsorted
#' @aliases assert_is_unsorted
#' @export is_unsorted
#' @export assert_is_unsorted
NULL

# From assertive.properties utils.R

#' Get the dimensions of an object
#' 
#' See \code{\link[assertive.properties]{DIM}}.
#' @name DIM
#' @export 
NULL

#' Get the number of elements
#' 
#' See \code{\link[assertive.properties]{n_elements}}.
#' @name n_elements
#' @export 
NULL
