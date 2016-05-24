#' query \code{dimVar}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{dimVar}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item varName: variable name of the variable from which \code{object} was calculated
#' \item codesOriginal: original codes (as specified by the user)
#' \item codesDefault: calculated, default codes
#' \item codesMinimal: all codes required to calculate the complete hierarchy (no sub-totals)
#' \item levels: level-structure of the dimensional variable
#' \item structure: vector showing how many digits in the default codes are required for each level
#' \item dims: list showing the complete hierarchy of the dimensional variable
#' \item dups: vector of duplicated codes
#' \item dupsUp: vector of codes that are the 'upper' levels to which the codes in \code{dups} correspond
#' \item hasDuplicates: does the dimensional variable has codes that can be (temporarily) removed
#' \item nrLevels: the total number of levels of a dimensional variable
#' \item minimalCodesDefault: the standardized codes of the minimal set of required level-codes
#'
#' @return information from objects of class \code{dataObj} depending on argument \code{type}
#' \itemize{
#' \item a list if argument \code{type} matches 'dims'
#' \item numeric vector if argument \code{type} matches 'levels' or 'nrLevels'
#' \item character vector if argument \code{type} matches 'codesOriginal', 'codesDefault', 'vName', 'dups', 'dupsUp' or 'minimalCodesDefault'
#' \item logical vector of length 1 if argument \code{type} matches 'hasDuplicates'
#' \item a logical vector if argument \code{type} matches 'codesMinimal'}
#'
#' @export
#' @docType methods
#' @rdname get.dimVar-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("get.dimVar", function(object, type) { standardGeneric("get.dimVar") })

#' modify \code{dimVar}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{dimVar}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item hasDefaultCodes: calculates if a vector of codes (specified by argument \code{input}) corresponds to default codes in \code{object}
#' \item matchCodeOrig: obtain default|standard codes for a vector of original codes specified by argument \code{input}
#' \item matchCodeDefault: obtain original codes for a vector of default|standard codes specified by argument \code{input}
#' \item standardize: perform standardization of level-codes (temporarily removing duplicates,..)
#' \item requiredMinimalCodes: calculate a set of minimal codes that are required to calculate a specific (sub)total specified by argument \code{input}
#' @param input a character vector
#'
#' @return information from \code{object} depending on \code{type}
#' \itemize{
#' \item a character vector if type matches 'matchCodeOrig', 'matchCodeDefault', 'standardize' or 'requiredMinimalCodes'
#' \item a logical vector of length 1 if type matches 'hasDefaultCodes' being TRUE if argument \code{input} are default codes and FALSE otherwise
#' }
#'
#' @export
#' @docType methods
#' @rdname calc.dimVar-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("calc.dimVar", function(object, type, input) { standardGeneric("calc.dimVar") })

#' initialize \code{dimVar}-object
#'
#' @param input a list with 2 elements}
#' \itemize{
#' \item first element: either an object of class 'matrix' or a data.frame or a link to a file. The input data need to be in a specific format (2 columns) with the first column defining the level-structure and the second column defining the level-codes.
#' \item second element: a character vector of length 1 specifying a variable name
#'
#' @export
#' @docType methods
#' @rdname init.dimVar-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric("init.dimVar", function(input) { standardGeneric("init.dimVar") })

# get-methods
setGeneric("g_varname", function(object) { standardGeneric("g_varname") })
setGeneric("g_original_codes", function(object) { standardGeneric("g_original_codes") })
setGeneric("g_default_codes", function(object) { standardGeneric("g_default_codes") })
setGeneric("g_minimal_codes", function(object) { standardGeneric("g_minimal_codes") })
setGeneric("g_levels", function(object) { standardGeneric("g_levels") })
setGeneric("g_structure", function(object) { standardGeneric("g_structure") })
setGeneric("g_dims", function(object) { standardGeneric("g_dims") })
setGeneric("g_dups", function(object) { standardGeneric("g_dups") })
setGeneric("g_dups_up", function(object) { standardGeneric("g_dups_up") })
setGeneric("g_has_dups", function(object) { standardGeneric("g_has_dups") })
setGeneric("g_nr_levels", function(object) { standardGeneric("g_nr_levels") })
setGeneric("g_minimal_default_codes", function(object) { standardGeneric("g_minimal_default_codes") })

# calc-methods
setGeneric("c_has_default_codes", function(object, input) { standardGeneric("c_has_default_codes") })
setGeneric("c_match_orig_codes", function(object, input) { standardGeneric("c_match_orig_codes") })
setGeneric("c_match_default_codes", function(object, input) { standardGeneric("c_match_default_codes") })
setGeneric("c_standardize", function(object, input) { standardGeneric("c_standardize") })
setGeneric("c_required_minimal_codes", function(object, input) { standardGeneric("c_required_minimal_codes") })

