#' query \code{simpleTriplet}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{simpleTriplet}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item rowInd: extract all row-indices
#' \item colInd: extract all column-indices
#' \item values: extract all values
#' \item nrRows: return the number of rows of the input object
#' \item nrCols: return the number of columns of the input object
#' \item nrCells: return the number of cells (different from 0!)
#' \item duplicatedRows: return a numeric vector showing indices of duplicated rows
#' \item transpose: transpose input \code{object} and return the transposed matrix
#' \item getRow: return a specific row of input \code{object}
#' \item getCol: return a specific column of input \code{object}
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item type == 'getRow': input is a list of length 1
#' \itemize{
#' \item first element: numeric vector of length 1 defining index of row that is to be returned }
#' \item type == 'getCol': input is a list of length 1
#' \itemize{
#' \item first element: numeric vector of length 1 defining index of column that is to be returned }
#' \item else: input is not used at all (empty list)
#' @return information from \code{object} depending on \code{type}
#' \itemize{
#' \item a numeric vector if type matches 'rowInd', 'colInd', 'values', 'nrRows', 'nrCols', 'nrCells' or 'duplicatedRows'
#' \item an object of class \code{simpleTriplet} if type matches 'transpose', 'getRow' or 'getCol'
#' }
#'
#' @export
#' @docType methods
#' @rdname get.simpleTriplet-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('get.simpleTriplet', function(object, type, input) {standardGeneric('get.simpleTriplet')})

#' modify \code{simpleTriplet}-objects depending on argument \code{type}
#'
#' @param object an object of class \code{simpleTriplet}
#' @param type a character vector of length 1 defining what to calculate|return|modify. Allowed types are:}
#' \itemize{
#' \item removeRow: remove a row with given index from \code{object}
#' \item removeCol: remove a column with given index from \code{object}
#' \item addRow: add a row to \code{object}
#' \item addCol: add a column to \code{object}
#' \item modifyRow: change specified row of \code{object}
#' \item modifyCol: change specified column of \code{object}
#' \item modifyCell: change specified cell of \code{object}
#' \item bind: bind two objects of class \code{simpleTriplet} together
#' @param input a list depending on argument \code{type}.}
#'
#' \itemize{
#' \item type==removeRow: input is a list of length 1
#' \itemize{
#' \item first element: numeric vector of length 1 defining the index of the row that should be removed }
#' \item type==removeCol: input is a list of length 1
#' \itemize{
#' \item first element: numeric vector of length 1 defining the index of the column that should be removed }
#' \item type==addRow: input is a list of length 2
#' \itemize{
#' \item first element: numeric vector of column-indices
#' \item second element: numeric vector defining the cell-values of the row that will be added }
#' \item type==addCol: input is a list of length 2
#' \itemize{
#' \item first element: numeric vector of row-indices
#' \item second element: numeric vector defining the cell-values of the column that will be added }
#' \item type==modifyRow: input is a list of length 3
#' \itemize{
#' \item first element: numeric vector of length 1 specifying the the row-index of the row that will be modified
#' \item second element: numeric vector specifying the column-indices that should be modified
#' \item third element: numeric vector defining values that should be set in the given row }
#' \item type==modifyCol: input is a list of length 3
#' \itemize{
#' \item first element: numeric vector specifying the row-indices that should be modified
#' \item second element: numeric vector of length 1 specifying the the column-index of the column that will be modified
#' \item third element: numeric vector defining values that should be set in the given column }
#' \item type==modifyCell: input is a list of length 3
#' \itemize{
#' \item first element: numeric vector of length 1 defining the column-index
#' \item second element: numeric vector of length 1 defining the row-index
#' \item third element: numeric vector of length 1 holding the value that should be set in the given cell }
#' \item type==bind: input is a list of length 2
#' \itemize{
#' \item first element: an object of class \code{simpleTriplet}
#' \item second argument: is a logical vector of length 1 being TRUE if a 'rbind' or 'FALSE' if a 'cbind' should be done }
#'
#' @return an object of class \code{simpleTriplet}
#'
#' @export
#' @docType methods
#' @rdname calc.simpleTriplet-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('calc.simpleTriplet', function(object, type, input) {standardGeneric('calc.simpleTriplet')})

#' initialize \code{simpleTriplet}-objects depending on argument \code{type}
#'
#' init.simpleTriplet should be used to create objects of class \code{simpleTriplet}.
#' It is possible to create an object from class \code{simpleTriplet} from an existing matrix (using type=='simpleTriplet').
#' A positive (or negative) identity matrix stored as an object of class \code{simpleTriplet} can be created by specifying type=='simpleTripletDiag'.
#'
#' @param type a character vector of length 1 defining what|how to initialize. Allowed types are:}
#' \itemize{
#' \item simpleTriplet: a simple triplet matrix
#' \item simpleTripletDiag: identity matrix
#' @param input a list depending on argument \code{type}.}
#' \itemize{
#' \item type == 'simpleTriplet': input is a list of length 1
#' \itemize{
#' \item first element: object of class 'matrix' }
#' \item type == 'simpleTripletDiag': input is a list of length 2
#' \itemize{
#' \item first element: numeric vector of length 1 defining the desired number of rows of the identiy matrix
#' \item second element: logical vector of length 1 being TRUE if a positive and FALSE if a negative identity matrix should be returned }
#' @return an object of class \code{simpleTriplet}
#'
#' @export
#' @docType methods
#' @rdname init.simpleTriplet-method
#'
#' @note internal function
#' @author Bernhard Meindl \email{bernhard.meindl@@statistik.gv.at}
setGeneric('init.simpleTriplet', function(type, input) {standardGeneric('init.simpleTriplet')})

# get-methods
setGeneric("g_row_ind", function(object) { standardGeneric("g_row_ind") })
setGeneric("g_col_ind", function(object) { standardGeneric("g_col_ind") })
setGeneric("g_values", function(object) { standardGeneric("g_values") })
setGeneric("g_nr_rows", function(object) { standardGeneric("g_nr_rows") })
setGeneric("g_nr_cols", function(object) { standardGeneric("g_nr_cols") })
setGeneric("g_nr_cells", function(object) { standardGeneric("g_nr_cells") })
setGeneric("g_duplicated_rows", function(object) { standardGeneric("g_duplicated_rows") })
setGeneric("g_transpose", function(object) { standardGeneric("g_transpose") })
setGeneric("g_row", function(object, input) { standardGeneric("g_row") })
setGeneric("g_col", function(object, input) { standardGeneric("g_col") })

# calc-methods
setGeneric("c_remove_row", function(object, input) { standardGeneric("c_remove_row") })
setGeneric("c_remove_col", function(object, input) { standardGeneric("c_remove_col") })
setGeneric("c_add_row", function(object, input) { standardGeneric("c_add_row") })
setGeneric("c_add_col", function(object, input) { standardGeneric("c_add_col") })
setGeneric("c_modify_row", function(object, input) { standardGeneric("c_modify_row") })
setGeneric("c_modify_col", function(object, input) { standardGeneric("c_modify_col") })
setGeneric("c_modify_cell", function(object, input) { standardGeneric("c_modify_cell") })
setGeneric("c_bind", function(object, input) { standardGeneric("c_bind") })

