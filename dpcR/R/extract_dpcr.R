#extract single panel from dpcr object


#' Extract Digital PCR Experiment
#' 
#' Extract panel(s) or experiment(s) from a matrix while preserving all other
#' attributes.
#' 
#' The \code{extract_dpcr} function allows to choose one or more panels from an
#' object of the \code{\linkS4class{adpcr}} or \code{\linkS4class{ddpcr}} class
#' and save it without changing other attributes. It is the most recommended
#' method of extracting a subset from an array of panels, because it preserves
#' class and structure of the object in contrary to standard operator
#' \link[base]{Extract}.
#' 
#' @param input object of the class \code{\linkS4class{adpcr}} or
#' \code{\linkS4class{ddpcr}}.
#' @param id vector of indices or names of runs.
#' @return The object of the input's class (\code{\linkS4class{adpcr}} or
#' \code{\linkS4class{ddpcr}}).
#' @note The standard \code{\link[base]{Extract}} operator \code{x[i]} treats
#' dpcr objects as \code{matrix} and extracts values without preserving other
#' attributies of the object.
#' @author Michal Burdukiewicz.
#' @seealso Opposite function: \code{\link{bind_dpcr}}
#' @keywords manip extract panel
#' @examples
#' 
#' #sample extracting
#' panels <- sim_adpcr(10, 40, 1000, pos_sums = FALSE, n_panels = 50)
#' single_panel <- extract_dpcr(panels, 5)
#' random_three <- extract_dpcr(panels, sample.int(ncol(panels), 3))
#' all_but_one <- extract_dpcr(panels, -5)
#' 
#' #the same for fluorescence data
#' fluos <- sim_ddpcr(10, 40, 1000, pos_sums = FALSE, n_exp = 50, 
#'                    fluo = list(0.1, 0))
#' single_fluo <- extract_dpcr(fluos, 5)
#' 
#' 
#' @export extract_dpcr
extract_dpcr <- function(input, id) {
  if (!(class(input) %in% c("adpcr", "ddpcr")))
    stop("Input must have 'adpcr' or 'ddpcr' class.")
  
  #when id is a column name
  if(!is.numeric(id)) 
    id <- which(colnames(input) %in% id) 

  selected <- input[, id, drop = FALSE]
  
  #because when id is single negative value, usually the
  #result has more than one column
  if (length(id) == 1 && id > 0) {
    selected <- matrix(selected, ncol = 1)
    colnames(selected) <- colnames(input)[id]
  }
  
  result <- input
  slot(result, ".Data") <- selected
  slot(result, "n") <- slot(input, "n")[id]
  slot(result, "exper") <- droplevels(slot(input, "exper")[id])
  slot(result, "replicate") <- droplevels(slot(input, "replicate")[id])
  slot(result, "assay") <- droplevels(slot(input, "assay")[id])
  
  
  #in case of tnp type extract also columns names
  if (class(input) == "adpcr") {
    slot(result, "panel_id") <- droplevels(slot(input, "panel_id")[id])
    if (slot(input, "type") == "tnp")
      slot(result, "col_names") <- slot(input, "col_names")[id]
  }
  
  result
}
