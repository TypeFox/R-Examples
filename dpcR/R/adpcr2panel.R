#' Convert adpcr object to array
#' 
#' Converts \code{\linkS4class{adpcr}} object into the list of array-like matrices.
#' 
#' @param input object of the \code{\linkS4class{adpcr}} class.
#' @param use_breaks if \code{TRUE}, input is cutted into intervals using the 
#' \code{breaks} slot. If \code{FALSE}, the \code{integer} class of the input is
#' preserved. Ignored if data has \code{"np"} type (see possible types of 
#' \code{\linkS4class{adpcr}} objects).
#' @return A named list of length equal to the number of arrays in the \code{input}. 
#' Each element is a single array in matrix-like form, where dimensions are set 
#' exactly as in case of the real plate. Names of the list corresponds to the names 
#' of assays (\code{"tnp"} data) or runs (any other type of \code{\linkS4class{adpcr}} 
#' data).
#' 
#' The matrices contain values from array, either integers (when \code{use_break} is
#' \code{FALSE}) or characters (when \code{use_break} is \code{TRUE}).
#' @author Michal Burdukiewicz.
#' @export
#' @keywords manip
#' @examples 
#' #generate data
#' ttest <- sim_adpcr(m = 400, n = 765, times = 20, pos_sums = FALSE, 
#'                    n_panels = 3)
#' #convert object into three arrays
#' arrays <- adpcr2panel(ttest)
#' length(arrays)
#' #see the first array
#' arrays[[1]]
#' 
#' #convert the object using breaks
#' arrays <- adpcr2panel(ttest, use_breaks = TRUE)
#' arrays[[1]]
#' 
adpcr2panel <- function(input, use_breaks = FALSE) {
  if (class(input) == "adpcr") {
    if (!(slot(input, "type") %in% c("nm", "np", "tnp", "fluo", "ct")))
      stop("Input must contain data of type 'nm', 'np', 'tnp', 'fluo' or 'ct'.") 
  } else {
    stop("Input must have the 'adpcr' class")
  }
  
  nx_a <- length(slot(input, "col_names"))
  ny_a <- length(slot(input, "row_names"))
  
  #in case of tnp, we analyze all experiments (all columns)
  #in the all other case, we analyze only a single value of n
  #assumption - number of experiments in each panel is the same
  len_n <- ifelse(slot(input, "type") == "tnp", table(slot(input, "panel_id"))[1], slot(input, "n"))
  
  if (len_n != nx_a * ny_a)
    stop(paste0("Input length (", len_n, ") differs from the array size (", 
                nx_a * ny_a, ")."))
  
  #apply in case input contains more than 1 array
  #here list of?
  
  array_data <- lapply(levels(slot(input, "panel_id")), function(single_level) {
    #data for a single assay
    single_panel <- extract_dpcr(input, which(slot(input, "panel_id") == single_level))
    # Use breaks points to split input 
    
    if (slot(input, "type") == "np")
      use_breaks = FALSE
    
    if(use_breaks)
      single_panel <- as.character(cut(as.vector(single_panel), breaks = slot(input, "breaks"), 
                                       include.lowest = TRUE, right = FALSE, dig.lab = 5))
    
    matrix(single_panel, ncol = nx_a, 
           dimnames = list(slot(input, "row_names"), slot(input, "col_names")))
  })
  
  if(slot(input, "type") == "tnp") {
    names(array_data) <- levels(slot(input, "panel_id"))
  } else {
    names(array_data) <- colnames(input)
  }
  
  array_data
}