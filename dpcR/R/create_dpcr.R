#' Create dpcR object
#' 
#' Creates \code{\linkS4class{adpcr}} and \code{\linkS4class{ddpcr}}
#' objects from data.
#' 
#' @details 
#' This constructor function assists in creation of objects used by other 
#' functions of the package. It is also responsible for checking the correctness of
#' arguments.
#' 
#' A warning is prompted whenever any of arguments is converted to other type.
#' 
#' @note 
#' \code{create_dpcr} is a preferred to calling directly 
#' \code{\link[methods]{new}}.
#' 
#' Currently only end-point measurements are supported.
#' 
#' @param data a \code{numeric} vector or matrix of data from dPCR
#' experiments. Data frames will be converted to matrices.
#' @param n \code{integer} equal to number of partitions.
#' @param exper The id of experiments.
#' @param replicate The id of technical replicates.
#' @param assay The name or id of assays.
#' @param threshold \code{numeric} value giving the threshold above which
#' droplet is counted as positive.  Ignored if \code{adpcr} is \code{TRUE}.
#' @param breaks \code{numeric} vector giving the number of intervals into
#' which \code{data} should be cut. Ignored if \code{adpcr} is \code{FALSE}.
#' @param type Object of class \code{"character"} defining type of data. Could
#' be \code{"nm"} (number of molecules per partition), \code{"tnp"} (total
#' number of positive wells in the panel), \code{"fluo"} (fluorescence), \code{"np"} 
#' (status (positive (1) or negative(0)) of each droplet) or\code{"ct"} 
#' (threshold cycle).
#' @param adpcr \code{logical}. If \code{TRUE}, function creates
#' \code{\linkS4class{adpcr}} object. If \code{FALSE}, function creates
#' \code{\linkS4class{ddpcr}} object.
#' @param col_names \code{character} vector of column names in array. Ignored if not
#' \code{adcpr}.
#' @param row_names \code{character} vector of row names in array. Ignored if not
#' \code{adcpr}.
#' @param panel_id \code{factor} vector of panel IDs (or names). Ignored if not
#' \code{adcpr}.
#' @return An \code{\linkS4class{adpcr}} or \code{\linkS4class{ddpcr}} object.
#' @author Michal Burdukiewicz, Stefan Roediger.
#' @keywords ddPCR adPCR
#' @examples
#' 
#' # Droplet digital PCR example
#' sample_runs <- matrix(rpois(60, lambda = 1.5), ncol = 2)
#' ddpcr1 <- create_dpcr(sample_runs[,1], n = 30L, 
#' threshold = 1, type = "nm", adpcr = FALSE)
#' ddpcr2 <- create_dpcr(sample_runs[,2], n = 30L, 
#' threshold = 1, type = "nm", adpcr = FALSE)
#' plot_vic_fam(ddpcr1, ddpcr2)
#' 
#' # Array digital PCR example
#' sample_adpcr <- create_dpcr(rpois(765, lambda = 0.8), n = 765L, 
#' 			    type = "nm", adpcr = TRUE)
#' plot_panel(sample_adpcr, 45, 17)
#' 
#' @export create_dpcr

create_dpcr <- function(data, n, exper = "Experiment 1", 
                        replicate = NULL, assay = "Unknown", type, threshold = NULL,
                        breaks = NULL, adpcr, col_names = NULL, row_names = NULL, panel_id = NULL) {
  if(adpcr) {
    create_adpcr(data = data, n = n, exper = exper, 
                 replicate = replicate, assay = assay, type = type, breaks = breaks,
                 col_names = col_names, row_names = row_names, panel_id = panel_id)
  } else {
    create_ddpcr(data = data, n = n, exper = exper, replicate = replicate, 
                 assay = assay, type = type, threshold = threshold)
  }
}