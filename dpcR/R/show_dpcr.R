#' Methods for Function \code{show}
#' 
#' Expands function \code{\link[methods]{show}} allowing showing objects of the
#' class \code{\linkS4class{adpcr}} or \code{\linkS4class{ddpcr}}.
#' 
#' 
#' @name show-methods
#' @aliases show-methods show-methods show,dpcr-method show,dpcr-method show
#' @docType methods
#' @param object an object of class \code{\linkS4class{dpcr}}.
#' @author Michal Burdukiewicz.
#' @keywords methods utilities
#' @examples
#' 
#' #array dpcr
#' ptest <- sim_adpcr(400, 765, 5, FALSE, n_panels = 1)
#' show(ptest)
#' 
#' #multiple experiments
#' ptest <- sim_adpcr(400, 765, 5, FALSE, n_panels = 5)
#' show(ptest)
#' 
#' #droplet dpcr - fluorescence
#' dropletf <- sim_ddpcr(7, 20, times = 5, fluo = list(0.1, 0))
#' show(dropletf)
#' 
#' #droplet dpcr - number of molecules
#' droplet <- sim_ddpcr(7, 20, times = 5)
#' show(droplet)
#' 
NULL

# Special method declared to hide slots other than .Data
setMethod("show", signature(object = "dpcr"), function(object) {
  if(nrow(slot(object, ".Data")) > 5) {
    print(slot(object, ".Data")[1L:5, , drop = FALSE])
    cat(paste0("\n", nrow(slot(object, ".Data")) - 5, " data points ommited."))
  } else {
    print(slot(object, ".Data"))
  }
  cat(paste0("\nData type: '", slot(object, "type"), "'"))
})
