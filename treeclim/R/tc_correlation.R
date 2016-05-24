##' (bootstrapped) correlation function
##' 
##' See documentation of dcc for details.
##' @param chrono a tree-ring chronology
##' @param climate data.frame with climate parameters
##' @param ci numeric: p-level for confidence interval (must be in
##' c(0.1, 0.05, 0.01)
##' @param boot character, which bootstrapping method should be used
##' (or "none")
##' @param p probability for rgeom, that determines distribution of
##' sampling blocks for stationary bootstrap scheme
##' @keywords internal
tc_correlation <- function(chrono, climate, ci, boot) {

  vnames <- climate$names
  n <- length(chrono)
  m <- dim(climate$aggregate)[2]
  
  boot_data <- init_boot_data(as.matrix(climate$aggregate),
                              chrono, 1000, boot)
  
  if (boot %in% c("stationary", "std", "dendroclim")) {
    param_matrix <- .Call("treeclim_corfun", boot_data$climate,
                          boot_data$chrono, PACKAGE = "treeclim")$coef
    
    out <- ptest(param_matrix, ci, NULL, "range")
  }
  
  if (boot == "exact") {
    param_matrix <- .Call("treeclim_corfunexact", boot_data$climate,
                          boot_data$chrono, chrono, PACKAGE = "treeclim")$coef
    
    out <- ptest(param_matrix[,2:1001], ci, param_matrix[,1], "weibull")
  }
  
  if (boot == "none") {
    param_matrix <- .Call("treeclim_corfunnoboot", boot_data$climate,
                          boot_data$chrono, PACKAGE = "treeclim")$coef
    
    out <- ptest(param_matrix, ci, NULL, "none")  
  }
  
  rownames(out) <- vnames
  attributes(out)$npar <- attributes(climate$aggregate)$npar
  attributes(out)$vnames <- vnames
  
  ## include information for pretty printing and assemble output data.frame
  
  out <- list(
    result = data.frame(
      id = climate$pretty_names$id,
      varname = climate$pretty_names$varname,
      month = climate$pretty_names$month_label,
      out
    ),
    ac = boot_data$ac)
  
  class(out$result) <- c("tc_coef", "data.frame")
  out
}
