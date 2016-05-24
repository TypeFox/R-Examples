#' Interactive table of results
#' 
#' A summary of past results
#' @inheritParams plot_past
#' @export
#' @examples 
#' ## Need the DT package
#' ## View all results for prog test
#' get_datatable_past()
#' ## View matrix_fun test
#' get_datatable_past("matrix_fun")
#' ## View matrix_fun test - only BLAS results
#' get_datatable_past("matrix_fun", blas_optimize=TRUE)
get_datatable_past = function(test_group=c("prog", "matrix_fun", "matrix_cal", 
                                           "read", "write"), 
                              byte_optimize=NULL, blas_optimize=NULL) {
  if(!requireNamespace("DT", quietly = TRUE))
    stop("Install DT package to use datatable")

  results = select_results(test_group, byte_optimize, blas_optimize)
  results$time = signif(results$time, 4)
  
  ## Format
  results = results[,c("rank", "time", "cpu", "byte_optimize", "blas_optimize", "sysname", "test_group")]
  colnames(results) = c("Rank", "Time (sec)", "CPU", "Byte Compile", "BLAS Opt", "OS", "Test")
  results = results[,c(TRUE, TRUE, TRUE, is.null(byte_optimize), is.null(blas_optimize), 
                       TRUE, length(unique(results$Test)) > 1)]
  DT::datatable(results, rownames=FALSE) 
}
