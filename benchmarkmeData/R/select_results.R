#' Convenience function for selecting
#' 
#' Selects, aggregates over the \code{past_results} data set.
#' @return A data frame
#' @inheritParams plot_past
#' @export
#' @examples 
#' select_results("prog", byte_optimize=TRUE)
select_results = function(test_group, byte_optimize=NULL, blas_optimize=NULL) {
  ## Load past data
  tmp_env = new.env()
  data(past_results, package="benchmarkmeData", envir = tmp_env)
  results = tmp_env$past_results
  
  if(!is.null(byte_optimize)) {
    if(byte_optimize) {
      results = results[results$byte_optimize > 0.5,]
    } else {
      results = results[results$byte_optimize < 0.5,]
    }
  }
  
  if(!is.null(blas_optimize)) {
    results = results[results$blas_optimize==blas_optimize,]
  }
  
  if(length(test_group) > 1) test_group = test_group[1]
  if(test_group %in% c("read", "write")) test_group = paste0(test_group, c(5, 50, 200))
  
  results = results[results$test_group %in% test_group,]
  
  ## Aggregate over test
  ## Ensure that we have timings for all required tests.
  results = aggregate(time ~ id + byte_optimize + cpu + date + sysname + blas_optimize + test_group, 
                      data=results, 
                      FUN=sum)
  results = results[!is.na(results$time), ]
  results = results[order(results$time), ]
  results$test_group = factor(results$test_group, levels=test_group)
  results = results[order(results$time), ]
  results$rank = 1:nrow(results)
  
  results
}