tag_results = function(results){
  if(!requireNamespace("benchmarkme", quietly = TRUE)) 
    stop("Install benchmarkme package to use shine.")

  if(is.null(results)) return(NULL)
  
  no_of_reps = nrow(results)/length(unique(results$test))
  times = tapply(results[,3], results[5], sum)/no_of_reps
  
  data.frame(id = -1, data = Sys.Date(), 
             time = as.vector(times), 
             test_group = names(times),
             cpu = benchmarkme::get_cpu()$model_name,
             ram = as.numeric(benchmarkme::get_ram()), 
             byte_optimize = benchmarkme::get_byte_compiler(),
             r_major = R.version$major, r_minor= R.version$minor,
             sysname = as.vector(Sys.info()["sysname"]),
             release= as.vector(Sys.info()["release"]),
             blas_optimize = is_blas_optimize(results), 
             row.names=1:length(times))
}

#' Shiny interface
#'
#' This function uses shiny to explore past benchmarks.
#' If you pass a benchmark, say by running \code{benchmark_std()} your
#' particular machine will be highlighted with a red dot.
#'
#' @param results Benchmark results, probably obtained via the 
#' \code{benchmarkme} package.
#' @examples 
#' \dontrun{
#' shine()
#' ## Example results
#' data(sample_results, package="benchmarkme")
#' shine(sample_results)
#' }
#' @export
shine = function(results=NULL){
  if(!requireNamespace("shiny", quietly = TRUE)) 
    stop("Install shiny package to use shine.")
  if(!requireNamespace("ggplot2", quietly = TRUE))
    stop("Install ggplot2 package to use shine.")
  if(!is.null(results))
    .bme_env$results = tag_results(results)  
  
  appDir = system.file("shinyExamples", "plotting", package = "benchmarkmeData")
  if(nchar(appDir) == 0) {
    stop("Could not find example directory. Try reinstalling `benchmarkmeData`.", call. = FALSE)
  }
  shiny::runApp(appDir)
}
