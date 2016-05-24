#' @include proto.R
{}



#' Benchmark experiment warehouse
#'
#' \code{warehouse} is the constructor function for a benchmark experiment
#' warehouse.
#'
#' @param datasets Names of the datasets
#' @param B Number of benchmark runs
#' @param algorithms Names of the candidate algorithms
#' @param performances Names of the performance measures
#' @param characteristics Names of the dataset characteristics
#' @param tests Names of the monitored test measures
#' @return Proto object with different views (see Details).
#'
#' @details
#'   A benchmark experiment warehouse collects all data during the
#'   execution of a benchmark experiment; see \code{\link{benchmark}}.
#'   Different views (based on the collected data) provide cleaned
#'   parts of the data for further analyses.
#'
#'   Implemented views:
#'   \enumerate{
#'   \item \code{viewAlgorithmPerformance()}: returns a data frame (S3
#'   class \code{AlgorithmPerformance}) with columns \code{samples,
#'   datasets, algorithms, performances} (factors with the corresponding
#'   levels) and the column \code{value} with the corresponding
#'   performance value.
#'
#'   \item \code{viewDatasetCharacterization()}: returns a data frame
#'   (S3 class \code{DatasetCharacterization}) with columns \code{samples,
#'   datasets, characteristics, value}.
#'
#'   \item \code{viewDatasetBasisCharacterization()}: returns a data
#'   frame (S3 class \code{DatasetBasisCharacterization}) with columns
#'   \code{datasets, characteristics, value}.
#'
#'   \item \code{viewTestResult()}: returns a data frame (S3 class
#'   \code{TestResult}) with columns \code{samples, datasets, tests, value}.
#'   }
#'
#' @seealso \code{\link{benchmark}}, \code{\link{as.warehouse}}
#'
#' @aliases AlgorithmPerformance DatasetCharacterization
#'   DatasetBasisCharacterization TestResult
#'
#' @import proto
#' @importFrom reshape melt
#' @export
warehouse <- function(datasets, B,
                      algorithms = NULL,
                      performances = NULL,
                      characteristics = NULL,
                      tests = NULL) {

  if ( length(datasets) != length(B) )
    B <- rep(B, length(datasets))


  a <- mapply(DatasetList, datasets, B,
              MoreArgs = list(algorithms = algorithms,
                              performances = performances,
                              characteristics = characteristics,
                              tests = tests),
              SIMPLIFY = FALSE)
  names(a) <- datasets


  ## Proto object and default data views:
  a <- proto(data = a)

  a$meta <- list(datasets = datasets, B = B,
                 algorithms = algorithms,
                 performances = performances,
                 characteristics = characteristics,
                 tests = tests,
                 algorithm_colors = default_colors(algorithms = algorithms))


  if ( !is.null(algorithms) & !is.null(performances) ) {
    setViewAlgorithmPerformance(a)

    if ( !is.null(tests) )
      setViewTestResult(a)
  }

  if ( !is.null(characteristics) ) {
    setViewDatasetCharacterization(a)
    setViewDatasetBasisCharacterization(a)
  }


  structure(a, class = c("warehouse", class(a)))
}



#' @S3method print warehouse
print.warehouse <- function(x, ...) {
  out <- matrix(NA_character_, nrow = 2, ncol = 2)
  out[1, ] <- c("Algorithms", paste(x$meta$algorithms, collapse = ", "))
  out[2, ] <- c("Datasets", paste(x$meta$datasets, collapse = ", "))

  out2 <- x$meta[4:6]
  out2 <- !sapply(out2, is.null)

  cat("Warehouse object\n")
  writeLines(formatDL(out, style = "list"))
  cat("Data:\n")
  print(out2)

  invisible(x)
}



### Data views: ######################################################

setViewAlgorithmPerformance <- function(object) {

  object$viewAlgorithmPerformance <- function(.,
                                              datasets = NULL,
                                              algorithms = NULL,
                                              performances = NULL) {

    if ( is.null(datasets) )
      datasets <- .$meta$datasets

    if ( is.null(algorithms) )
      algorithms <- .$meta$algorithms

    if ( is.null(performances) )
      performances <- .$meta$performances


    view <- lapply(.$data[datasets],
                   function(ds)
                   ds$AlgorithmPerformance[,
                                           algorithms,
                                           performances,
                                           drop = FALSE])
    attr(view, "varname") <- "datasets"

    view <- melt(view)
    view$datasets <- as.factor(view$datasets)
    view$samples <- as.factor(view$samples)

    view <- view[, c("samples", "datasets", "algorithms", "performances", "value")]


    structure(view, class = c("AlgorithmPerformance", class(view)),
              algorithm_colors = .$meta$algorithm_colors)
  }

  invisible(NULL)
}



setViewDatasetCharacterization <- function(object) {

  object$viewDatasetCharacterization <- function(.,
                                                 datasets = NULL,
                                                 characteristics = NULL,
                                                 basis = TRUE) {

    if ( is.null(datasets) )
      datasets <- .$meta$datasets

    if ( is.null(characteristics) )
      characteristics <- .$meta$characteristics


    view <- lapply(.$data[datasets],
                   function(ds)
                   ds$DatasetCharacterization[,
                                              characteristics,
                                              drop = FALSE])
    attr(view, "varname") <- "datasets"

    view <- melt(view)
    view$datasets <- as.factor(view$datasets)
    view$samples <- as.factor(view$samples)

    view <- view[, c("samples", "datasets", "characteristics", "value")]

    if ( basis ) {
        basis <- .$viewDatasetBasisCharacterization(datasets = datasets,
                                                    characteristics = characteristics)
        basis$samples <- "basis"
        view <- rbind(view, basis)
    }


    structure(view, class = c("DatasetCharacterization", class(view)))
  }

  invisible(NULL)
}



setViewDatasetBasisCharacterization <- function(object) {

  object$viewDatasetBasisCharacterization <- function(.,
                                                      datasets = NULL,
                                                      characteristics = NULL) {

    if ( is.null(datasets) )
      datasets <- .$meta$datasets

    if ( is.null(characteristics) )
      characteristics <- .$meta$characteristics


    view <- lapply(.$data[datasets],
                   function(ds)
                   ds$DatasetBasisCharacterization[,
                                                   characteristics,
                                                   drop = FALSE])
    attr(view, "varname") <- "datasets"

    view <- melt(view)
    view$datasets <- as.factor(view$datasets)

    view <- view[, c("datasets", "characteristics", "value")]


    structure(view, class = c("DatasetBasisCharacterization", class(view)))
  }

  invisible(NULL)
}



setViewTestResult <- function(object) {

  object$viewTestResult <- function(.,
                                    datasets = NULL,
                                    tests = NULL) {

    if ( is.null(datasets) )
      datasets <- .$meta$datasets

    if ( is.null(tests) )
      tests <- .$meta$tests


    view <- lapply(.$data[datasets],
                   function(ds)
                   ds$TestResult[,
                                 tests,
                                 drop = FALSE])
    attr(view, "varname") <- "datasets"

    view <- melt(view)
    view$datasets <- as.factor(view$datasets)

    view <- view[, c("samples", "datasets", "tests", "value")]


    structure(view, class = c("TestResult", class(view)))
  }

  invisible(NULL)
}



### Internal data structures: ########################################

WarehouseArray <- function(B, ..., class) {
  d <- list(...)

  dim <- c(B, sapply(d, length))
  dimnames <- c(list(samples = NULL), d)

  a <- array(NA_integer_, dim = dim, dimnames = dimnames)

  structure(a, class = c(class, class(a)))
}



AlgorithmPerformanceArray <- function(B, algorithms, performances) {
  WarehouseArray(B, algorithms = algorithms, performances = performances,
                 class = "AlgorithmPerformanceArray")
}



DatasetCharacterizationArray <- function(B, characteristics) {
  WarehouseArray(B, characteristics = characteristics,
                 class = "DatasetCharacterizationArray")
}



TestResultArray <- function(B, tests) {
  WarehouseArray(B, tests = tests,
                 class = "TestResultArray")
}



DatasetList <- function(dataset, B,
                        algorithms = NULL,
                        performances = NULL,
                        characteristics = NULL,
                        tests = NULL) {

  a <- list()

  if ( !is.null(algorithms) && !is.null(performances) ) {
    a$AlgorithmPerformance <- AlgorithmPerformanceArray(B, algorithms,
                                                        performances)
  }

  if ( !is.null(characteristics) ) {
    a$DatasetCharacterization <- DatasetCharacterizationArray(B, characteristics)
    a$DatasetBasisCharacterization <- DatasetCharacterizationArray(1, characteristics)
  }

  if ( !is.null(tests) ) {
    a$TestResult <- TestResultArray(B, tests)
  }


  structure(a, class = c("DatasetList", class(a)),
            dataset = dataset)
}



### Internal functions: ##############################################

default_colors <- function(n = length(algorithms), algorithms = NULL) {
  # Based on ggplot2:::ScaleHue
  h <- c(0, 360) + 15
  l <- 65
  c <- 100

  start <- 1
  direction <- -1

  rotate <- function(x) (x + start) %% 360 * direction

  if ( (diff(h) %% 360) < 1 ) {
    h[2] <- h[2] - 360 / n
  }

  structure(grDevices::hcl(h = rotate(seq(h[1], h[2], length = n)),
                           c = c, l = l),
            names = algorithms)
}


