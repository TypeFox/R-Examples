#install.packages("doParallel", repos="http://cran.rstudio.com/") # Library to parallelize the execution of the program.

#library(doParallel)

#' Parallel computation of the lower/upper bounds of the uncertainty interval of the Rao-Stirling diversity index.
#'
#' This function allows the computation of the lower/upper bounds of the uncertainty interval of the Rao-Stirling index (Calatrava et al., 2016) in parallel threads.
#' It includes the parallel computation of the permissible disciplines (i.e., function \code{\link{PruneDisciplines}}).
#' The use of this function is recommended for an efficient computation of the lower and upper bounds of the uncertainty interval of the Rao-Stirling index.
#' The computation of the lower bound is an NP-hard problem.
#' Because the computation of the lower bound might require long computing times, this function creates a log file 'parallel-bounds-log.txt' in the user's workspace.
#' The content of the log file is written during the execution of the function and indicates number of publications that have been processed.
#'
#' @param bound String that indicates which index to compute.
#' Two values are valid: \emph{upper} and \emph{lower}.
#' @param count.matrix Vector or matrix that contains the counts of references to different disciplines of a single publication (a vector) or of several publications (a matrix).
#' If it is a vector its length is equal to the total number of disciplines.
#' In case it is a matrix its dimensions are \emph{n} x \emph{m}, being n the total number of disciplines and m the number of publications for which the lower/upper bound will be computed.
#' @param uncat.refs Number of uncategorized references of a publication (a number) or several publications (a vector).
#' @param similarity A positive semi-definite matrix that encodes the similarity between disciplines, as explain in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The self-similarities (i.e. values in the diagonal) have to be 1.
#' @param pruning Logical value that indicates whether the set of permissible disciplines will be calculated and used to avoid improbable redistributions of disciplines.
#' This argument is optional and leaving it unspecified ignores the pruning of unlikely disciplines in the redistribution.
#' @param tolerance A real number in the interval [0,1].
#' This argument modulates the similarity between disciplines with which the strictness of the pruning of unlikely disciplines is controlled.
#' A value of 0 allows all disciplines to participate in the redistribution process.
#' A value of 1 permits no tolerance.
#' This argument is optional and leaving it unspecified deactivates tolerances.
#' @param redistribution.limit A positive integer that limits the number of disciplines that each uncategorized reference can have redistributed.
#' This argument is optional and leaving it unspecified will set the redistribution.limit to default.
#' @param threads A positive number that specifies the number of parallel threads that will be executed.
#' This argument should be set according to the number of processor core in the CPU of the user.
#' This argument is optional and leaving it unspecified will set the number of threads to default.
#' @param max.batch.size A positive integer that sets the size of the batch of candidates that is computed at once.
#' This positive value determines the quantity of allocated memory and has to be reduced if corresponding errors arise.
#' This argument is optional and leaving it unspecified sets it to a default value.
#' @return The lower or the upper bound/s of the uncertainty interval of the Rao-Stirling index of one publication (an integer) or several publications (a vector).
#' @section Warning:
#' This function solves a computationally intensive optimization problem. In order to reduce the search space it is recommended to provide the function with the vector of permissible disciplines and redistribution limit.
#' @examples
#' #Load data
#' data(pubdata1)
#'
#' #Get upper bound indices of the uncertainty interval of the Rao-Stirling diversity index.
#' ParallelBoundIndices("upper", pd1.count.matrix, pd1.uncat.refs, pd1.similarity, TRUE, 0.233, 4, 2)
#'
#' #Get lower bound indices of the uncertainty interval of the Rao-Stirling diversity index.
#' ParallelBoundIndices("lower", pd1.count.matrix, pd1.uncat.refs, pd1.similarity, TRUE, 0.233, 4, 2)
#'
#' #When many references of a publication are uncategorized, a warning message is displayed
#' #to inform the user. Such cases require longer computation times.

#'
#' @references
#' Calatrava Moreno, M. C., Auzinger, T. and Werthner, H. (2016) On the uncertainty of interdisciplinarity measurements due to incomplete bibliographic data. Scientometrics. DOI:10.1007/s11192-016-1842-4
#'
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @import doParallel foreach parallel
#' @export
ParallelBoundIndices <- function(bound,
                                 count.matrix,
                                 uncat.refs,
                                 similarity,
                                 pruning = TRUE,
                                 tolerance = 1,
                                 redistribution.limit = 4,
                                 threads = 1,
                                 max.batch.size = 131072) {
  # Error handling.
  if (!((bound == "lower") || (bound == "upper"))) {
    stop("Invalid argument value 'bound': only 'lower' and 'upper' are valid values.")
  }
  if (is.vector(count.matrix)) {
    if ((is.numeric(uncat.refs) && (length(uncat.refs)!=1))) {
      stop("Dimensions of 'count.matrix' and 'uncat.refs' do not match.")
    }
  } else if (is.matrix(count.matrix)) {
    if (!is.vector(uncat.refs)) {
      stop ("Invalid dimensions argument 'uncat.refs'.")
    }
    if (length(uncat.refs) != dim(count.matrix)[2]) {
      stop("Dimensions of 'count.matrix' and 'uncat.refs' do not match.")
    }
  } else {
    stop("Invalid argument type 'count.matrix'.")
  }
  if (!((pruning == TRUE) || (pruning == FALSE))) {
    stop("Invalid argument value 'pruning': only logic values are valid.")
  }
  if (is.nan(tolerance) || tolerance < 0 || tolerance > 1) {
    stop("Argument 'tolerance' is out of range.")
  }
  if (is.nan(redistribution.limit) || (redistribution.limit < 0)) {
    stop("Argument 'redistribution.limit' is out of range.")
  }
  if (is.nan(threads) || threads < 0) {
    stop("Argument 'threads' is out of range.")
  }

  # number of cores to be used
  cl <- makeCluster(threads)
  # make a cluster of cores available
  registerDoParallel(cl)

  bound.indices <- c()

  # Initialization of log file
  cat(paste("Log thread processing of ",bound, "index","\n"), file="parallel-bounds-log.txt")

  # Obtain the number of papers for which the index will be calculated
  if (is.vector(count.matrix)) {
    num.papers <- 1
  } else if (is.matrix(count.matrix)) {
    num.papers <- dim(count.matrix)[2]
  }

  i <- 0
  bound.indices <- foreach (i=1:num.papers,
                            .export=ls(envir=globalenv()),
                            .combine=c,
                            #.noexport=c("count.matrix", "similarity", "uncat.refs"),
                            #.verbose = TRUE,
                            .packages=c("igraph","polynom","methods","gmp","quadprog","digest","iterpc")
  ) %dopar% {

    cat(paste("Starting iteration ",i,"\n"), file="parallel-bounds-log.txt", append=TRUE)

    # Obtain the count vector and number of uncategorized references for each paper
    if (is.vector(count.matrix)) {
      known.ref.counts <- count.matrix
      uncat.ref.count <- uncat.refs
    } else if (is.matrix(count.matrix)) {
      known.ref.counts <- as.vector(count.matrix[,i])
      uncat.ref.count <- uncat.refs[i]
    }

    if ((all(known.ref.counts == 0)) && (uncat.ref.count == 0)) {
      bound.index <- NA
    } else {
      # Obtain a logic vector of known referenced disciplines
      logic.known.ref.counts <- known.ref.counts > 0

      # Select if pruning will be applied
      if (pruning) {
        permissible.disciplines <- PruneDisciplines(logic.known.ref.counts, tolerance, similarity)
      } else {
        permissible.disciplines <- NULL
      }

      if (bound == "lower") {
        bound.index <- LowerIndexBound(known.ref.counts, uncat.ref.count, similarity, permissible.disciplines, redistribution.limit, max.batch.size)
      } else if (bound == "upper") {
        bound.index <- UpperIndexBound(known.ref.counts, uncat.ref.count, similarity, permissible.disciplines, redistribution.limit)
      }
    }
  }

  # free the cores
  stopCluster(cl)

  return (bound.indices)
}


#' Rao-Stirling diversity index based on the counts of cited disciplines.
#'
#' This function calculates the Rao-Stirling diversity index of one or several publications, based on the count of citations of the publication(s) to different disciplines.
#'
#' @param count.matrix Vector or matrix that contains the counts of references to different disciplines of a single publication (vector) or of several publications (matrix).
#' If count.matrix is a vector its length is equal to the total number of disciplines.
#' In case it is a matrix its dimensions are \emph{n} x \emph{m}, being n the total number of disciplines and m the number of publications for which the lower/upper bound will be computed.
#' @param similarity A positive semi-definite matrix that encodes the similarity between disciplines, as explain in Porter and Rafols (2009).
#' The dimensions of this matrix are \emph{n} x \emph{n}, being \emph{n} the total number of disciplines.
#' The number of rows and the number columns of this matrix need to be equal to the number of rows of \code{count.matrix}.
#' The self-similarities (i.e. the diagonal elements) have to be 1.
#' @return The Rao-Stirling diversity index of one or several publications.
#' @examples
#' #Load data
#' data(pubdata1)
#'
#' #Get Rao-Stirling diversity index of all publications in the dataset
#' RaoStirling(pd1.count.matrix, pd1.similarity)
#'
#' #Get Rao-Stirling diversity index of one publication of the dataset
#' RaoStirling(pd1.count.matrix[,2], pd1.similarity)
#'
#' @references
#' Porter, A. and Rafols, I. (2009) Is science becoming more interdisciplinary? Measuring and mapping six research fields over time. Scientometrics, Vol. 81, No. 3 (719-745). DOI:10.1007/s11192-008-2197-2
#' @export
RaoStirling <- function (count.matrix, similarity) {

  # Obtain the number of papers for which the index will be calculated
  if (is.vector(count.matrix)) {
    num.papers <- 1
  } else if (is.matrix(count.matrix)) {
    num.papers <- dim(count.matrix)[2]
  }

  rao.stirling <- c()

  for (i in 1:num.papers) {
    # Call the function that computes the Rao-Stirling index for a single paper.
    if (is.vector(count.matrix)) {
      if (sum(count.matrix) == 0) {
        rao.stirling <- c(rao.stirling, NA)
      } else {
        rao.stirling <- c(rao.stirling, RaoCounts (count.matrix, similarity))
      }
    # Call the function that computes the Rao-Stirling index for a set of papers.
    } else if (is.matrix(count.matrix)) {
      if (sum(count.matrix[,i]) == 0) {
        rao.stirling <- c(rao.stirling, NA)
      } else {
        rao.stirling <- c(rao.stirling, RaoCounts (count.matrix[,i], similarity))
      }
    }
  }
  return (rao.stirling)
}

