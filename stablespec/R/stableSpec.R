#' Search stable specifications (structures) of constrained structural equation models.
#' @title Stable specifications of constrained structural equation models.
#' @param theData a data frame containing the data to which the model will be
#' be fit. If argument \code{longitudinal} is \code{TRUE}, the data frame
#' should be reshaped such that the first \code{n} data points contain the
#' relations that occur in the first two time slices \code{t_0} and \code{t_1}.
#' The next \code{n} data points contain the relations that occur in
#' time slices \code{t_1} and \code{t_2}. The \code{i-th} subset of \code{n}
#' data points contain the relations in time slices \code{t_i-1} and \code{t_i}.
#' @param nSubset number of subsets to draw. In practice, it is suggested
#' to have at least 25 subsets. The default is 10.
#' @param iteration number of iterations/generations for NSGA-II.
#' @param nPop population size (number of models) in a generation.
#' The default is 50.
#' @param mutRate mutation rate. The default is 0.075.
#' @param crossRate crossover rate. The default is 0.85.
#' @param longitudinal \code{TRUE} for longitudinal data,
#' and \code{FALSE} for cross-sectional data.
#' @param numTime number of time slices. If a cross-sectional data then
#' it is 1. The default is 1.
#' @param seed integer vector representing seeds that are used to subsample data.
#' The default is an integer vector with range \code{100:1000} with length
#' equal to \code{nSubset}.
#' @param co whether to use \code{"covariance"} or \code{"correlation"} matrix.
#' The default is \code{"covariance"}.
#' @param consMatrix \code{m by 2} binary \code{\link{matrix}}
#' representing constraint/prior knowledge,
#' where \code{m} is the number of constraint. For example, known that
#' variables 2 and 3 do not cause variable 1, then
#' \code{constraint <- matrix(c(2, 1, 3, 1), 2, 2, byrow=TRUE))} will be
#' the constraint matrix. If \code{NULL}, then it is assumed
#' that there is no constraint.
#' @param threshold threshold of stability selection. The default is 0.6.
#' @param toPlot if \code{TRUE} a plot of inferred causal model is generated,
#' otherwise a graph object is returned. The default is \code{TRUE}.
#' @return a list of the following elements:
#' \itemize{
#' \item \code{listofFronts} is a \code{\link{list}} of optimal models for
#' the whole range of model complexity of all subsets.
#' \item \code{causalStab} is a \code{\link{list}} of causal path stability
#' for the whole range of model complexity
#' \item \code{causalStab_l1} is a \code{\link{list}} of
#' causal path stability of length 1
#' for the whole range of model complexity
#' \item \code{edgeStab} is a \code{\link{list}} of edge stability
#' for the whole range of mdoel complexity
#' \item \code{relCausalPath} is \code{n by n} \code{\link{matrix}} of
#' relevant causal path,
#' where \code{n} is the number of variables. Each positive element
#' \code{i,j} represents the stability of causal path
#' from \code{i} to \code{j}.
#' \item \code{relCausalPath_l1} is \code{n by n} \code{\link{matrix}}
#' of relevant causal path with length 1, where \code{n} is
#' the number of variables. Each positive element \code{i,j}
#' represents the stability of causal path
#' from \code{i} to \code{j} with length 1.
#' \item \code{relEdge} is \code{n by n} \code{\link{matrix}} of relevant edge,
#' where \code{n} is the number of variables. Each positive element
#' \code{i,j} represents the stability of edge
#' between \code{i} to \code{j}.
#' \item If argument \code{toPlot = TRUE}, then a plot of inferred causal model
#' is generated. Otherwise an object of graph is returned.
#' An arc represents a causal path, and an (undirected)
#' edge represents strong association where the direction is undecidable.
#' \item \code{allSeed} is an integer vector representing seeds that are used in
#' subsampling data. This can be used to replicate the result
#' in next computation.
#' }
#'
#' @examples
#' # Cross-sectional data example.
#' # with a data set about patients with ADHD.
#' # Detail about the data set can be found in the documentation.
#' # As an example, we only run one subset.
#' the_data <- adhd
#' numSubset <- 1
#' num_iteration <- 5
#' num_pop <- 10
#' mut_rate <- 0.075
#' cross_rate <- 0.85
#' longi <- FALSE
#' num_time <- 1
#' the_seed <- NULL
#' the_co <- "covariance"
#' # assummed that nothing causing variable Gender
#' cons_matrix <- matrix(c(2, 1, 3, 1, 4, 1, 5, 1, 6, 1), 5, 2, byrow=TRUE)
#' th <- 0.1
#' to_plot <- FALSE
#'
#' result_adhd <- stableSpec(theData=the_data, nSubset=numSubset,
#' iteration=num_iteration,
#' nPop=num_pop, mutRate=mut_rate, crossRate=cross_rate,
#' longitudinal=longi, numTime=num_time, seed=the_seed,
#' co=the_co, consMatrix=cons_matrix, threshold=th, toPlot=to_plot)
#'
#' @author Ridho Rahmadi \email{r.rahmadi@cs.ru.nl}, Perry Groot, Tom Heskes
#' @details This function performs exploratory search over
#' recursive (acyclic) SEM models.
#' Models are scored along two objectives: the model fit and
#' the model complexity. Since both objectives are often conflicting
#' we use NSGA-II to search for Pareto optimal models. To handle the
#' instability of small finite data samples, we repeatedly subsample
#' the data and select those substructures that are both stable and
#' parsimonious which are then used to infer a causal model.
#' @references
#' Rahmadi, R., Groot, P., Heins, M., Knoop, H., & Heskes, T. (2015).
#' Causality on Cross-Sectional Data: Stable Specification Search in
#' Constrained Structural Equation Modeling. arXiv preprint arXiv:1506.05600.
#'
#' John Fox, Zhenghua Nie and Jarrett Byrnes (2015). sem:
#' Structural Equation Models. R package version 3.1-6.
#' https://CRAN.R-project.org/package=sem
#'
#' Ching-Shih Tsou (2013). nsga2R: Elitist Non-dominated Sorting
#' Genetic Algorithm based on R. R package version 1.0.
#' https://CRAN.R-project.org/package=nsga2R
#'
#' Kalisch, M., Machler, M., Colombo, D., Maathuis, M. H., &
#' Buehlmann, P. (2012). Causal inference using graphical models
#' with the R package pcalg.
#' \emph{Journal of Statistical Software}, 47(11), 1-26.
#'
#' Meinshausen, N., & Buehlmann, P. (2010). Stability selection.
#' \emph{Journal of the Royal Statistical Society:
#' Series B (Statistical Methodology)}, 72(4), 417-473.
#'
#' Deb, K., Pratap, A., Agarwal, S., and Meyarivan, T. (2002),
#' A fast and elitist multiobjective genetic algorithm: NSGA-II,
#' \emph{IEEE Transactions on Evolutionary Computation}, 6(2), 182-197.
#'
#' Chickering, D. M. (2002). Learning equivalence classes of
#' Bayesian-network structures. \emph{The Journal of
#' Machine Learning Research}, 2, 445-498.
#'
#' @importClassesFrom graph graphNEL
#' @importFrom graphics axis lines par plot title
#' @importFrom methods as
#' @importFrom stats cor cov runif
#' @importFrom utils head tail
#' @export
stableSpec <- function(theData = NULL,
                       nSubset = NULL,
                       iteration = NULL,
                       nPop = NULL,
                       mutRate = NULL,
                       crossRate = NULL,
                       longitudinal = NULL,
                       numTime = NULL,
                       seed = NULL,
                       co = NULL,
                       consMatrix = NULL,
                       threshold = NULL,
                       toPlot = NULL) {

  # to check arguments
  # argument data
  if(!is.null(theData)) { # if data is supplied
    if (is.numeric(theData) && !(is.matrix(theData))) {
      stop("Data should be either a data frame or a matrix of numerical values.")
    } else if (!(is.numeric(theData)) && is.data.frame(theData)) {
      if (any(sapply(theData, is.numeric) == FALSE)) {
        stop("Data should be either a data frame or a matrix of numerical values.")
      }
    } else if (!is.numeric(theData)) {
      stop("Data should be either a data frame or a matrix of numerical values.")
    }
  } else { # if not supplied
    stop("Data cannot be missing")
  }


  # arguments nSubset, iteration, nPop, mutRate, crossRate, iteration
  if (!is.null(nSubset)) {
    if (!is.numeric(nSubset) || is.matrix(nSubset)) {
      stop("Argument nSubset should be positive numeric, e.g., 10.")
    }
  } else {
    nSubset <- 10
  }

  if (!is.null(nPop)) {
    if (!(is.numeric(nPop)) || is.matrix(nPop)) {
      stop("Argument nPop should be positive numeric, e.g., 50.")
    }
  } else {
    nPop <- 50
  }

  if (!is.null(iteration)) {
    if (!is.numeric(iteration) || is.matrix(iteration)) {
      stop("Argument iteration or NSGA-II generations should be positive numeric, e.g., 20.")
    }
  } else {
    iteration <- 20
  }

  if (!is.null(numTime)) {
    if (!is.numeric(numTime) || is.matrix(numTime)) {
      stop("Argument numTime should be positive numeric, e.g., 1 if cross-sectional data.")
    }
  } else {
    numTime <- 1
  }

  if (!is.null(mutRate)) {
    if (!is.numeric(mutRate) || is.matrix(mutRate)) {
      stop("Argument mutRate should be positive numeric, e.g., 0.075.")
    }
  } else {
    mutRate <- 0.075
  }

  if (!is.null(crossRate)) {
    if (!is.numeric(crossRate) || is.matrix(crossRate)) {
      stop("Argument crossRate should be positive numeric, e.g., 0.85.")
    }
  } else {
    crossRate <- 0.85
  }

  if (!is.null(threshold)) {
    if (!is.numeric(threshold) || is.matrix(threshold)) {
      stop("Argument threshold should be positive numeric, e.g., 0.6.")
    }
  } else {
    threshold <- 0.6
  }

  # arguments nTime and longitudinal
  if (!is.null(longitudinal)) {
    if (!is.logical(longitudinal)) {
      stop("Argument longitudinal should be either logical TRUE or FALSE.")
    }
  } else {
    stop("Argument longitudinal cannot be missing.")
  }

  if (!longitudinal && numTime > 1) {
    stop("Cross-sectional data should have only one time slice, e.g., numTime = 1")
  } else if (longitudinal && numTime == 1) {
    stop("Longitudinal data should have more than one time slices, e.g., numTime = 2, with two time slices.")
  }

  # argument consMatrix
  if(!is.null(consMatrix)) {
    if(!is.matrix(consMatrix)) {
      stop("The constraints should be formed in a matrix.")
    }
  }

  # argument co
  if(!is.null(co)) {
    if(!is.character(co)) {
      stop("Argument co should be a vector of characters, e.g., either covariance or correlation.")
    } else {
      covMatrix <- c("covariance", "correlation")
      if (!co %in% covMatrix) {
        stop("Argument co should be either covariance or correlation matrix.")
      }
    }
  } else {
    co <- "covariance"
  }

  # argument toPlot
  if (!is.null(toPlot)) {
    if (!is.logical(toPlot)) {
      stop("Argument toPlot should be either logical TRUE or FALSE.")
    }
  } else {
    toPlot <- TRUE
  }


  if (!is.null(seed)) {
    if (!is.numeric(seed) || is.matrix(seed)) {
      stop("Argument seed should be numeric vector.")
    }
  } else {
    seed <- sample(100:1000, nSubset)
  }

  #get the optimal models from the whole range of model complexities
  optimal_models <- optimalModels(theData, nSubset, iteration, nPop,
                                  mutRate, crossRate, longitudinal,
                                  numTime, seed, co, consMatrix)


  #doIt <- TRUE
  #while(doIt) {
    #print("Now the relevant specifications are computed....")

    #compute the stability of structures
    stabRes <- structureStab(optimal_models$listOfFronts,
                             optimal_models$string_size,
                             optimal_models$num_var,
                             longitudinal,
                             consMatrix)


    #relevant structures
    rel_struct <-
      relevantStructure(optimal_models$listOfFronts,
                        threshold, stabRes$causalStab,
                        stabRes$causalStab_l1, stabRes$edgeStab,
                        optimal_models$string_size)

    #doIt <- FALSE
  #}

  #return output
  if (toPlot) {
    return(list(listOfFronts=optimal_models$listOfFronts,
                causalStab=stabRes$causalStab,
                causalStab_l1=stabRes$causalStab_l1,
                edgeStab=stabRes$edgeStab,
                relCausalPath=rel_struct$relCausalPath,
                relCausalPath_l1=rel_struct$relCausalPathL1,
                relEdge=rel_struct$relEdge,
                graph=Rgraphviz::renderGraph(rel_struct$graph),
                allSeed=seed))
  } else {
    return(list(listOfFronts=optimal_models$listOfFronts,
                causalStab=stabRes$causalStab,
                causalStab_l1=stabRes$causalStab_l1,
                edgeStab=stabRes$edgeStab,
                relCausalPath=rel_struct$relCausalPath,
                relCausalPath_l1=rel_struct$relCausalPathL1,
                relEdge=rel_struct$relEdge,
                graph=rel_struct$graph,
                allSeed=seed))

  }
}
