#' scenario: Construct reduced trees with a predefined nodal structures
#'

#' The \code{\link{buildtree}} function uses the neural gas method to generate a scenario tree of predefined nodal structure. The \code{\link{checktree}} function plots a scenario tree structure as defined by a nodal parition matrix.
#' @docType package
#' @name scenario
#' @references Latorre, J.M., Cerisola, S. and Ramos, A. (2007) Clustering algorithms for scenario tree generation: Application to natural hydro flows, European Journal of Operational Research, 181, 1339-1353.
#' @references Xu, B., Zhong, P.A., Zambon, R.C., Zhao, Y., Yeh, W. (2015) Scenario tree reduction in stochastic programming with recourse for hydropower operations, Water Resources Research, 51, 6359-6380.
#' @references Dupacova, Jitka, Giorgio Consigli, and Stein W. Wallace. "Scenarios for multistage stochastic programs." Annals of operations research 100.1-4 (2000): 25-53.
#' @examples # TEST BY GENERATING SCENARIOS FROM KNOWN CENTROIDS AND THEN
#' # COMPARING FIT BETWEEN THE GENERATED TREE AND INTIAL CENTROIDS.
#'
#' # 1. Generate scenarios with known centroids:
#'
#' centroids <- cbind(c(0,2,3), c(0,2,1), c(0,-2,-3),c(0,-2,-1))
#' matplot(centroids, type="l", lwd = 3, col = "black", lty = 3)
#' scenarios <- matrix(rep(centroids,5), ncol=20) + matrix(rnorm(60,0,0.25),ncol=20)
#' matlines(scenarios, col = "grey")
#'
#'
#' # 2. Assign and check nodal structure for tree:
#'
#' treeStruct <- rbind(c(1,1,1,1),
#'                     c(2,2,5,5),
#'                     c(3,4,6,7))
#' checktree(treeStruct)
#'
#'
#' # 3. Build scenario tree:
#'
#' tree <- buildtree(scenarios, treeStruct, jMax = 1000)
#'
#'
#' #4. Compare original centroids
#'
#' matlines(centroids,lwd = 3, col = "black", lty = 3)
#' # Improved convergence is achieved by increasing the number of iterations, jMax.

NULL
