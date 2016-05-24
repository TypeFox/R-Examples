#' @title Build a scenario tree with a predefined nodal structure.
#' @description Uses the neural gas method to build a scenario tree.
#' @param x               Matrix of initial scenarios, realizations or ensemble members. Each column stores a scenario, with number of rows equal to number of time steps.
#' @param treeStruct      Matrix defining the nodal structure of the tree (see example). This is a scenario tree nodal partition matrix.
#' @param lambda_0        Upper neighbourhood range parameter. Default = 10.
#' @param lambda_f        Lower neighborhood range paramger. Default = 0.01.
#' @param e_0             Upper adaptation step parameter. Default = 0.5.
#' @param e_f             Lower adaptation step parameter. Default = 0.05.
#' @param jMax            Number of iterations. Default = 40000.
#' @param plot            logical. If TRUE (the default) the final tree is plotted.
#' @return Returns a list object containing the initial input scenarios, the input scenarios tree structure, the values of the final reduced scenario tree, and the tree branch probabilities at the end nodes.
#' @references Xu, B., Zhong, P.A., Zambon, R.C., Zhao, Y., Yeh, W. (2015) Scenario tree reduction in stochastic programming with recourse for hydropower operations, Water Resources Research, 51, 6359-6380.
#' @references Dupacova, Jitka, Giorgio Consigli, and Stein W. Wallace. "Scenarios for multistage stochastic programs." Annals of operations research 100.1-4 (2000): 25-53.
#' @examples # Generate some 25 random realizations of length 4 and reduce to scenario tree.
#' scenarios <- matrix(rnorm(100),ncol=25)
#' treeStruct <- rbind(c(1, 1, 1, 1, 1),
#'                     c(2, 2, 7, 7, 11),
#'                     c(3, 5, 8, 8, 12),
#'                     c(4, 6, 9, 10, 13)
#'                     )
#' tree <- buildtree(scenarios, treeStruct, jMax = 1000)
#' @importFrom graphics matplot matlines
#' @export
buildtree <- function(x, treeStruct, lambda_0 = 10, lambda_f = 0.01,
                      e_0 = 0.5, e_f = 0.05, jMax = 40000, plot = TRUE){


  ## -------- SET UP FUNCTIONS --------##
  getEucDist <- function(scenario, member){
    EucDist <- sqrt(sum((scenario - member)^2))
    return(EucDist)
  }
  step_size <- function(e_0, e_f, jMax, j){
    e_0 * ( (e_f / e_0) ^ (j / jMax))
  }
  getLambda_j <- function(lambda_f, lambda_0, jMax, j){
    lambda_j <- lambda_0 * ( (lambda_f/lambda_0) ^ (j / jMax))
    return(lambda_j)
  }
  getH <- function(O, lambda_j){
    h <- exp(-O/lambda_j)
    return(h)
  }

  ## -------- INITIALIZE TREE NODES -------- ##
  numScenarios <- ncol(treeStruct)
  tree <- x[,sample(ncol(x), numScenarios)]    # select x at random
  for (i in 1:max(treeStruct)){ #max(treeStruct) gives the number of nodes in the tree
    tree[which(treeStruct == i)] <- mean(tree[which(treeStruct == i)])
  }

  ## -------- ITERATE NODES -------- ##
  for (j in 1:jMax){
    randMember <- x[,sample(ncol(x), 1)] # Sample a single member from the x
    Distances <- apply(tree, 2, getEucDist, member = randMember)
    O <- rank(Distances)
    for (node in 1:max(treeStruct)){
      adapt <- getH(O[which(treeStruct==node, arr.ind=TRUE)[,2]],getLambda_j(lambda_f, lambda_0, jMax, j))
      diff <- randMember[which(treeStruct==node, arr.ind = TRUE)[,1]] - tree[which(treeStruct==node)][1]
      delta_node <- step_size(e_0, e_f, jMax, j) * mean(adapt * diff)
      tree[which(treeStruct==node)] <- tree[which(treeStruct==node)] + delta_node
    }
  }


  ## -------- CALCULATE BRANCH PROBABILITIES -------- ##
  branchProbs <- vector("integer", ncol(tree))
  for (m in 1:ncol(x)){
    dist <- apply(tree, 2, getEucDist, member = x[,m])
    branchProbs[which.min(dist)] <- branchProbs[which.min(dist)] + 1 / ncol(x)
  }

  ## -------- OUTPUT -------- ##
  if (plot) {
    matplot(x, type = "l", col = "grey", lty = 2, ylab = "Disturbance", xlab = "Time step")
    matlines(tree, pch = 3, lty = 1)
  }
  output <- list(treeStruct,tree, branchProbs)
  names(output) <- c("tree_structure", "tree_values", "branch_probabilities")
  return(output)
}
