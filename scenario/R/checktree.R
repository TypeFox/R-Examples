#' @title Check the tree structure implied by a scenario tree nodal partition matrix.
#' @description Returns a plot showing the nodal structure (not values) of a scenario tree defined using a scenario tree nodal partition matrix.
#' @param treeStruct      Matrix defining the nodal structure of the tree.
#' @return Returns a plot of the scenario tree structure implied by the input nodal partition matrix.
#' @references Dupacova, Jitka, Giorgio Consigli, and Stein W. Wallace. "Scenarios for multistage stochastic programs." Annals of operations research 100.1-4 (2000): 25-53.
#' @examples
#' treeStruct <- rbind(c(1, 1, 1, 1, 1),
#'                     c(2, 2, 7, 7, 11),
#'                     c(3, 5, 8, 8, 12),
#'                     c(4, 6, 9, 10, 13)
#'                     )
#' checktree(treeStruct)
#' @importFrom graphics matplot axis
#' @export
checktree <- function(treeStruct){
  centroids <- matrix(NA, nrow = nrow(treeStruct), ncol = ncol(treeStruct))
  centroids[nrow(treeStruct),] <- seq(-1,1,length.out = ncol(treeStruct))
  for (t in (nrow(treeStruct) - 1):1){
    for (i in unique(treeStruct[t,])){
      x <- which(treeStruct == i, arr.ind = TRUE)
      x[,1] <- x[,1] + 1
      centroids[which(treeStruct == i)] <- mean(centroids[mean(x[,1]),][x[,2]])
    }
  }
  matplot(centroids, type = "b", lwd=2, pch = 15,
          xlab = "Time step", xaxt="n",
          main = "Scenario tree structure",
          yaxt="n", ylab = "")
  axis(1, at = 1:nrow(treeStruct), labels = 1:nrow(treeStruct))
}
