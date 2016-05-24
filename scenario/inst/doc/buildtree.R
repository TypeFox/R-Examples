## ---- echo=FALSE, fig.show='hold', fig.width=7, fig.height=5-------------
centroids <- cbind(c(0,2,3), c(0,2,1), c(0,-2,-3),c(0,-2,-1))
matplot(centroids, type = "b", lwd=2, pch = 15,
        xlab = "Time step", xaxt="n",
        xlim=c(0.8,3.2), main = "Simple scenario tree structure",
        yaxt="n", ylab = "", ylim = c(-3.3,3.3))
axis(1, at=1:3, labels=1:3)
text(3.2,3,expression('s'["1"]*''), cex = 1.5)
text(3.2,1,expression('s'["2"]*''), cex = 1.5)
text(3.2,-1,expression('s'["3"]*''), cex = 1.5)
text(3.2,-3,expression('s'["4"]*''), cex = 1.5)
text(0.85,0.6,expression('s'["1,1"]*''))
text(0.85,0.2,expression('s'["2,1"]*''))
text(0.85,-0.2,expression('s'["3,1"]*''))
text(0.85,-0.6,expression('s'["4,1"]*''))
text(2,1.7,expression('s'["1,2"]*''))
text(2,1.3,expression('s'["2,2"]*''))
text(2,-1.3,expression('s'["3,2"]*''))
text(2,-1.7,expression('s'["4,2"]*''))
text(3,2.7,expression('s'["1,3"]*''))
text(3,1.3,expression('s'["2,3"]*''))
text(3,-1.3,expression('s'["3,3"]*''))
text(3,-2.7,expression('s'["4,3"]*''))

## ---- echo=FALSE, fig.show='hold', fig.width=7, fig.height=5-------------
matrix(c("S_1,1", "S_1,2", "S_1,3", "S_2,1", "S_2,2", "S_2,3", "S_3,1", "S_3,2", "S_3,3", "S_4,1", "S_4,2", "S_4,3" ), ncol=4)

## ---- echo=FALSE, fig.show='hold', fig.width=7, fig.height=5-------------
rbind(c(1,1,1,1),
      c(2,2,5,5),
      c(3,4,6,7))

## ---- fig.show='hold', fig.width=7, fig.height=5-------------------------
treeStruct <- rbind(c(1,1,1,1),
                    c(2,2,5,5),
                    c(3,4,6,7))
scenario::checktree(treeStruct)

## ---- fig.show='hold', fig.width=7, fig.height=5-------------------------
treeStruct <- rbind(c(1,1,1,1,1,1),
                    c(2,2,5,5,8,11),
                    c(3,4,6,7,9,12))
scenario::checktree(treeStruct)

## ---- fig.show='hold', fig.width=7, fig.height=5-------------------------
known_tree <- cbind(c(0,2,3),
                   c(0,2,1),
                   c(0,-2,-3),
                   c(0,-2,-1)
                   )
# now add some noise to the known tree...
realizations <- matrix(rep(known_tree,5), ncol=20) + matrix(rnorm(60,0,1),ncol=20)
matplot(realizations, lty=2, col = "grey", type="l",
        ylab = "Disturbance", xaxt = "n", main = "Initial realizations")
axis(1, at=1:3, labels=1:3)

## ---- echo=FALSE, fig.show='hold', fig.width=7, fig.height=5-------------
treeStruct <- rbind(c(1,1,1,1),
                    c(2,2,5,5),
                    c(3,4,6,7))
numScenarios <- ncol(treeStruct)
tree <- realizations[,sample(ncol(realizations), numScenarios)]
for (i in 1:max(treeStruct)){
  tree[which(treeStruct == i)] <- mean(tree[which(treeStruct == i)])
}
matplot(realizations, lty=2, col = "grey", type="l",
        ylab = "Disturbance", xaxt = "n", main = "Initial tree node positions")
axis(1, at=1:3, labels=1:3)
matlines(tree, pch = 3, lty = 1)

## ---- echo=FALSE, fig.show='hold', fig.width=7, fig.height=5-------------
matplot(realizations, lty=2, col = "lightgrey", type="l",
        ylab = "Disturbance", xaxt = "n", main = "Randomly selected realization")
axis(1, at=1:3, labels=1:3)
matlines(tree, pch = 3)
lines(realizations[,1], lwd = 3, lty = 2)

## ---- echo=FALSE, fig.show='hold', fig.width=7, fig.height=5-------------
getEucDist <- function(scenario, member){
  EucDist <- sqrt(sum((scenario - member)^2))
  return(EucDist)
}
Euclidean_distances <- apply(tree, 2, getEucDist, member = realizations[,1])
Rank <- rank(Euclidean_distances)
cbind(Euclidean_distances,Rank)

## ---- fig.show='hold', fig.width=7, fig.height=5-------------------------
scenario::buildtree(realizations, treeStruct, jMax = 1000)
matlines(known_tree, lwd = 2, lty = 2, col = "black")

