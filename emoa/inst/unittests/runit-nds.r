##
## runit-nds.r - Pareto dominance stuff
##

points <- matrix(c(1.0, 0.0, 0.0,
                   0.0, 1.0, 0.0,
                   0.0, 0.0, 1.0,
                   0.5, 0.5, 0.5,
                   0.5, 0.6, 0.6,
                   0.6, 0.5, 0.6,
                   0.6, 0.6, 0.5,
                   0.8, 0.8, 0.8),
                 ncol=8, byrow=FALSE)

nd <- c(TRUE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
no <- c(1, 1, 1, 1, 2, 2, 2, 3)

test.is_dominated <- function() {
  k <- nrow(points)
  n <- ncol(points)
  ## Check for different permutations of the rows and columns of
  ## points.
  for (i in 1:10) {
    o <- sample(1:n)
    p <- sample(1:k)
    m <- points[p,o]
    res <- is_dominated(m)
    checkEquals(res, !nd[o])
  }  
}

test.is_maximally_dominated <- function() {
  k <- nrow(points)
  n <- ncol(points)
  ## Check for different permutations of the rows and columns of
  ## points.
  for (i in 1:10) {
    o <- sample(1:n)
    p <- sample(1:k)
    m <- points[p,o]
    res <- is_maximally_dominated(m)
    checkEquals(res, max(no[o]) == no[o])
  }  
}

test.nds_rank <- function() {
  k <- nrow(points)
  n <- ncol(points)
  ## Check for different permutations of the rows and columns of
  ## points.
  for (i in 1:10) {
    o <- sample(1:n)
    p <- sample(1:k)
    checkEquals(nds_rank(points[p,o]), no[o])
  }  
}

test.nds_rank.args <- function() {
  checkException(nds_rank("a"))
  checkException(nds_rank(1))
  checkException(nds_rank(list(1, 2, 3)))
  checkException(nds_rank(data.frame(x=1:10)))
  checkException(nds_rank(points, partial="a"))
}

## Bug fixed i
test.single_nds <- function() {
  checkEquals(dim(nondominated_points(points[,-(1:3)])), c(3, 1))
}

##test_dominates_op <- function() {
##  x <- c(1, 2, 1)
##  y <- c(2, 1, 2)
##  z <- c(0, 1, 0)
##  checkEquals(x %dominates% y, FALSE)
##  checkEquals(y %dominates% x, FALSE)
##  checkEquals(z %dominates% x, FALSE)
##  checkEquals(z %dominates% y, FALSE)
##}
