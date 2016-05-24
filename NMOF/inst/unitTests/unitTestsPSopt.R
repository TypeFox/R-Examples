## PSopt
test.PSopt <- function() {
    ## test function: PS should find minimum
    OF <- tfRosenbrock
    size <- 3L ## define dimension
    algo <- list(printBar = FALSE,
                 printDetail = FALSE,
                 nP = 50L, nG = 1000L,
                 c1 = 0.0, c2 = 1.5,
                 iner = 0.8, initV = 0.50, maxV = 50,
                 min = rep(-50, size),
                 max = rep( 50, size))

    sol <- PSopt(OF = OF, algo = algo)
    checkEquals(sol$OFvalue, 0)

    ## exception: wrong size of initP
    algo$initP <- array(0, dim = c(20L,20L))
    checkException(res <- PSopt(OF = OF, algo), silent = TRUE)
    algo$initP <- function()
        array(0, dim = c(5,20))
    checkException(res <- PSopt(OF = OF, algo), silent = TRUE)
    algo$initP <- NULL


    ## check if Fmat/xlist are returned
    ## ...if FALSE
    algo <- list(printBar = FALSE,
                 printDetail = FALSE,
                 nP = 50L, nG = 100L,
                 c1 = 0.0, c2 = 1.5,
                 iner = 0.8, initV = 0.50, maxV = 50,
                 min = rep(-50, size),
                 max = rep( 50, size),
                 storeF = FALSE,
                 storeSolutions = FALSE)
    sol <- PSopt(OF = OF, algo = algo)
    checkTrue(is.na(sol$Fmat))
    checkTrue(is.na(sol$xlist))
    checkEquals(length(sol$Fmat), 1L)
    checkEquals(length(sol$xlist), 1L)

    ## ...if TRUE
    algo$storeF <- TRUE
    algo$storeSolutions <- TRUE
    sol <- PSopt(OF = OF, algo = algo)
    checkEquals(names(sol$xlist), c("P","Pbest", "initP"))
    checkEquals(dim(sol$xlist[[c(1L, algo$nG)]]),
                c(length(algo$min),algo$nP))
    checkEquals(dim(sol$xlist[[c(2L, algo$nG)]]),
                c(length(algo$min),algo$nP))
    ## xlist has only three elements
    checkException(sol$xlist[[c(4L, algo$nG)]], silent = TRUE)
    ## xlist[[i]] stores only algo$nG elements
    checkException(sol$xlist[[c(2L, algo$nG + 1L)]], silent = TRUE)

    ## a trivial problem: 'recover' a vector X
    X <- 1:10 - 5
    OF <- function(x, X)
        sum(abs(x - X))
    algo <- list(nP = 200, nG = 1000, c2=1, c1=1,
                 iner = 0.9,
                 min = rep(-4, length(X)),
                 max = rep( 4,  length(X)),
                 minmaxConstr = FALSE,
                 printBar = FALSE, printDetail = FALSE,
                 storeF = FALSE, storeSolutions = FALSE)
    sol <- PSopt(OF, algo, X)
    checkEquals(round(sol$xbest,3), X)
    algo$minmaxConstr <- TRUE
    sol <- PSopt(OF, algo, X)
    checkEquals(round(sol$xbest,3), c(-4, -3, -2, -1, 0, 1, 2, 3, 4, 4))

    ## ERROR: vectorised comp and wrong dim
    X <- 1:10 - 5
    OF <- function(x, X)  ## correct
        colSums(abs(x - X))
    OF <- function(x, X) {  ## WRONG
        res <- colSums(abs(x - X))
        res <- res[-1]
    }
    algo <- list(nP = 20, nG = 2,
                 min = rep(-3, length(X)),
                 max = rep( 3,  length(X)),
                 loopOF = FALSE, loopRepair = FALSE, loopPen = FALSE,
                 printBar = FALSE, printDetail = FALSE,
                 storeF = FALSE, storeSolutions = FALSE)
    checkException(sol <- PSopt(OF, algo, X), silent = TRUE)

    OF <- function(x, X)  ## correct
        colSums(abs(x - X))
    algo$repair <- function(x, X) {
        x[,-1]
    }
    checkException(sol <- PSopt(OF, algo, X), silent = TRUE)

}
