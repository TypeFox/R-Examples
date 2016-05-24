##
## Unit testing of Geometric Program
test.GP <- function(){
    F0 <- matrix(c(3, -2, -1, 0, 1, -3), nrow = 3, ncol = 2, byrow = TRUE)
    g0 <- log(c(0.44, 10, 0.592))
    F1 <- matrix(c(-1, 3), nrow = 1, ncol = 2, byrow = TRUE)
    g1 <- log(8.62)
    ans <- gp(F0, g0, FList = list(F1), gList = list(g1))
    checkTrue(ans$status == "optimal")
    return()
}
