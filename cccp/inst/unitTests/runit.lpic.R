##
## Unit testing of Linear Program with inequality constraints
test.LPIC <- function(){
    ## First example
    q <- c(-4, -5)
    G <- matrix(c(2, 1, -1, 0,
                  1, 2, 0, -1),
                nrow = 4, ncol = 2)
    h <- c(3, 3, 0, 0)
    nno1 <- nnoc(G = G, h = h)
    ans <- cccp(q = q, cList = list(nno1), optctrl = ctrl())
    checkTrue(ans$status == "optimal")
    ## Second example
    q <- c(2, 1)
    G <- matrix(c(-1, 1,
                  -1, -1,
                  0, -1,
                  1, -2),
                nrow = 4, ncol = 2, byrow = TRUE)
    h <- matrix(c(1, -2, 0, 4))
    nno1 <- nnoc(G = G, h = h)
    ans <- cccp(q = q, cList = list(nno1))
    checkTrue(ans$status == "optimal")
    return()
}
