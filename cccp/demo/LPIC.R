##
## Demo for solving a Linear Programs with (linear) inequality constraints
## (Example taken from cvxopt's userguide)
##
## First example
## Creating LP
q <- c(-4, -5)
G <- matrix(c(2, 1, -1, 0,
              1, 2, 0, -1),
            nrow = 4, ncol = 2)
h <- c(3, 3, 0, 0)
nno1 <- nnoc(G = G, h = h)
## Using main function of package
ans <- cccp(q = q, cList = list(nno1), optctrl = ctrl())
ans
getx(ans)
## Second example
## Creating LP
q <- c(2, 1)
## linear constraints
G <- matrix(c(-1, 1,
              -1, -1,
              0, -1,
              1, -2),
            nrow = 4, ncol = 2, byrow = TRUE)
h <- matrix(c(1, -2, 0, 4))
nno1 <- nnoc(G = G, h = h)
ans <- cccp(q = q, cList = list(nno1))
ans
getx(ans)
