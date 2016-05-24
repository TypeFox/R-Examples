##
##  n e w t o n _ c o t e s . R  Newton-Cotes Formulas
##


cotes <- function(f, a, b, n = 20, nodes = 5, ...) {
    if (nodes < 2 || nodes > 8)
        stop("Number of nodes, 'nodes', must be between 2 and 8.")
    if (n < nodes)
        stop("Argument 'n' must be greater or equal to the number of nodes.")

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    N <- (nodes-1) * ceiling(n/(nodes-1))
    N1 <- N+1
    x<- linspace(a,b,N1)
    h <- x[2] - x[1]
    g <- f(x)
    if (length(g) !=  N1) {
        g <- numeric(N1)
        for (i in 1:N1)  g[i] <- f(x[i])
    }

    endpts <- g[1] + g[N1]

    switch(nodes - 1,
        # (2) Trapezoidal rule
        Q <- (h/2) * (endpts + 2*sum(g[2:N])),

        # (3) Simpson-Kepler rule
        Q <- (h/3) * (endpts + 4*sum(g[seq(2,N,2)]) + 2*sum(g[seq(3,N,2)])),

        # (4) Simpson's 3/8 rule
        Q <- (3*h/8) * (endpts + 3*sum(g[seq(2,N,3)] + g[seq(3,N,3)]) +
              2*sum(g[seq(4,N,3)])),

        # (5) Boole's 4/90 rule
        Q <- (4*h/90) * (7*endpts + 32*sum(g[seq(2,N,4)]) +
              12*sum(g[seq(3,N,4)]) + 32*sum(g[seq(4,N,4)]) +
              14*sum(g[seq(5,N,4)])),

        # (6) Five-point rule  
        Q <- (5*h/288) * (19*endpts + 75*sum(g[seq(2,N,5)] + g[seq(5,N,5)]) +
              50*sum(g[seq(3,N,5)] + g[seq(4,N,5)]) + 38*sum(g[seq(6,N,5)])),
        
        # (7) Weddle rule
        Q <- (6*h/840) * (41*endpts + 216*sum(g[seq(2,N,6)] + g[seq(6,N,6)]) +
              27*sum(g[seq(3,N,6)] + g[seq(5,N,6)]) + 272*sum(g[seq(4,N,6)]) +
              82*sum(g[seq(7,N,6)])),
        
        # (8) (no name)
        Q <- (7*h/17280) * (751*endpts + 3577*sum(g[seq(2,N,7)] + g[seq(7,N,7)]) +
              1323*sum(g[seq(3,N,7)] + g[seq(6,N,7)]) + 2989*sum(g[seq(4,N,7)] +
              g[seq(5,N,7)]) + 1502*sum(g[seq(8,N,7)]))
        )

    return(Q)
}
