#Jana Fruth (2016)

support <- function(fun, d, xi = 1, h = 0.01, n = 5000, n.points = 50, q, q.arg,...)
{
  # distribution defaults
  if (is.null(q)) {
    q <- rep("qunif", d)
  }
  else if (length(q) == 1) {
    q <- rep(q, d)
  }
  if (is.null(q.arg)) {
    q.arg <- rep(list(list()), d)
  }
  else if (FALSE %in% sapply(q.arg, is.list)) {
    q.arg <- rep(list(q.arg), d)
  }
  
  # define support points
  trager <- do.call(q[xi], c(list(c(0,1)),q.arg[[xi]]))
  points <- seq(trager[1], trager[2],, n.points)
  
  main <- numeric(n.points)
  total <- numeric(n.points)
  
  for (p in 1:n.points){
    # mc-sample
    
    X01 <- matrix(runif(n * d), ncol = d)
    X <- matrix(, n, d)
    for (j in 1:d) X[, j] <- do.call(q[j], c(list(p = X01[,j]), q.arg[[j]]))
    
    Xplus <- X
    Xplus[,xi] <- points[p] + h
    Xminus <- X
    Xminus[,xi] <- points[p] - h
    
    # finite differences
    main[p] <- (mean((fun(Xplus) - fun(Xminus))/(2*h)))^2
    total[p] <- mean(((fun(Xplus) - fun(Xminus))/(2*h))^2)
    
  }
  col1 <- "lightskyblue1"
  colT <- "lightskyblue4"
  plot(points, main, type = "n", xlab="distribution support", ...)
  polygon(c(points[1],points,points[n.points]), c(0,total,0),col = colT)
  polygon(c(points[1],points,points[n.points]), c(0,main,0),col = col1)
  return(list(main=main, total=total))
}