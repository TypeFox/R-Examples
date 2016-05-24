binom.blaker.acc.single.p <- function(x,n,p,type="orig",acc.tol=1e-10,output="acc",maxiter=100) {
# Reduce calculations to p <= x/n and search, thus,
# below just for the LEFT interval limits.
  if (p > x/n) {
    x <- n - x
    p <- 1 - p
  }
# "Ordinary" acceptability at p.
  p1.p <- pbinom(x-1,n,p,lower.tail=FALSE)
  q1.p <- qbinom(p1.p,n,p)-1
  a1.p <- min(p1.p + pbinom(q1.p,n,p),1)
# "Unimodalization"
  if (type == "unimod" && q1.p >= 0) {
    upper <- p
    a1.upp <- a1.p
    lower <- 0
    a1.low <- 1
#
#   Search for the first discontinuity point p* of the
#   acceptability function to the left of p.
#   Continue as long as there is a chance to find there
#   a higher acceptability value than at p.
#
#   More detailed explanation:
#
#   Let f_{ij}(p) = P_p(X >= i) + P_p(X <= j),
#   a strictly quasiconvex function on [0, 1] for i > j+1.
#   Acceptability equals f_{ij}(p) for some i > j+1
#   on intervals of continuity except those where equal to 1.
#   As soon as a p' with f_{ij}(p') <= f_{ij}(p) is found
#   to the left of p*, it becomes clear, due to quasiconvexity,
#   that the acceptability in p* is not bigger than in p,
#   and its maximum on [0, p] lies in p.
#
    iter <- 0
#   In 1.0-4, ... >= acc.tol replaced with ... > acc.tol
    while (a1.low > a1.p && a1.low - a1.upp > acc.tol) {
      iter <- iter + 1
      if (iter > maxiter) {
        warning("Convergence not attained after ",maxiter, 
                             " iterations for n = ",n,", x = ",x,", p = ",p,
                             ", and acceptability tolerance limit of ",acc.tol)
        break
      }
      mid <- (lower+upper)/2
      p1.mid <- pbinom(x-1,n,mid,lower.tail=FALSE)
      p2.mid <- pbinom(q1.p,n,mid)
      a1.mid <- p1.mid + p2.mid
      if (p1.mid < p2.mid) {
        lower <- mid
        a1.low <- a1.mid
      }
      else {
        upper <- mid
        a1.upp <- a1.mid
      }
    }
    a1.p <- max(a1.p,a1.low)
  }
  if (output == "acc") {
    return(a1.p)
  }
  else if (output == "both") {
    return(c(a1.p,q1.p))
  }
  else if (output == "q1") {
    return(q1.p)
  }
}
