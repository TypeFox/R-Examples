##' @name ZDT1
##' @aliases ZDT2
##' @aliases ZDT3
##' @aliases ZDT4
##' @aliases ZDT6
##' @aliases P1
##' @aliases P2
##' @aliases MOP2
##' @aliases MOP3
##' @aliases DTLZ1
##' @aliases DTLZ2
##' @aliases DTLZ3
##' @aliases DTLZ7
##' 
##' @description Multi-objective test functions.
##' @title Test functions of x
##'
##' @param x matrix specifying the location where the function is to be evaluated, one point per row,
##' @param nobj optional argument to select the number of objective for the DTLZ test functions.
##'
## ' @note \code{funs} is a generic name for the functions documented.
##' 
##' @details These functions are coming from different benchmarks:
##' the \code{ZDT} test problems from an article of E. Zitzler et al., \code{P1} from the thesis of J. Parr and \code{P2}
##' from an article of Poloni et al. . \code{MOP2} and \code{MOP3} are from Van Veldhuizen and \code{DTLZ} functions are from Deb et al. . \cr \cr
##' 
##' Domains (sometimes rescaled to \code{[0,1]}):
##' \itemize{
##' \item \code{ZDT1-6}: \code{[0,1]^d} 
##' \item \code{P1}, \code{P2}: \code{[0,1]^2} 
##' \item \code{MOP2}: \code{[0,1]^d}
##' \item \code{MOP3}: \code{[-3,3]}, tri-objective, 2 variables
##' \item \code{DTLZ1-3,7}: \code{[0,1]^d}, m-objective problems, with at least \code{d>m} variables. 
##' }

##' 
##' @return Matrix of values corresponding to the objective functions, the number of colums is the number of objectives.
##' 
##' @references
##' 
##' J. M. Parr (2012), \emph{Improvement Criteria for Constraint Handling and Multiobjective Optimization}, University of Southampton, PhD thesis. 
##' 
##' C. Poloni, A. Giurgevich, L. Onesti, V. Pediroda (2000), Hybridization of a multi-objective genetic algorithm, a neural network and a classical optimizer for a complex design problem in fluid dynamics, \emph{Computer Methods in Applied Mechanics and Engineering}, 186(2), 403-420.
##' 
##' E. Zitzler, K. Deb, and L. Thiele (2000), Comparison of multiobjective evolutionary
##' algorithms: Empirical results, \emph{Evol. Comput.}, 8(2), 173-195.
##' 
##' K. Deb, L. Thiele, M. Laumanns and E. Zitzler (2002), Scalable Test Problems for Evolutionary Multiobjective Optimization, 
##' \emph{IEEE Transactions on Evolutionary Computation}, 6(2), 182-197.
##' 
##' D. A. Van Veldhuizen, G. B. Lamont (1999), Multiobjective evolutionary algorithm test suites, \emph{In Proceedings of the 1999 ACM symposium on Applied computing}, 351-357.
##' 
##' 
##' 
##' @rdname TestFunctions
##' @export
##' @examples 
##' # ----------------------------------
##' # 2-objectives test problems
##' # ---------------------------------- 
##' 
##' plotParetoGrid("ZDT1", n.grid = 21)
##' 
##' plotParetoGrid("ZDT2", n.grid = 21)
##' 
##' plotParetoGrid("ZDT3", n.grid = 21)
##' 
##' plotParetoGrid("ZDT4", n.grid = 21)
##' 
##' plotParetoGrid("ZDT6", n.grid = 21)
##' 
##' plotParetoGrid("P1", n.grid = 21)
##' 
##' plotParetoGrid("P2", n.grid = 21)
##' 
##' plotParetoGrid("MOP2", n.grid = 21)
##' 
##' @rdname TestFunctions
##' @export
ZDT1 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    n <- ncol(x)
    g <- 1+rowSums(x[,2:n, drop = FALSE])*9/(n-1)
    return(cbind(x[,1],g*(1-sqrt(x[,1]/g))))
  }

##' @rdname TestFunctions
##' @export
ZDT2 <- 
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    n <- ncol(x)
    g <- 1+rowSums(x[,2:n, drop = FALSE])*9/(n-1)
    return(cbind(x[,1],g*(1-(x[,1]/g)^2)))
  }

##' @rdname TestFunctions
##' @export
ZDT3 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    n <- ncol(x)
    g <- 1+rowSums(x[,2:n, drop = FALSE])*9/(n-1)
 
    return(cbind(x[,1], g*(1 - sqrt(x[,1]/g) - x[,1]/g*sin(10*pi*x[,1]))))
  }

##' @rdname TestFunctions
##' @export
ZDT4 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    n <- ncol(x)
    
    g <- 1+10*(n-1) + rowSums((x[,2:n, drop = FALSE]*10-5)^2-10*cos(4*pi*(x[,2:n, drop = FALSE]*10-5)))
    return(cbind(x[,1], g*(1 - sqrt(x[,1]/g))))
  }

##' @rdname TestFunctions
##' @export
ZDT6 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    n <- ncol(x)
    f1 <- 1-exp(-4*x[,1])*(sin(6*pi*x[,1]))^6
    g <- 1+9*(1/(n-1)*rowSums(x[,2:n, drop = FALSE]))^(0.25)
    return(cbind(f1,g*(1-(f1/g)^2)))
  }

##' @rdname TestFunctions
##' @export
P1 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, nrow = 1) 
    }
    b1<-15*x[,1]-5
    b2<-15*x[,2]
    return(cbind((b2-5.1*(b1/(2*pi))^2+5/pi*b1-6)^2 +10*((1-1/(8*pi))*cos(b1)+1),
                 -sqrt((10.5-b1)*(b1+5.5)*(b2+0.5)) - 1/30*(b2 -5.1*(b1/(2*pi))^2-6)^2 - 1/3*((1-1/(8*pi))*cos(b1)+1)
                 ) 
    )
  }

##' @rdname TestFunctions
##' @export
P2 <-
  function(x){
    if(is.null(dim(x))){
      x <- matrix(x, 1) 
    }
    xmod <- x*2*pi - pi
    ap <- matrix(c(0.5,1.5,1,2),2)
    bp <- matrix(c(-2,-1,-1.5,-0.5),2)
    A1 <- ap[1,1]*sin(1)+bp[1,1]*cos(1)+ap[1,2]*sin(2)+bp[1,2]*cos(2)
    A2 <- ap[2,1]*sin(1)+bp[2,1]*cos(1)+ap[2,2]*sin(2)+bp[2,2]*cos(2)
    
    B1 <- ap[1,1]*sin(xmod[,1])+bp[1,1]*cos(xmod[,1])+ap[1,2]*sin(xmod[,2])+bp[1,2]*cos(xmod[,2])
    B2 <- ap[2,1]*sin(xmod[,1])+bp[2,1]*cos(xmod[,1])+ap[2,2]*sin(xmod[,2])+bp[2,2]*cos(xmod[,2])
    F1 <- 1+(A1-B1)^2+(A2-B2)^2
    F2 <- (xmod[,1]+3)^2+(xmod[,2]+1)^2
    return(cbind(-F1,-F2))
  }

## ' @rdname TestFunctions
## ' @export
GSP <- function(x, gamma=1)
{
  N <- nrow(x)
  n <- ncol(x)
  
  if (is.null(N))
  { N <- 1
    n <- length(x)
    x <- matrix(x, nrow=1,ncol=n)
  }
  
  obj <- matrix(rep(0, 2*N), nrow=N, ncol=2)
  alpha <- 1/(2*gamma)
  
  obj[,1] <- (1/(n^alpha)) * (rowSums(x^2)^alpha)
  obj[,2] <- (1/(n^alpha)) * (rowSums((1-x)^2)^alpha)
  
  return(obj)
}

##' @rdname TestFunctions
##' @export
MOP2 <- function(x)
{
  xmod <- x*4 - 2
  if (is.null(nrow(x)))
  { 
    n <- length(xmod)
    y1 <- 1 - exp(-sum((xmod - 1/sqrt(n))^2) )
    y2 <- 1 - exp(-sum((xmod + 1/sqrt(n))^2) )
    Y <- matrix(c(y1,y2),2,1)
  } else
  {
    n <- ncol(xmod)
    y1 <- 1 - exp(-rowSums((xmod - 1/sqrt(n))^2) )
    y2 <- 1 - exp(-rowSums((xmod + 1/sqrt(n))^2) )
    Y <- cbind(y1,y2)
  }
  
  return(Y)
}

##' @rdname TestFunctions
##' @export
MOP3 <- function(x){
  if(is.null(dim(x))){
    x <- matrix(x, 1) 
  }
  # One point per row
  f1 <- 0.5 * (x[,1]^2 + x[,2]^2) + sin(x[,1]^2 + x[,2]^2)
  f2 <- (3 * x[,1] - 2 * x[,2] + 4)^2/8 + (x[,1] - x[,2] + 1)^2/27  + 15
  f3 <- 1/(x[,1]^2 + x[,2]^2 + 1) - 1.1 * exp(-x[,1]^2 - x[,2]^2)
  f <- cbind(f1, f2, f3)
}

##' @rdname TestFunctions
##' @export 
DTLZ1 <- function(x, nobj = 3){
  
  if(is.null(dim(x))){
    x <- matrix(x, 1) 
  }
  n <- ncol(x)
  
  y <- matrix(x[,1:(nobj-1)], nrow(x))
  z <- matrix(x[,nobj:n], nrow(x))
  
  g <- 100 * (n - nobj + 1 + rowSums((z-0.5)^2 - cos(20 * pi *(z - 0.5))))

  tmp <- t(apply(y, 1, cumprod))
  tmp <- cbind(t(apply(tmp, 1, rev)) ,1)
  
  tmp2 <- cbind(1, t(apply(1-y, 1, rev)))
  
  f <- tmp * tmp2 * 0.5* (1 + g)
  return(f)
}

##' @rdname TestFunctions
##' @export 
DTLZ2 <- function(x, nobj = 3){
  if(is.null(dim(x))){
    x <- matrix(x, 1) 
  }
  n <- ncol(x)
  
  y <- matrix(x[,1:(nobj-1)], nrow(x))
  z <- matrix(x[,nobj:n], nrow(x))
  
  g <- rowSums((z-0.5)^2)
  
  #   tmp <- c(rev(cumprod(cos(y * pi/2))), 1)
  #   tmp2 <- c(1, rev(sin(y * pi/2)))
  tmp <- t(apply(cos(y * pi/2), 1, cumprod))
  tmp <- cbind(t(apply(tmp, 1, rev)), 1)
  
  tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
  
  f <- tmp * tmp2 * (1 + g)
  
}


##' @rdname TestFunctions
##' @export 
DTLZ3 <- function(x, nobj = 3){
  if(is.null(dim(x))){
    x <- matrix(x, 1) 
  }
  n <- ncol(x)
  
  y <- matrix(x[,1:(nobj-1)], nrow(x))
  z <- matrix(x[,nobj:n], nrow(x))
  
  g <- 100 * (n - nobj + 1 + rowSums((z-0.5)^2 - cos(20 * pi *(z - 0.5))))
  
  tmp <- t(apply(cos(y * pi/2), 1, cumprod))
  tmp <- cbind(t(apply(tmp, 1, rev)), 1)
  
  tmp2 <- cbind(1, t(apply(sin(y * pi/2), 1, rev)))
  
  f <- tmp * tmp2 * (1 + g)
                     
}

##' @rdname TestFunctions
##' @export 
DTLZ7 <- function(x, nobj = 3){
  if(is.null(dim(x))){
    x <- matrix(x, 1) 
  }
  n <- ncol(x)
  
  y <- matrix(x[,1:(nobj-1)], nrow(x))
  z <- matrix(x[,nobj:n], nrow(x))
  
  g <- 1 + 9 * rowSums(z/(1:(n - nobj + 1)))
  
  tmp <- cbind(y, 1)
  
  tmp2 <- cbind(matrix(1,nrow(x),nobj - 1), (1 + g) * (nobj -  rowSums(y/(1+g) * (1 + sin(3 * pi * y)))))

#   tmp2 <- c(rep(1,(nobj - 1)), (1 + g) * (nobj -  sum(y/(1+g) * (1 + sin(3 * pi * y)))))
  
  f <- tmp * tmp2
  
}
