setClass(Class="test.problem",
         representation=representation(
           name="character",
           f="function",
           grad="function",
           n="integer",
           maxf="integer",
           objective="numeric",
           ntest="integer",
           lower="numeric",
           upper="numeric"))

setMethod(f="show",
          signature="test.problem",
          definition=function(object) {
            cat("Test problem: ",object@name,"\n",sep="")
            cat("Dimension: ",object@n,"\n",sep="")
            cat("Objective: ",object@objective,"\n",sep="")
            cat("Max function eval.: ",object@maxf,"\n",sep="")
            cat("Number of tests: ",object@ntest,"\n",sep="")
            invisible(object)
          })

test.problem <- function(name, n.test=100, dim, maxf, objective, lower, upper) {
  parabola <- function(x) {
    return(sum(x*x))
  }
  parabola.grad <- function(x) {
    return(2*x)
  }
  parabola.dim <- 30
  parabola.maxf <- 50000
  parabola.objective <- 1e-8
  parabola.lower <- rep(-100,parabola.dim)
  parabola.upper <- rep(100,parabola.dim)

  griewank <- function(x) {
    return(sum(x*x)/4000-prod(cos(x/sqrt(1:n)))+1)
  }
  griewank.grad <- function(x) {
    z <- cos(x/sqrt(1:n))
    dz <- -sin(x/sqrt(1:n))/sqrt(1:n)
    return(x/2000-dz*sapply(1:n,function(i) prod(z[-i])))
  }
  griewank.dim <- 2
  griewank.maxf <- 30000
  griewank.objective <- 1e-3
  griewank.lower <- rep(-600,griewank.dim)
  griewank.upper <- rep(600,griewank.dim)
  
  rosenbrock <- function(x) {
    t0 <- x+1
    t1 <- t0[2:n]-t0[1:(n-1)]*t0[1:(n-1)]
    return(1e2*sum(t1*t1) + sum(x[1:(n-1)]*x[1:(n-1)]))
  }
  rosenbrock.grad <- function(x) {
    t0 <- x+1
    t1 <- t0[2:n]-t0[1:(n-1)]*t0[1:(n-1)]
    return(c(-4e2*t1*t0[1:(n-1)]+2*t0[1:(n-1)],0)+c(0,200*t1))
  }
  rosenbrock.dim <- 20
  rosenbrock.maxf <- 100000
  rosenbrock.objective <- 1e-4
  rosenbrock.lower <- rep(-10,rosenbrock.dim)
  rosenbrock.upper <- rep(10,rosenbrock.dim)

  rastrigin <- function(x) {
    return(10*n+sum(x^2-10*cos(2*pi*x)))
  }
  rastrigin.grad <- function(x) {
    return(2*x+20*pi*sin(2*pi*x))
  }
  rastrigin.dim <- 2
  rastrigin.maxf <- 3000
  rastrigin.objective <- 0
  rastrigin.lower <- rep(-5.12,rastrigin.dim)
  rastrigin.upper <- rep(5.12,rastrigin.dim)
  
  ackley <- function(x) {
    return(-20*exp(-.2*sqrt(sum(x^2)/n))-exp(sum(cos(2*pi*x))/n)+20+exp(1))
  }
  ackley.grad <- function(x) {
    return((4*exp(-.2*sqrt(sum(x^2)/n))/sqrt(sum(x^2)*n))*x+(exp(sum(cos(2*pi*x))/n)*2*pi/n)*sin(2*pi*x))
  }
  ackley.dim <- 10
  ackley.maxf <- 5000
  ackley.objective <- 1e-4
  ackley.lower <- rep(-32,ackley.dim)
  ackley.upper <- rep(32,ackley.dim)

  name <- match.arg(name,c("parabola","griewank","rosenbrock","rastrigin",
                           "ackley"),FALSE)
  n <- ifelse(missing(dim),eval(as.name(paste(name,"dim",sep="."))),dim)
  objective <- ifelse(missing(objective),
                      eval(as.name(paste(name,"objective",sep="."))),
                      objective)
  e <- new.env()
  assign("n",n,envir=e)
  f <- eval(as.name(name))
  environment(f) <- e
  g <- eval(as.name(paste(name,"grad",sep=".")))
  if (is.null(g)) {
    g <- function(x,...) { require(numDeriv); return(grad(f,x,...)) }
  } else {
    environment(g) <- e
  }
  maxf <- ifelse(missing(maxf),eval(as.name(paste(name,"maxf",sep="."))),maxf)
  lower <- ifelse(missing(lower),
                  eval(as.name(paste(name,"lower",sep="."))),lower)
  upper <- ifelse(missing(upper),
                  eval(as.name(paste(name,"upper",sep="."))),upper)
  return(new(Class="test.problem",
             name=name,
             f=f,
             grad=g,
             n=as.integer(n),
             maxf=as.integer(maxf),
             objective=as.double(objective),
             ntest=as.integer(n.test),
             lower=as.double(lower),
             upper=as.double(upper)))
}

