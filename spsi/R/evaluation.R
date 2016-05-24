# This function as an input takes univariate version of object 'spline' and returns executable function.
# der.x = vector or numeric specifying which derivatives are required. At the moment only original function
# (der.x = 0) or first derivative (der.x = 1) are supported.

sps_eval_u <- function (x, spline, der.x=0)
{

  der.x <- ifelse(is.null(der.x),0,der.x)
  N <- length(der.x)
  names <- c('f','fx')

  values <- matrix(ncol=N, nrow = length(x))

  for(l in 1:length(x)) {

    i <- SEARCH(spline$x,x[l])
    Bx <- Base(x[l],spline$x[i],
               spline$x[i+1],
               spline$deg[i],
               spline$k,der.x)
    for(k in 1:N)
      values[l,k] <- patch_U(spline$y,spline$fx,i,Bx[[der.x[k]+1]])
  }

  if(N == 1) as.vector(values)
  else{
    result <- lapply(seq_len(ncol(values)), function(i) values[,i])
    names(result) <- names[1 + der.x]
    result
  }
}

# This function as an input takes univariate version of object 'spline' and returns executable function.
# der.x, der.y = vectors or numerics specifying which derivatives are required. At the moment only original function
# (der.x = 0, der.y=0), first partial derivative (der.x = 1. der.y=0 or der.x=0, der.y=1) and mixed partial derivative
# (der.x = 1, der.y = 1) are supported.

sps_eval_bi = function(x, y, spline, der.x=0, der.y=0, grid)
{

  stopifnot(class(x) == class(y))

  if(grid && is.vector(x)){
    mesh <- mesh(x,y)
    x <- mesh$x
    y <- mesh$y
  }
  else
    stopifnot(all(dim(as.matrix(y)) == dim(as.matrix(x))))

  if(is.null(der.x)) der.x <- 0
  if(is.null(der.y)) der.y <- 0

  stopifnot(length(der.x) == length(der.y))

  names <- c('f','fx','fy','fxy')

  N <- length(der.x)

  to.matrix <- ifelse(is.matrix(x), T, F)

  if(is.matrix(x)){
    numcol <- ncol(x)
    x <- as.vector(x)
    y <- as.vector(y)
  }

  values <- matrix(ncol=N, nrow = length(x))

  for(l in 1:length(x)){
    i <- SEARCH(spline$x, x[l])
    j <- SEARCH(spline$y, y[l])

    Bx <- Base(x[l], spline$x[i],
               spline$x[i+1], spline$deg.x[i],
               spline$k, der.x)
    By <- Base(y[l], spline$y[j],
               spline$y[j+1], spline$deg.y[j],
               spline$k, der.y)

    for (k in 1:N)
      values[l,k] <- patch(spline$z, spline$fx,
                           spline$fy, spline$fxy,
                           i, j, Bx[[der.x[k]+1]],
                           By[[der.y[k]+1]])
  }
  if(to.matrix){
    if(N == 1)
      matrix(values, ncol = numcol)
    else{
      result <- list()
      for(k in 1:N)
        result[[k]] <-  matrix(values[,k], ncol = numcol)
      names(result) <- names[1 + der.x + 2 * der.y]
      result
    }
  }
  else {
    if(N == 1)
      as.vector(values)
    else {
      result <- lapply(seq_len(ncol(values)), function(i) values[,i])
      names(result) <- names[1 + der.x + 2 * der.y]
      result
    }
  }
}


# This function evaluates one-dimensional basis functions. First 4 piecewise-linear
# functions li, i=1,2,3,4, are evaluated. These functions have knots at x0,t1,t2,x1
# and are uniquely defined by the following properties:
# l1(x0) = 1; l1(x1) = Dl1(x0) = Dl1(x1) = 0
# l2(x1) = 1; l2(x0) = Dl2(x0) = Dl2(x1) = 0
# Dl3(x0) = 1, l2(x0) = l2(x1) = Dl3(x1) = 0
# Dl4(x1) = 1; l4(x0) = l4(x1) = Dl4(x0) = 0
#
# Function bi is a Berstain polynomials of degree n of function li, i=1,2,3,4;
# Berstain polynomial preserves unique property of piecewise linear functions
#
#           n
# bi(x) = (SUM n C nu * li(x(nu)) * (x-x0)^nu * (x1-x)^(n-nu))/(x1-x0)^n
#          nu=0

# Function uses the fact that b1(x) = 1-b2(x); b3(x) = (x-x0) - b4(x) - (x1-x0)*b2(x)
#
# b1_prime, i = 1,2,3,4 are first derivatives of corresponding functions b and are
# computed if user requires derivative of spline.
#
# Function is written in terms of x but can be equally used in y-dimension.
Base <- function(x, x0, x1, n, k, der.x)
{
  nu <- 1:n
  t1 <- x0 + k * (x1 - x0) / n
  t2 <- x1 - k * (x1 - x0) / n
  x_nu <- x0 + nu * (x1 - x0) / n

  l2 <- ifelse(x_nu <= t1, 0,
               ifelse(x_nu < t2, (x_nu-t1) / (t2-t1),
                      1))
  l4 <- ifelse(x_nu <= t1, 0,
               ifelse(x_nu < t2, (t2-x1)*(x_nu-t1) / (t2-t1),
                      x_nu - x1))

  BT <- choose(n,nu)
  result <- list()
  if(any(der.x == 0)) {
    b2 <- sum(BT * l2 * (x - x0)^nu * (x1 - x)^(n - nu))/(x1 - x0)^n
    b4 <- sum(BT * l4 * (x - x0)^nu * (x1 - x)^(n - nu))/(x1 - x0)^n
    b1 <- 1 - b2
    b3 <- (x - x0) - b4 - (x1 - x0) * b2
    result[[1]] <- c(b1,b2,b3,b4)
  }

  if(any(der.x == 1)) {
    b2_prime <- sum(BT * l2 * (x - x0)^(nu - 1) * (x1 - x)^(n - nu - 1) *
                      (nu * x1 + (n - nu) * x0 - n * x))/(x1 - x0)^n
    b4_prime <- sum(BT * l4 * (x - x0)^(nu - 1) * (x1 - x)^(n - nu - 1) *
                      (nu * x1 + (n - nu) * x0 - n*x))/(x1 - x0)^n
    b1_prime <- - b2_prime
    b3_prime <- 1 - b4_prime - (x1 - x0) * b2_prime
    result[[2]] <- c(b1_prime, b2_prime, b3_prime, b4_prime)

  }
  result
}



# This function returns index of grid point x_grid s.t x_grid <= x_tab unless
# x_tab falls outside grid, then function returns index of last point
SEARCH <- function(X, x)
{

  ind <- ifelse(x < min(X), 1,
                ifelse(x >= max(X), length(X) - 1,
                       max(which(X <= x))))
}


# This function evaluates the spline at given point. It makes use of the unique
# properties of basis functions b so that values of function and its derivatives,
# whether provided by user or computed in program, is preserved on the corners of
# rectangle.
patch <-  function(Z, Fx, Fy, Fxy, i, j, bx, by)
{
  a <- 0:1
  sum(Z[i+a,j] * bx[1+a] * by[1],
      Fx[i+a,j] * bx[3+a] * by[1],
      Fy[j,i+a] * bx[1+a] * by[3],
      Fxy[i+a,j] * bx[3+a] * by[3],
      Z[i+a,j+1] * bx[1+a] * by[2],
      Fx[i+a,j+1] * bx[3+a] * by[2],
      Fy[j+1,i+a] * bx[1+a] * by[4],
      Fxy[i+a,j+1] * bx[3+a] * by[4]
  )
}

patch_U <- function(Z,fx,i,bx)
{

  Z[i] * bx[1] +
  Z[i+1] * bx[2] +
  fx[i] * bx[3] +
  fx[i+1] * bx[4]
}
