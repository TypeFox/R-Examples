
#to meet the standard 'fn' argument and specific name arguments, 

#constraint function for MBBEFD(a,b)
constrmbbefd <- function(x, fix.arg, obs, ddistnam)
{
  res <- x[1]*(1-x[2]) #a*(1-b) >= 0
  names(res) <- NULL
  res
}
#domain : (a,b) in (-1, 0) x (1, +Inf)
constrmbbefd1 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1]+1, -x[1], x[2]-1, x[1]*(1-x[2])) #-1 < a < 0, b > 1, a*(1-b) >= 0
  names(res) <- NULL
  res
}
constrmbbefd1jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 4, 2)
  j[1,] <- c(1, 0)
  j[2,] <- c(0, -1)
  j[3,] <- c(1, 0)
  j[4,] <- c(-x[2], -x[1])
  dimnames(j) <- NULL
  j
}  

#domain : (a,b) in (0, +Inf) x (0, 1)
constrmbbefd2 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1], x[1], 1-x[2], x[1]*(1-x[2])) #0 < a , 0 < b < 1, a*(1-b) >= 0
  names(res) <- NULL
  res
}
constrmbbefd2jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 4, 2)
  j[1,] <- c(1, 0)
  j[2,] <- c(0, 1)
  j[3,] <- c(0, -1)
  j[4,] <- c(-x[2], -x[1])
  dimnames(j) <- NULL
  j
}  


#constraint function for MBBEFD(g,b)
constrMBBEFD <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1]-1, x[2]) #g >= 1, b > 0
  names(res) <- NULL
  res
}
#domain : (g,b) in (1, +Inf) x (1, +Inf) with gb > 1
constrMBBEFD1 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1]-1, x[2]-1, x[1]*x[2]-1) #g > 1, b > 1, gb > 1
  names(res) <- NULL
  res
}
constrMBBEFD1jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 3, 2)
  j[1,] <- c(1, 0)
  j[2,] <- c(0, 1)
  j[3,] <- c(x[2], x[1])
  dimnames(j) <- NULL
  j
}

#domain : (g,b) in (1, +Inf) x (0, 1) with gb < 1
constrMBBEFD2 <- function(x, fix.arg, obs, ddistnam)
{
  res <- c(x[1]-1, 1-x[2], x[2], 1-x[1]*x[2]) #g > 1, 1 > b > 0, gb < 1
  names(res) <- NULL
  res
}
constrMBBEFD2jac <- function(x, fix.arg, obs, ddistnam)
{
  j <- matrix(0, 4, 2)
  j[1,] <- c(1, 0)
  j[2,] <- c(0, -1)
  j[3,] <- c(0, 1)
  j[4,] <- c(-x[2], -x[1])
  dimnames(j) <- NULL
  j
}

