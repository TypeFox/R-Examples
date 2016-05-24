######################### Code Source ##########################
#
# All Code and Comments Below are
# Copyright Marc Schwartz
# e-mail: marc_schwartz@me.com
# This code is made available under the GNU Public License V2.0
# This is free software and comes with ABSOLUTELY NO WARRANTY.

# R Code to calculate various Measures of
# Association for m x n tables
# References include:
# 1. Agresti A. (2002): Categorical Data Analysis, Second Edition, J. Wiley and Sons
# 2. Stokes M., Davis C. & Koch G. (1997): Categorical Data Analysis Using the SAS System, SAS Institute
# 3. Liebetrau A.M. (1983): Measures of Association (Sage University Papers Series on Quantitative Applications
#    in the Social Sciences, Series no. 07-032), Sage Publications
# 4. SAS Institute (1999): SAS/STAT User's Guide V8, SAS Institute
# 5. SPSS, Inc. (2003): SPSS 11.5 Statistical Algorithms
#    (http://www.spss.com/tech/stat/Algorithms/11.5/crosstabs.pdf)
# 6. Sheskin DJ. (2004): Handbook of Parametric and Nonparametric Statistical Procedures, Chapman & Hall/CRC


# MOST MEASURES TAKE A 2 DIMENSIONAL TABLE/MATRIX "x" AS AN ARGUMENT

# See the 'vcd' CRAN package for some examples and code
# on calculations and p values


#### Work horses ####

# Calculate CONcordant Pairs in a table
# cycle through x[r, c] and multiply by
# sum(x elements below and to the right of x[r, c])
# x = table
concordant <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  # get sum(matrix values > r AND > c)
  # for each matrix[r, c]
  mat.lr <- function(r, c)
  { 
    lr <- x[(r.x > r) & (c.x > c)]
    sum(lr)
  }

  # get row and column index for each
  # matrix element
  r.x <- row(x)
  c.x <- col(x)

  # return the sum of each matrix[r, c] * sums
  # using mapply to sequence thru each matrix[r, c]
  sum(x * mapply(mat.lr, r = r.x, c = c.x))
}


# Calculate DIScordant Pairs in a table
# cycle through x[r, c] and multiply by
# sum(x elements below and to the left of x[r, c])
# x = table
discordant <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  # get sum(matrix values > r AND < c)
  # for each matrix[r, c]
  mat.ll <- function(r, c)
  { 
    ll <- x[(r.x > r) & (c.x < c)]
    sum(ll)
  }

  # get row and column index for each
  # matrix element
  r.x <- row(x)
  c.x <- col(x)

  # return the sum of each matrix[r, c] * sums
  # using mapply to sequence thru each matrix[r, c]
  sum(x * mapply(mat.ll, r = r.x, c = c.x))
}


# Calculate Pairs tied on Rows
# cycle through each row of x and multiply by
# sum(x elements to the right of x[r, c])
# x = table
ties.row <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  total.pairs <- 0

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  for (r in 1:rows)
  {
    for (c in 1:(cols - 1))
    {
      total.pairs <- total.pairs + (x[r, c] * sum(x[r, (c + 1):cols]))
    }
  }

  total.pairs
}

# Calculate Pairs tied on Columns
# cycle through each col of x and multiply by
# sum(x elements below x[r, c])
# x = table
ties.col <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  total.pairs <- 0

  rows <- dim(x)[1]
  cols <- dim(x)[2]

  for (c in 1:cols)
  {
    for (r in 1:(rows - 1))
    {
      total.pairs <- total.pairs + (x[r, c] * sum(x[(r + 1):rows, c]))
    }
  }

  total.pairs
}

#### Phi ####

# Calculate Phi Coefficient
# x = table
calc.phi <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  phi <- sqrt(chisq.test(x, correct = FALSE)$statistic / sum(x))
  as.numeric(phi)
}

####  Contingency Coefficient (Pearson's C) ####
# Calculate Contingency Coefficient (Pearson's C)
# and Sakoda's Adjusted Pearson's C
# x = table
calc.cc <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  # CC - Pearson's C
  chisq <- chisq.test(x, correct = FALSE)$statistic
  C <- sqrt(chisq / (chisq + sum(x)))

  # Sakoda's adjusted Pearson's C
  k <- min(dim(x))
  SC <- C / sqrt((k - 1) / k)

  CClist <- list(as.numeric(C), as.numeric(SC))
  names(CClist) <- c("Pearson.C", "Sakoda.C")

  CClist
}


#### Tshuprow's T ####

# Calculate Tshuprow's T
# Not meaningful for non-square tables
# For 2 x 2 tables T = Phi
# x = table
calc.TT <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  TT <- sqrt(chisq.test(x, correct = FALSE)$statistic /
       (sum(x) * sqrt((dim(x)[1] - 1) * (dim(x)[2] - 1))))

  as.numeric(TT)
}


#### Cramer's V ####

# For 2 x 2 tables V = Phi
# x = table
calc.CV <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  CV <- sqrt(chisq.test(x, correct = FALSE)$statistic /
       (sum(x) * min(dim(x) - 1)))

  as.numeric(CV)
}


#### Lambda ####

#' Calculate Lambda for nominal data tables.
#' @param x     A table object.
#' @return A named list with the three values:
#'  \item{lambda.cr}{The row variable is used as independent, the column variable as dependent variable.}
#'  \item{lambda.rc}{The column variable is used as independent, the row variable as dependent variable.}
#'  \item{lambda.symmetric}{Symmetric Lambda (the mean of both above).}
#' @export
#' @author Marc Schwartz, Mark Heckmann
#' @note The code for the calculation was supplied by Marc Schwartz (under GPL 2). 
#' Checked against SPSS results.
#' @examples {
#'   
#' }
#' 
nom.lambda <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  SumRmax <- sum(apply(x, 1, max))
  SumCmax <- sum(apply(x, 2, max))
  MaxCSum <- max(colSums(x))
  MaxRSum <- max(rowSums(x))
  n <- sum(x)

  L.CR <- (SumRmax - MaxCSum) / (n - MaxCSum)
  L.RC <- (SumCmax - max(rowSums(x))) / (n - MaxRSum)
  L.S <- (SumRmax + SumCmax - MaxCSum - MaxRSum) /
        ((2 * n) - MaxCSum - MaxRSum)

  res <- list(L.CR, L.RC, L.S)
  names(res) <- c("lambda.cr", "lambda.rc", "lambda.symmetric")
  class(res) <- "nom.lambda"
  res
}


#' Print method for class nom.lambda
#' 
#' @param x Object to be printed.
#' @param digits Number of decimal places to round to.
#' @param ...  Not evaluated.
#' @export
#' @keywords internal
#' @method print nom.lambda
#' 
print.nom.lambda <- function(x, digits=3, ...)
{ 
  cat("Lambda:\n")
  cat("\tColumns dependent:", formatC(x$lambda.cr, digits, format="f"), "\n")
  cat("\tRows dependent:", formatC(x$lambda.rc, digits, format="f"), "\n")  
  cat("\tSymmetric:", formatC(x$lambda.symmetric, digits, format="f"), "\n")  
}


#### Uncertainty Coefficient ####

#' Calculate the Uncertainty Coefficient (Theil's U)
#' 
#' @param x     A table object.
#' @return A named list with the three values:
#'  \item{ucc.cr}{The row variable is used as independent, the column variable as dependent variable.}
#'  \item{uc.rc}{The column variable is used as independent, the row variable as dependent variable.}
#'  \item{uc.symmetric}{Symmetric uncertainty coefficient.}
#' @export
#' @author Marc Schwartz, Mark Heckmann
#' @note The code for the calculation was supplied by Marc Schwartz (under GPL 2). 
#' Note: Asymmetric formulae denomiators corrected on May 4, 2007
#' thanks to Antti Arppe. Checked against SPSS results.
#' @examples {
#'   
#' }
#' 
nom.uncertainty <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  SumR <- rowSums(x)
  SumC <- colSums(x)
  n <- sum(x)

  HY <- -sum((SumC / n) * log(SumC / n))
  HX <- -sum((SumR / n) * log(SumR / n))
  HXY <- -sum((x / n) * log(x / n))

  UC.RC <- (HX + HY - HXY) / HX
  UC.CR <- (HY + HX - HXY) / HY
  UC.S <- 2 * (HX + HY - HXY) / (HX + HY)

  res <- list(UC.RC, UC.CR, UC.S)
  names(res) <- c("uc.rc", "uc.cr", "uc.symmetric")
  class(res) <- "nom.uncertainty"
  res
}


#' Print method for class nom.uncertainty
#' @param x Object to be printed.
#' @param digits Number of decimal places to round to.
#' @param ...  Not evaluated.
#' @export
#' @keywords internal
#' @method print nom.uncertainty
#' 
print.nom.uncertainty <- function(x, digits=3, ...)
{ 
  cat("Uncertainty coefficient:\n")
  cat("\tColumns dependent:", formatC(x$uc.cr, digits, format="f"), "\n")
  cat("\tRows dependent:", formatC(x$uc.rc, digits, format="f"), "\n")  
  cat("\tSymmetric:", formatC(x$uc.symmetric, digits, format="f"), "\n")  
}



#### Gamma (Goodman and Kruskal) ####

#' Calculate Goodman-Kruskal gamma for ordinal data tables.
#' @param x     A table object.
#' @return The gamma value.
#' @export
#' @author Marc Schwartz, Mark Heckmann
#' @note The code for the calculation was supplied by Marc Schwartz (under GPL 2).
#' Checked against SPSS results.
#' @examples {
#'  # TODO
#' }
#'
ord.gamma <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)

  gamma <- (c - d) / (c + d)

  class(gamma) <- "ord.gamma"
  gamma
}


#' Print method for class ord.gamma
#' @param x Object to be printed.
#' @param digits Number of decimal places to round to.
#' @param ...  Not evaluated.
#' @export
#' @keywords internal
#' @method print ord.gamma
#' 
print.ord.gamma <- function(x, digits=3, ...)
{
  cat("Goodman-Kruskal Gamma:", formatC(x, digits, format="f"), "\n")
}


#### Tau ####

# Calculate Kendall-Stuart Tau-c
# x = table
tau.a <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)
  n <- sum(x)
  # check if table is quadratic otherwise issue warning ?
  
  tau.a <- (c + d) / (n*(n-1) / 2)
  tau.a
}


# Calculate Kendall's Tau-b for ordinal data tables.
tau.b <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)

  # An alternative computation is:
  #Tr <- ties.row(x)
  #Tc <- ties.col(x)
  #KTb <- (c - d) / sqrt((c + d + Tc) * (c + d + Tr))

  # The "preferred" computation is:
  n <- sum(x)
  SumR <- rowSums(x)
  SumC <- colSums(x)

  tau.b <- (2 * (c - d)) / sqrt(((n ^ 2) - (sum(SumR ^ 2))) * ((n ^ 2) - (sum(SumC ^ 2))))
  tau.b
}


# Calculate Kendall-Stuart Tau-c
# x = table
tau.c <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)
  m <- min(dim(x))
  n <- sum(x)

  tau.c <- (m * 2 * (c - d)) / ((n ^ 2) * (m - 1))

  tau.c
}


#' Calculate Kendall's Tau statistics for ordinal data tables 
#' (Tau-b and Tau-c).
#' @param x     A table object.
#' @return A named list with the three values:
#'  \item{tau.a}{Tau-a satistic (for quadratic tables only)}
#'  \item{tau.b}{Tau-b satistic}
#'  \item{tau.c}{Kendall-Stuart Tau-c satistic}
#' @export
#' @author Marc Schwartz, Mark Heckmann
#' @note The code for the calculation was supplied by Marc Schwartz (under GPL 2)
#' @examples {
#'  # TODO
#' }
ord.tau <- function(x){
  l <- list(#tau.a=tau.a(x),  
            tau.b=tau.b(x),  
            tau.c=tau.c(x))  
  class(l) <- "ord.tau"
  l
}


#' Print method for class ord.tau
#' @param x Object to be printed.
#' @param digits Number of decimal places to round to.
#' @param ...  Not evaluated.
#' @export
#' @keywords internal
#' @method print ord.tau
#' 
print.ord.tau <- function(x, digits=3, ...){
  cat("Kendall's (and Stuart's) Tau statistics")
  #cat("\n\tTau-a:", formatC(x$tau.a, digits, format="f"))
  cat("\n\tTau-b:", formatC(x$tau.b, digits, format="f"))
  cat("\n\tTau-c:", formatC(x$tau.c, digits, format="f"))  
}


#### Somer's d ####

#' Calculate Somers' d for ordinal data tables.
#' @param x     A table object.
#' @return Kendall's Tau-b value.
#' @return A named list with the three values:
#'  \item{sd.cr}{The row variable is used as independent, the column variable as dependent variable.}
#'  \item{sd.rc}{The column variable is used as independent, the row variable as dependent variable.}
#'  \item{sd.symmetric}{Symmetric Somers' d.}
#' @export
#' @author Marc Schwartz, Mark Heckmann
#' @note The code for the calculation was supplied by Marc Schwartz (under GPL 2)
#' @examples {
#'  # TODO
#' }
#'
#'
ord.somers.d <- function(x)
{
  x <- matrix(as.numeric(x), dim(x))
  
  c <- concordant(x)
  d <- discordant(x)
  n <- sum(x)
  SumR <- rowSums(x)
  SumC <- colSums(x)

  sd.cr <- (2 * (c - d)) / ((n ^ 2) - (sum(SumR ^ 2)))
  sd.rc <- (2 * (c - d)) / ((n ^ 2) - (sum(SumC ^ 2)))
  sd.s <- (2 * (c - d)) / ((n ^ 2) - (((sum(SumR ^ 2)) + (sum(SumC ^ 2))) / 2))

  res <- list(sd.cr, sd.rc, sd.s)
  names(res) <- c("sd.cr", "sd.rc", "sd.symmetric")
  class(res) <- "ord.somersd"
  res
}


#' Print method for class somersd
#' @param x Object to be printed.
#' @param digits Number of decimal places to round to.
#' @param ...  Not evaluated.
#' @export
#' @keywords internal
#' @method print ord.somersd
#' 
print.ord.somersd <- function(x, digits=3, ...)
{
  cat("Somers' d:\n")
  cat("\tColumns dependent:", formatC(x$sd.cr, digits, format="f"), "\n")
  cat("\tRows dependent:", formatC(x$sd.rc, digits, format="f"), "\n")  
  cat("\tSymmetric:", formatC(x$sd.symmetric, digits, format="f"), "\n")  
}




#### Cochran's Q ####

# Test for proportions in dependent samples
# a k > 2 generalization of the mcnemar test
# 'mat' is a matrix, where:
# each row is a subject
# each column is the 0/1 result of a test condition
cochranq.test <- function(mat)
{
  k <- ncol(mat)

  C <- sum(colSums(mat) ^ 2)
  R <- sum(rowSums(mat) ^ 2)
  T <- sum(rowSums(mat))

  num <- (k - 1) * ((k * C) - (T ^ 2))
  den <- (k * T) - R

  Q <- num / den
  
  df <- k - 1
  names(df) <- "df"
  names(Q) <- "Cochran's Q"
  
  p.val <- pchisq(Q, df, lower.tail = FALSE)

  QVAL <- list(statistic = Q, parameter = df, p.value = p.val,
               method = "Cochran's Q Test for Dependent Samples",
               data.name = deparse(substitute(mat)))
  class(QVAL) <- "htest"
  return(QVAL)
}



#### Eta ####

#' Eta coefficient for nominal/interval data.
#' @param x Independent nominal variable (factor or numeric).
#' @param y Dependent interval variable (numeric).
#' @param breaks If \code{x} is interval data the \code{breaks} argument can be
#'   specified to classify the data. \code{breaks} is passed on to the function 
#'   \code{\link{cut}}.
#' @param na.rm Logical. Indicating if \code{NA} values are removed.
#' @return Eta coefficient
#' @export
#' @author Mark Heckmann
#' @examples 
#' attach(d.eta)     # using d.eta dataset
#' eta(x1, y)
#' 
#' # removing missing data
#' eta(c(x1, 2), c(NA, y), na.rm=TRUE)   # NA added to y to show NA behaviour
#'    
#' # classify interval data x
#' eta(x, y, breaks=c(1, 4, 7,10))
#' # visualize classication
#' plot(x, y)
#' abline(v=c(1, 4, 7,10))
#'  
#' # setting number of breaks for classification
#' eta(x, y, breaks=7)   
#'
eta <- function(x, y, breaks=NULL, na.rm=FALSE)
{
  if (is.factor(x))                           # convert factor to numeric in case it is a factor
    x <- as.numeric(x)                          
  if (! (is.vector(x) & is.vector(y)) ) 
      stop("'x' must be vectors or a factor, 'y' must be a vector", call.=FALSE)
  if (na.rm) {
    i <- !(is.na(x) | is.na(y))                 # listwise deletion for NAs
    x <- x[i]
    y <- y[i]
  }
  if (!is.null(breaks))                       # if x is interval data, breaks can be specified
    x <- droplevels(cut(x, breaks=breaks))    # as order does not matter drop levels
  l <- split(y, x)  
  .eta(l)
}


# x: list with y values for each category
# code adapted from: http://stackoverflow.com/questions/3002638/eta-eta-squared-routines-in-r
#
.eta <- function(x) 
{
  y <- unlist(x)                    # list with values by category
  mg <- sapply(x, mean)             # group means
  ng <- sapply(x, length)           # group size
  mtot <- mean(y)                   # total mean
  ssb <- sum(ng * (mg - mtot)^2)    # SSb
  sst <- sum((y - mtot)^2)          # SSt
  sqrt(ssb/sst)
}


