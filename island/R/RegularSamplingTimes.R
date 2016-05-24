# This is the script for the functions and generics for the analyses for regular
# sampling times.

# I added cetotrans, changes, cilike, cilprobs, ciplevel, ciplevel2,
# ciplevel2b, lfun, lfun2, lfunr, lfunrsolver, like, multirates, rates, rlevel,
# rlevel2, rlinterval, tplevel.

# Revisar cilike,

#' c/e rates for a regular sampling scheme
#'
#' \code{regular_sampling_scheme} calculates colonization and extinction rates
#' for a community or groups in a community.
#'
#' @param x A single dataset.
#' @param vector A vector indicating the columns with presence-absence data.
#' @param level The name of the column with groups to calculate their c_e pair.
#' @param otus Minimal number of rows for each group
#' @param int Accuracy to calculate the c_e pairs with.
#' @examples regular_sampling_scheme(alonso[[1]],3:6)
#' regular_sampling_scheme(alonso[[1]],3:6,"Guild",5)
#' regular_sampling_scheme(alonso[[1]],3:6,"Guild",5, 0.001)
#' @return A dataframe with colonization and extinction rates along with their
#'   associated transition probabilities or their lower and upper confidence
#'   intervals, for each group if specified.
#' @seealso \code{\link{irregular_single_dataset}},
#'   \code{\link{irregular_multiple_datasets}}
#' @export
regular_sampling_scheme <- function(x, vector, level = NULL, otus = NULL,
                                    int = NULL) {

  if (is.null(level)) {
    rates(changes(x,vector))
  } else {
    if (is.null(int)) {
      rlevel2(x, level, vector, otus)
    } else {
      rlinterval(x, level, vector, otus, int)
    }
  }
}


changes <- function(x, vector) {

  ### This function identifies and counts the kind of transition produced in a
  ### dataframe "x" among every column in "vector". The output is a dataframe
  ### with the number of extinctions, colonizations and mantained absences or
  ### presences in the dataframe.

  N01<-0
  N10<-0
  N00<-0
  N11<-0
  resultado<-data.frame()
  for (i in 1:nrow(x)){
    for (j in vector) {
      if (utils::tail(vector,1)==j) break
      if (x[i,j]<x[i,j+1]) N10<-N10+1
      if (x[i,j]>x[i,j+1]) N01<-N01+1
      if (x[i,j]==x[i,j+1] && x[i,j]==0) N00<-N00+1
      if (x[i,j]==x[i,j+1] && x[i,j]==1) N11<-N11+1
      if (i!=nrow(x)) next
      resultado<-rbind(c(N01,N10,N00,N11))
    }
  }
  resultado
}

lfun <- function(c, e, x) {

  ### This function calculates the likelihood of a given set of transitions
  ### x for the rates c and e.

  (1 - ( (c / (e+c)) * (1 - exp( - (e + c))))) ^ x[1, 3] * ( (c / (e + c)) *
  (1 - exp( - (e + c)))) ^ x[1, 2] * ( (e / (e + c)) * (1 - exp( - (e + c)))) ^
  x[1, 1] * (1 - ( (e / (e + c)) * (1 - exp( - (e + c))))) ^ x[1, 4]
}

lfun2 <- function(c, e, x) {

  ### This function calculates the log-likelihood of a given set of
  ### transitions x for the rates c and e.

  log(1 - ( (c / (e + c)) * (1 - exp( - (e + c))))) * x[1, 3] +
  log(1 - ( (e / (e + c)) * (1 - exp( - (e + c))))) * x[1, 4] +
  log( (c / (e + c)) * (1 - exp( - (e + c)))) * x[1, 2] +
  log( (e / (e + c)) * (1 - exp( - (e + c)))) * x[1, 1]

}

lfunr <- function(T01, T10, x) {

  ### This function calculates the log-likelihood of a given set of
  ### transitions x for the transition probabilities T01 and T10.

  log(1 - T10) * x[1, 3] + log(T10) * x[1, 2] + log(T01) * x[1, 1] +
  log(1 - T01) * x[1, 4]
}

lfunrsolver <- function(y, x) {

  ### This function calculates the loglikelihood of a given dataset of
  ### transitions x for the vector y of transition probabilities T01
  ### and T10.

  T01 <- y[1]
  T10 <- y[2]
  log(1 - T10) * x[1, 3] + log(T10) * x[1, 2] + log(T01) * x[1, 1] +
  log(1 - T01) * x[1, 4]
}

like <- function(x) {

  ### This function calculates the log-likelihood of a given set of
  ### transitions x for the c and e calculated with the transitions in
  ### x.

  out <- numeric()
  edivc1  <- x[1, 1] * (x[1, 2] + x[1, 3]) / (x[1, 2] * (x[1, 1] + x[1, 4]))
  eplusc1 <- - log( (x[1, 3] * x[1, 4] - x[1, 2] * x[1, 1]) /
                    ( (x[1, 2] + x[1, 3]) * (x[1, 1] + x[1, 4])))
  c <- eplusc1 / (edivc1 + 1)
  e <- eplusc1 / (1 + (1 / edivc1))
  T00 <- 1 - ( (c / (e + c)) * (1 - exp( - (e + c))))
  T11 <- 1 - ( (e / (e + c)) * (1 - exp( - (e + c))))
  T01 <- (e / (e + c)) * (1 - exp( - (e + c)))
  T10 <- (c / (e + c)) * (1 - exp( - (e + c)))
  out <- log(T00) * x[1, 3] + log(T10) * x[1, 2] + log(T01) * x[1, 1] +
         log(T11) * x[1, 4]
  out
}

rates <- function(x) {

  ### This function calculates the rates and probabilities (c, e, T00, T01,
  ### T11, T10) for the given set of transitions x.

  sol <- data.frame()
  edivc1 <- (x[1, 1] * (x[1, 2] + x[1, 3])) / (x[1, 2] * (x[1, 1] + x[1, 4]))
  eplusc1 <- - log( ( (x[1, 3] * x[1, 4]) - (x[1, 2] * x[1, 1])) /
                    ( (x[1, 2] + x[1, 3]) * (x[1, 1] + x[1, 4])))
  c <- eplusc1 / (edivc1 + 1)
  e <- eplusc1 / (1 + (1 / edivc1))
  T00 <- 1 - ( (c / (e + c)) * (1 - exp( - (e + c))))
  T01 <-       (e / (e + c)) * (1 - exp( - (e + c)))
  T11 <- 1 - ( (e / (e + c)) * (1 - exp( - (e + c))))
  T10 <-       (c / (e + c)) * (1 - exp( - (e + c)))

  sol <- cbind(c, e, T00, T01, T11, T10)
  colnames(sol) <- c("c", "e", "T00", "T01", "T11", "T10")
  sol

}

rlevel2 <- function(x, level, vector, otus) {

  ### Starting from a dataset x that contains transitions in the columns
  ### specified in vector, this function calculates the rates and
  ### probabilities (c, e, T00, T01, T11, T10)  for the groups in level
  ### that have more than a number otus of replicas.

  results <- data.frame()
  ## First, we make an index of the actual elements in "level".
  index <- levels(factor(x[, level]))
  ## Now, start the calculations.
  for (i in 1:length(index)) {
    y <- x[x[, level] == index[i], ]
    if (nrow(y) < otus) next
    rr <- rates(changes(y, vector))
    rr <- data.frame(index[i], nrow(y), rr)
    results <- rbind(results, rr)
  }
  colnames(results) <- c("Group", "N", "c", "e", "T00", "T01", "T11", "T10")
  results
}

rlinterval <- function(x, level, vector, otus, int) {

  ### Starting from a dataset x that contains transitions in the columns
  ### specified in vector, this function calculates the rates and its 95%
  ### confidence interval for the groups in level that have more than a
  ### number "otus" of replicas.

  out <- data.frame()
  out2 <- data.frame()
  ## First, we make an index of the actual elements in "level".
  index <- levels(factor(x[, level]))
  ## Now, start the calculations.
  for (i in 1:length(index)) {
    y <- x[x[, level] == index[i], ]
    if (nrow(y) < otus) next
    rr <- rates(changes(y, vector))
    c1 <- rr[1]
    e1 <- rr[2]
    if (is.na(c1) == T) next
    z <- changes(y, vector)

    llce <- lfun2(c1, e1, z)
    df <- data.frame(matrix(0, 1, 6))
    lli <- llce

    cup <- NULL
    clo <- NULL
    eup <- NULL
    elo <- NULL

    j <- 1
    while ( (llce-lli) < 2.1) {
      cup <- c1 + (int * j)
      lli <- lfun2(cup, e1, z)
      j <- j + 1
      if ( (llce - lli) < 2) next
      df[1, 3] <- cup
      break
    }
    j <- 1
    lli <- llce
    while ( (llce - lli) < 2.1) {
      clo <- c1 - (int * j)
      j <- j + 1
      lli <- lfun2(clo, e1, z)
      if ( (llce - lli) < 2) next
      df[1, 2] <- clo
      break
    }
    j <- 1
    lli <- llce
    while ( (llce - lli) < 2.1) {
      eup <- e1 + (int * j)
      lli <- lfun2(c1, eup, z)
      j <- j + 1
      if ( (llce - lli) < 2) next
      df[1, 6] <- eup
      break
    }
    j <- 1
    lli <- llce
    while ( (llce - lli) < 2.1) {
      elo <- e1 - (int * j)
      lli <- lfun2(c1, elo, z)
      j <- j + 1
      if ( (llce - lli) < 2) next
      df[1, 5] <- elo
      break
    }

    out <- data.frame(index[i], nrow(y), c1, df[1, 2], df[1, 3], e1, df[1, 5],
                      df[1, 6])
    colnames(out) <- c("Group", "N", "c", "clo", "cup", "e", "elo", "eup")
    out2 <- rbind(out2, out)


  }
  out2
}
