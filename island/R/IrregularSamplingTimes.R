# This is the script for the functions and generics for the analyses for regular
# sampling times.

# I added dtlsolver, dtlsolver2, dtlsolverdavid, incounts, levelwrap,
# levelwrap2, levelwrap2m, levelwrapb, times, wrapper, wrapper2, wrapperdavid,
# wrapperlike.

#' c/e rates for irregular samplings in a dataset
#'
#' \code{irregular_single_dataset} calculates colonization and extinction rates
#' in a single dataset.
#'
#'
#' @param dataframe A single dataframe.
#' @param vector A vector indicating the columns with presence-absence data.
#' @param c Tentative colonization rate.
#' @param e Tentative extinction rate.
#' @param column The name of the column with groups to calculate their c_e pair.
#' @param n Minimal number of rows for each group
#' @param int Accuracy to calculate the c_e pairs with.
#' @param assembly Logical indicating if the assembly starts from zero species or not.
#' @export
#' @seealso \code{\link{regular_sampling_scheme}},
#'   \code{\link{irregular_multiple_datasets}}
#' @note The columns with the presence-absence data should have the day of that
#'   sampling on the name of the column in order to calculate colonization and
#'   extinction.
#' @examples irregular_single_dataset(simberloff[[1]], 3:17, 0.001, 0.001)
#' irregular_single_dataset(simberloff[[1]], 3:17, column = "Tax. Unit 1",
#' 0.001, 0.001, 3)
#' \dontrun{
#' irregular_single_dataset(simberloff[[1]], 3:17, column = "Tax. Unit 1",
#' 0.001, 0.001, 3, 0.000001)
#' }
#' @return A dataframe with colonization and extinction rates and their upper
#'   and lower confidence interval, and if needed, the names of the groups to
#'   which colonization and extinction rates have been calculated.
irregular_single_dataset <- function(dataframe, vector, c, e, column = NULL,
                                     n = NULL, int = NULL, assembly = F) {

  if (is.null(column)) {
    wrapper(dataframe, vector, c, e, assembly)
  } else {
    if (is.null(int)) {
      levelwrap2(dataframe, vector, column, c, e, n, assembly)
    } else {
      wrapperlike(dataframe, vector, column, c, e, n, int, assembly)
    }
  }
}

#' c/e rates for irregular samplings in multiple datasets
#'
#' \code{irregular_multiple_datasets} calculates colonization and extinction
#' rates for data in several datasets.
#' @param list A list of dataframes.
#' @param vectorlist A list of vectors indicating the columns with
#'   presence-absence data.
#' @param c Tentative colonization rate.
#' @param e Tentative extinction rate.
#' @param column The name of the column with groups to calculate their c_e pair.
#' @param n Minimal number of rows for each group.
#' @param CI Logical. If TRUE, gives the confidence interval of the colonization
#'   and extinction rates.
#' @param assembly Logical indicating if the assembly starts from zero species or not.
#' @export
#' @seealso \code{\link{regular_sampling_scheme}},
#'   \code{\link{irregular_single_dataset}}
#' @note The columns with the presence-absence data should have the day of that
#'   sampling on the name of the column in order to calculate colonization and
#'   extinction.
#' @examples irregular_multiple_datasets(simberloff, list(3:17, 3:18, 3:17,
#' 3:19, 3:17, 3:16), 0.001, 0.001)
#' \dontrun{
#' irregular_multiple_datasets(simberloff, list(3:17, 3:18, 3:17, 3:19, 3:17,
#'  3:16), 0.001, 0.001, "Tax. Unit 1", n = 13)
#' irregular_multiple_datasets(simberloff, list(3:17, 3:18, 3:17, 3:19, 3:17,
#'  3:16), 0.001, 0.001, "Tax. Unit 1", n = 13, CI = TRUE)
#'  }
#' @return A dataframe with colonization and extinction rates and their upper
#'   and lower confidence interval, and if needed, the names of the groups to
#'   which colonization and extinction rates have been calculated.

irregular_multiple_datasets <- function(list,vectorlist, c, e, column=NULL,
                                        n=NULL, CI=FALSE, assembly = F) {

  if (is.null(column)) {
    wrapper2(list, vectorlist, c, e, assembly)
  } else {
    if (CI == F) {
      levelwrapb(list, vectorlist, column, c, e, n, assembly)
    } else {
      levelwrap(list, vectorlist, column, c, e, n, assembly)
    }
  }
}

dtlsolver <- function(x, t, n00, n10, n01, n11) {

  ### This function calculates the sum of the loglikelihood for a set x
  ### with two components, c and e, given a time series t and the number of
  ### transitions n00, n10, n01 and n11 in each time. Necessary for
  ### the solver.

  c   <- x[1]
  e   <- x[2]
  sum <- 0

  for (i in 1:length(t)) {
    ll <- n00[i] * log(1 - (c / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n11[i] * log(1 - (e / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n10[i] * log( (c / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n01[i] * log( (e / (e + c)) * (1 - exp( - (e + c) * t[i])))
    sum <- sum + ll
  }
  sum
}

dtlsolver2 <- function(x, t, abs, col, ext, pre) {

  ### Another version of the dtlsolver. Posibly, made to work with multiple
  ### discontinous times.

  c   <- x[1]
  e   <- x[2]
  res <- 0


  for (i in 1:length(t)) {
    sum   <- 0
    times <- t[[i]]
    n00   <- abs[[i]]
    n10   <- col[[i]]
    n01   <- ext[[i]]
    n11   <- pre[[i]]
    for (j in 1:length(times)) {
      ll <- n00[j] * log(1 - (c / (e + c)) * (1 - exp( - (e + c) * times[j]))) +
            n11[j] * log(1 - (e / (e + c)) * (1 - exp( - (e + c) * times[j]))) +
            n01[j] * log( (e / (e + c)) * (1 - exp( - (e + c) * times[j]))) +
            n10[j] * log( (c / (e + c)) * (1 - exp( - (e + c) * times[j])))
      sum <- sum + ll
    }
    res <- res + sum
  }
  res
}

dtlsolverdavid <- function(x, t, n00, n10, n01, n11) {

  ### Dtlsolver made with suggestions of David. I think it doesn't work well.

  c   <- x[1]
  e   <- x[2]
  sum <- 0

  for (i in 1:length(t)) {
    a00 <- 0
    a   <- n00[i]

    while (a > 0) {
      a00 <- a00 + log(1 - (c / (e + c)) * (1 - exp( - (e + c) * t[i])))
      a <- a - 1
    }
    a01 <- 0
    b   <- n01[i]
    while (b > 0) {
      a01 <- a01 + log( (e / (e + c)) * (1 - exp( - (e + c) * t[i])))
      b <- b - 1
    }
    a10 <- 0
    d   <- n10[i]
    while (d > 0) {
      a10 <- a10 + log( (c / (e + c)) * (1 - exp( - (e + c) * t[i])))
      d <- d - 1
    }
    a11 <- 0
    f   <- n11[i]
    while (f > 0) {
      a11 <- a11 + log(1 - (e / (e + c)) * (1 - exp( - (e + c) * t[i])))
      f <- f - 1
    }
    ll  <- a00 + a01 + a10 + a11
    sum <- sum + ll
  }
  sum
}

incounts <- function(dataset, vector, assembly = FALSE) {
  x <- dataset
  resultado <- data.frame()
  ### This function returns, from a dataset with data of presence-absence
  ### in each column of vector and the next in the dataset, a dataframe
  ### with four columns, corresponding to N01, N10, N00 and N11 respectively,
  ### and each row corresponding with a different temporal change. This
  ### function works for the cases in which the first column of the vector
  ### matches with the presence-absence at time 0 and those times when the
  ### first column is not at time 0, if the names of the columns are numeric.
  if (assembly) {

    if (as.numeric(colnames(x)[vector])[1] != 0) {
      resultado <- c(0, colSums(dataset[vector[1]])[[1]], (nrow(x) -
                                          colSums(dataset[vector[1]])[[1]]), 0)
    }
    for (j in vector) {
      if (j == utils::tail(vector, 1)) break
      N10 <- 0; N01 <- 0; N00 <- 0; N11 <- 0
      for (i in 1:nrow(x)) {

        if (x[i, j]  < x[i, j + 1]) N10 <- N10 + 1
        if (x[i, j]  > x[i, j + 1]) N01 <- N01 + 1
        if (x[i, j] == x[i, j + 1] && x[i, j] == 0) N00 <- N00 + 1
        if (x[i, j] == x[i, j + 1] && x[i, j] == 1) N11 <- N11 + 1


      }
      resultado <- rbind(resultado, c(N01, N10, N00, N11))
    }
    resultado
  } else {
    for (j in vector) {
      if (j == utils::tail(vector, 1)) break
      N10 <- 0; N01 <- 0; N00 <- 0; N11 <- 0
      for (i in 1:nrow(x)) {

        if (x[i, j]  < x[i, j + 1]) N10 <- N10 + 1
        if (x[i, j]  > x[i, j + 1]) N01 <- N01 + 1
        if (x[i, j] == x[i, j + 1] && x[i, j] == 0) N00 <- N00 + 1
        if (x[i, j] == x[i, j + 1] && x[i, j] == 1) N11 <- N11 + 1


      }
      resultado <- rbind(resultado, c(N01, N10, N00, N11))
    }
    resultado
  }
}

index <- function(list, column) {
  out <- vector()
  for (i in 1:length(list)) {
    x <- list[[i]]
    index <- levels (factor(x[, column]))
    out <- c(out, index)
  }
  out <- levels(factor(out))
  out
}

levelwrap <- function(list, vectorlist, column, c, e, n, assembly) {

  ### This function wraps the process of obtaining c and e in multiple
  ### datasets with different times (with different intervals) for the
  ### different levels contained in the specified "column". It also gives the
  ### C.I. calculated with the hessian.

  groups <- index(list, column)
  out2 <- data.frame()

  for (k in 1:length(groups)) {
    min  <- 0
    time <- vector("list", length(list))
    abs  <- vector("list", length(list))
    pre  <- vector("list", length(list))
    col  <- vector("list", length(list))
    ext  <- vector("list", length(list))

    for (i in 1:length(list)) {
      x <- list[[i]]
      dataset <- x[x[, column] == groups[k], ]
      if (nrow(dataset) == 0) next
      min <- min + nrow(dataset)
      vector <- vectorlist[[i]]
      co <- data.frame()
      co <- incounts(dataset, vector, assembly)
      time[i] <- data.frame(times(dataset, vector, assembly))
      abs[i]  <- data.frame(co[, 3])
      pre[i]  <- data.frame(co[, 4])
      col[i]  <- data.frame(co[, 2])
      ext[i]  <- data.frame(co[, 1])
    }

    if (min < n) next
    r <- stats::optim(c(c, e), dtlsolver2, t = time, abs = abs, col = col, ext =
                        ext, pre = pre, control = list(maxit = 5000000,
                                                    fnscale = - 1), hessian = T)



    fisher_info <- solve(-r$hessian)

    prop_sigma  <- sqrt(diag(fisher_info))
    prop_sigma  <- diag(prop_sigma)
    upper <- c(r$par[1] + 1.96 * prop_sigma[1, 1],
               r$par[2] + 1.96 * prop_sigma[2, 2])
    lower <- c(r$par[1] - 1.96 * prop_sigma[1, 1],
               r$par[2] - 1.96 * prop_sigma[2, 2])
    interval <- data.frame(value = r$par, upper = upper, lower = lower)
    out <- data.frame(groups[k], interval[1, ], interval[2, ])
    colnames(out) <- c("Group", "c", "Up-c", "Lo-c", "e", "Up-e", "Lo-e")
    out2 <- rbind(out2, out)

  }
  out2
}

levelwrap2 <- function(dataset, vector, column, c, e, n, assembly) {

  ### This function wraps the process of obtaining c and e in a dataset with
  ### different times (with different intervals) for the different levels
  ### contained in the specified "column". It also gives the C.I. calculated
  ### with the hessian.

  groups <- levels(factor(dataset[, column]))
  out2   <- data.frame()
  out    <- data.frame()

  for (k in 1:length(groups)) {
    x  <- dataset[dataset[, column] == groups[k], ]
    if (nrow(x) <= n) next
    co <- incounts(x, vector, assembly)
    t  <- times(x, vector, assembly)
    r  <- stats::optim(c(c, e), dtlsolver, t = t, n00 = co[, 3], n10 = co[, 2],
                n01 = co[, 1],n11 = co[, 4], control = list(maxit = 5000000,
                                                    fnscale = - 1), hessian = T)
    fisher_info <- solve(-r$hessian)
    prop_sigma  <- sqrt(diag(fisher_info))
    prop_sigma  <- diag(prop_sigma)
    upper <- c(r$par[1] + 1.96 * prop_sigma[1, 1],
               r$par[2] + 1.96 * prop_sigma[2, 2])
    lower <- c(r$par[1] - 1.96 * prop_sigma[1, 1],
               r$par[2] - 1.96 * prop_sigma[2, 2])
    interval <- data.frame(value = r$par, upper = upper, lower = lower)

    out <- data.frame(groups[k], interval[1, ], interval[2, ])
    colnames(out) <- c("Group", "c", "Up-c", "Lo-c", "e", "Up-e", "Lo-e")
    out2 <- rbind(out2, out)

  }
  out2
}

levelwrap2m <- function(dataset, vector, column, c, e, n) {

  ### This function wraps the process of obtaining c and e in a dataset with
  ### different times (with different intervals) for the different levels
  ### contained in the specified "column".

  groups <- levels(factor(dataset[, column]))
  out2   <- data.frame()
  out    <- data.frame()

  for (k in 1:length(groups)) {
    x <- dataset[dataset[, column] == groups[k], ]
    if (nrow(x) <= n) next
    co <- incounts(x, vector)
    t  <- times(x, vector)
    r  <- stats::optim(c(c, e), dtlsolver, t = t, n00 = co[, 3], n10 = co[, 2],
                n01 = co[, 1], n11 = co[, 4],
                control = list(maxit = 5000000, fnscale = - 1), hessian = F)

    out <- data.frame(groups[k], r$par[1], r$par[2])
    colnames(out) <- c("Group", "c", "e")
    out2 <- rbind(out2, out)

  }
  out2
}

levelwrapb <- function(list, vectorlist, column, c, e, n, assembly) {

  ### This function wraps the process of obtaining c and e in multiple
  ### datasets with different times (with different intervals) for the
  ### different levels contained in the specified "column".

  groups <- index(list, column)
  out2   <- data.frame()

  for (k in 1:length(groups)) {
    min  <- 0
    time <- vector("list", length(list))
    abs  <- vector("list", length(list))
    pre  <- vector("list", length(list))
    col  <- vector("list", length(list))
    ext  <- vector("list", length(list))
    for (i in 1:length(list)) {
      x <- list[[i]]
      dataset <- x[x[, column] == groups[k], ]
      if (nrow(dataset) == 0) next
      min <- min + nrow(dataset)
      vector <- vectorlist[[i]]
      co <- data.frame()
      co <- incounts(dataset, vector, assembly)
      time[i] <- data.frame(times(dataset, vector, assembly))
      abs[i]  <- data.frame(co[, 3])
      pre[i]  <- data.frame(co[, 4])
      col[i]  <- data.frame(co[, 2])
      ext[i]  <- data.frame(co[, 1])
    }

    if (min < n) next
    r <- stats::optim(c(c, e), dtlsolver2, t = time, abs = abs, col = col, ext =
                        ext, pre = pre, control = list(maxit = 5000000,
                                                       fnscale = - 1))

    out <- data.frame(groups[k], r$par[1], r$par[2])
    colnames(out) <- c("Group", "c", "e")
    out2 <- rbind(out2, out)

  }
  out2
}

times <- function(dataset, vector, assembly = FALSE) {

  ### Given a dataset, it calculates the changes in time from column to
  ### column in the vector. The names of the columns need to be numbers. The
  ### output is a vector of increments in time from column to column.
  if (assembly)  {
    times <- as.numeric(colnames(dataset)[vector])
    if (times[1] != 0) {out <- c(1:length(times))} else {out <- c(1:(length(times) - 1))}
    if (times[1] != 0) {out[1] <- times[1]} else {out[1] <- times[2] - times[1]}
    for (i in 2:length(out)) {
      if (times[1] == 0) {
        out[i] <- times[i + 1] - times[i]
      } else {
        out[i] <- times[i] - times[i - 1]
      }
    }
    out
  } else {
    times <- as.numeric(colnames(dataset)[vector])
    out <- c(1:(length(times) - 1))
    for (i in 1:length(out)) {

      out[i] <- times[i + 1] - times[i]
    }
    out
  }

}

wrapper <- function(dataset, vector, c, e, assembly) {

  ################## REVISAR LA DESCRIPCION
  ### This function wraps the whole process of obtaining c and e from a
  ### dataset with discontinuous times.

  co <- incounts(dataset, vector, assembly)
  t <- times(dataset, vector, assembly)
  r <- stats::optim(c(c, e), dtlsolver, t = t, n00 = co[, 3], n10 = co[, 2],
             n01 = co[, 1], n11 = co[, 4],
             control = list(maxit = 5000000, fnscale = - 1), hessian = T)
  fisher_info <- solve(-r$hessian)
  prop_sigma <- sqrt(diag(fisher_info))
  prop_sigma <- diag(prop_sigma)
  upper <- c(r$par[1] + 1.96 * prop_sigma[1, 1],
             r$par[2] + 1.96 * prop_sigma[2, 2])
  lower <- c(r$par[1] - 1.96 * prop_sigma[1, 1],
             r$par[2] - 1.96 * prop_sigma[2, 2])
  interval <- data.frame(value = r$par, upper = upper, lower = lower)
  out <- data.frame(interval[1, ], interval[2, ])
  colnames(out) <- c("c", "Up-c", "Lo-c", "e", "Up-e", "Lo-e")
  out
}

wrapper2 <- function(list, vectorlist, c, e, assembly) {

  ### This function wraps the whole process of obtaining c and e from a
  ### list of datasets with discontinuous times, that could be different for
  ### each dataset.

  time <- vector("list", length(list))
  abs  <- vector("list", length(list))
  pre  <- vector("list", length(list))
  col  <- vector("list", length(list))
  ext  <- vector("list", length(list))

  for (i in 1:length(list)) {
    dataset <- list[[i]]
    vector <- vectorlist[[i]]
    co <- data.frame()
    co <- incounts(dataset, vector, assembly)
    time[i] <- data.frame(times(dataset, vector, assembly))
    abs[i]  <- data.frame(co[, 3])
    pre[i]  <- data.frame(co[, 4])
    col[i]  <- data.frame(co[, 2])
    ext[i]  <- data.frame(co[, 1])
  }

  r <- stats::optim(c(c, e), dtlsolver2, t = time, abs = abs, col = col, ext =
                      ext, pre = pre, control = list(maxit = 5000000, fnscale =
                                                       - 1), hessian = T)

  fisher_info <- solve(-r$hessian)
  prop_sigma  <- sqrt(diag(fisher_info))
  prop_sigma  <- diag(prop_sigma)
  upper <- c(r$par[1] + 1.96 * prop_sigma[1, 1],
             r$par[2] + 1.96 * prop_sigma[2, 2])
  lower <- c(r$par[1] - 1.96 * prop_sigma[1, 1],
             r$par[2] - 1.96 * prop_sigma[2, 2])
  interval <- data.frame(value = r$par, upper = upper, lower = lower)
  out <- data.frame(interval[1, ], interval[2, ])
  colnames(out) <- c("c", "Up-c", "Lo-c", "e", "Up-e", "Lo-e")
  out
}

wrapperdavid <- function(dataset, vector, c, e) {

  ### I think that this one calculates the loglikelihood of the data for a
  ### model with discontinous times given a specific c and e.

  ############# REVISAR

  times <- as.numeric(colnames(dataset)[vector])
  if (times[1] != 0) {
    out <- c(1:length(times))
  } else {
    out <- c(1:(length(times) - 1))
  }
  if (times[1] != 0) {
    out[1] <- times[1]
  } else {
    out[1] <- times[2] - times[1]
  }

  for (i in 2:length(out)) {
    if (times[1]==0) {
      out[i] <- times[i + 1] - times[i]
    } else {
      out[i] <- times[i] - times[i - 1]
    }
  }
  t <- out
  x <- dataset
  resultado <- data.frame()
  if (as.numeric(colnames(x)[vector])[1] != 0) {
    resultado <- c(0, colSums(dataset[vector[1]])[[1]],
                   (nrow(x) - colSums(dataset[vector[1]])[[1]]), 0)
  }
  for (j in vector) {
    if (j == utils::tail(vector, 1)) break
    N10 <- 0; N01 <- 0; N00 <- 0; N11 <- 0
    for (i in 1:nrow(x)) {

      if (x[i, j]  < x[i, j + 1]) N10 <- N10 + 1
      if (x[i, j]  > x[i, j + 1]) N01 <- N01 + 1
      if (x[i, j] == x[i, j + 1] && x[i, j] == 0) N00 <- N00 + 1
      if (x[i, j] == x[i, j + 1] && x[i, j] == 1) N11 <- N11 + 1
    }

        resultado <- rbind(resultado, c(N01, N10, N00, N11))
  }
  n00 <- resultado[, 3]
  n10 <- resultado[, 2]
  n01 <- resultado[, 1]
  n11 <- resultado[, 4]
  sum <- 0

  for (i in 1:length(t)) {
    ll <- n00[i] * log(1 - (c / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n11[i] * log(1 - (e / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n10[i] * log( (c / (e + c)) * (1 - exp( - (e + c) * t[i]))) +
          n01[i] * log( (e / (e + c)) * (1 - exp( - (e + c) * t[i])))
    sum <- sum + ll
  }
  sum
}

wrapperlike <- function(dataset, vector, column, c, e, n, int, assembly) {

  ### This function wraps the whole process of obtaining c and e from a
  ### dataset with discontinuous times for the different groups contained in
  ### the "column" with more than "n" otus. It obtains the 95% C.I. of the
  ### rates by the method of profile likelihood, and an "int" interval can be
  ### specified to calculate the profile.

  groups <- levels(factor(dataset[, column]))
  out2 <- data.frame()
  out  <- data.frame()
  for (k in 1:length(groups)) {

    x <- dataset[dataset[, column] == groups[k], ]
    if (nrow(x) <= n) next
    co <- incounts(x, vector, assembly)
    t  <- times(x, vector, assembly)
    r  <- stats::optim(c(c, e), dtlsolver, t = t, n00 = co[, 3], n10 = co[, 2],
                n01 = co[, 1], n11 = co[, 4],
                control = list(maxit = 5000000, fnscale = - 1), hessian = F)

    c1 <- r$par[1]
    e1 <- r$par[2]


    llce <- r$value
    df <- data.frame(matrix(0, 2, 3))
    lli <- llce

    df[1, 1] <- c1
    df[2, 1] <- e1
    cup <- NULL
    clo <- NULL
    eup <- NULL
    elo <- NULL
    i <- 1
    while ( (llce - lli) < 2.1) {
      cup <- c1 + (int * i)
      lli <- dtlsolver(c(cup, e1), t, n00 = co[, 3], n10 = co[, 2],
                       n01 = co[, 1], n11 = co[, 4])
      i   <- i + 1
      if ( (llce - lli) < 2) next
      df[1, 3] <- cup
      break
    }
    i <- 1
    lli <- llce
    while ( (llce - lli) < 2.1) {
      clo <- c1 - (int * i)
      i <- i + 1
      lli <- dtlsolver(c(clo, e1), t, n00 = co[, 3], n10 = co[, 2],
                       n01 = co[, 1], n11 = co[, 4])
      if ( (llce - lli) < 2) next
      df[1, 2] <- clo
      break
    }
    i <- 1
    lli <- llce
    while ( (llce - lli) < 2.1) {
      eup <- e1 + (int * i)
      lli <- dtlsolver(c(c1, eup), t, n00 = co[, 3], n10 = co[, 2],
                       n01 = co[, 1], n11 = co[, 4])
      i <- i + 1
      if ( (llce - lli) < 2) next
      df[2, 3] <- eup
      break
    }
    i <- 1
    lli <- llce
    while ( (llce - lli) < 2.1) {
      elo <- e1 - (int * i)
      lli <- dtlsolver(c(c1, elo), t, n00 = co[, 3], n10 = co[, 2],
                       n01 = co[, 1], n11 = co[, 4])
      i <- i + 1
      if ( (llce - lli) < 2) next
      df[2, 2] <- elo
      break
    }

    out <- data.frame(groups[k], nrow(x), r$par[1], df[1, 2], df[1, 3],
                      r$par[2], df[2, 2], df[2, 3])
    colnames(out) <- c("Group", "N", "c", "clo", "cup", "e", "elo", "eup")
    out2 <- rbind(out2, out)

  }
  out2

}
