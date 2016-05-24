
hegy.boot.pval <- function(x, model0, stats0, 
  deterministic = c(1,0,0), lag.method = c("fixed", "AIC", "BIC"), maxlag = 0, 
  byseason = FALSE, nb = 500, u = NULL, debug.tid = -1)
{
  lag.method <- match.arg(lag.method)
  S <- frequency(x)
  #nx <- length(x) - S

  e <- residuals(model0)

  arcoefs <- coef(model0)
  arcoefs <- arcoefs[grepl("Lag\\d{1,3}$", names(arcoefs))]
  narcoefs <- length(arcoefs)
  Spnlags <- S + narcoefs
  ns <- c(length(x) + S + narcoefs, length(e))

  if (!is.null(u))
  {
    if (nrow(u) != ns[1])
      stop("not enough rows in argument ", sQuote("u"))
    if (ncol(u) < nb)
      stop("not enough columns in argument ", sQuote("u"))
  }

  refnms <- names(stats0)
  bpvals <- rep(0, length(stats0))

  for (i in seq_len(nb))
  {
    # generate data
    # resample residuals and use them as innovations in the model under 
    # the null hypothesis (seasonal random walk possibly with serial correlation)

    if (byseason)
    {
      #e <- ts(e, frequency = S, start = time(x)[S+narcoefs+1])
      e <- ts(e, frequency = S, start = time(x)[narcoefs+1])
      le <- split(e, cycle(e))

      if (is.null(u))
      {
        if ((length(e) %% S) == 0)
        {
          tmp <- lapply(le, FUN = function(x) 
            sample(x, size = length(x) + 2, replace = TRUE))          
          be <- do.call("rbind", tmp)
          be <- c(be[cycle(x)[seq_len(S)],])

        } else { # (length(e) %% S) != 0)
          tmp <- vector("list", S)
          nse <- table(cycle(x)) + 1
          j <- 0
          for (k in cycle(x)[seq_len(S)+narcoefs])
          {
            j <- j + 1
            tmp[[j]] <- sample(le[[k]], size = nse[k], replace = TRUE)
          }
          be <- suppressWarnings(c(do.call("rbind", tmp)))
          be <- be[seq_len(length(be) - length(x) %% S)]
        }

      } else { #!is.null(u)
        be <- unlist(le)[u[,i]]
      }

      if (narcoefs == 0)
      {
        bx <- ts(filter(be[-seq_len(S)], filter = c(rep(0, S-1),1), 
          method = "rec", init = rev(be[seq_len(S)])), 
          frequency = S, start = start(x))
      } else { # narcoefs > 0
        m <- outer(c(1, -arcoefs), c(1, rep(0, S-1),-1))
        pcoefs <- as.vector(tapply(m, row(m) + col(m), sum))
        bx <- ts(filter(be[-seq_len(Spnlags)], filter = -pcoefs[-1], 
          method = "rec", init = rev(be[seq_len(Spnlags)])), 
          frequency = S, start = start(x))
      }

    } else # !byseason
    {
      if (is.null(u)) {
        be <- sample(e, size = ns[1], replace = TRUE)
      } else {
        be <- e[u[,i]]
      }
      if (narcoefs == 0)
      {
        bx <- ts(filter(be[-seq_len(S)], filter = c(rep(0, S-1),1), 
          method = "rec", init = rev(be[seq_len(S)])), 
          frequency = S, start = start(x))
      } else # narcoefs > 0
      {
        m <- outer(c(1, -arcoefs), c(1, rep(0, S-1),-1))
        pcoefs <- as.vector(tapply(m, row(m) + col(m), sum))      
        bx <- ts(filter(be[-seq_len(Spnlags)], filter = -pcoefs[-1], 
          method = "rec", init = rev(be[seq_len(Spnlags)])), 
          frequency = S, start = start(x))
      }
    }

    # HEGY statistics for the bootstrap replicate

    bres <- hegy.test(bx, deterministic = deterministic, 
      lag.method = lag.method, maxlag = maxlag, pvalue = "raw")

    # counter
    # increment those statistics that are more extreme (in the right tail) 
    # than the statistics obtained for the original data 

    for (j in seq_along(bpvals))
    {
      if (grepl("t\\_", refnms[j]))
      {
        if (bres$stat[j] <= stats0[j])
          bpvals[j] <- bpvals[j] + 1
      } else 
        if (bres$stat[j] >= stats0[j])
          bpvals[j] <- bpvals[j] + 1
    }

    if (debug.tid == i-1)
      return(bx)

  } # loop seq_len(nb)

  names(bpvals) <- refnms  
  bpvals/nb
}
