### File version $Id: mcgibbsit.R 55 2005-03-15 03:19:21Z r_burrows $

"mcgibbsit" <- function (data, q = 0.025, r = 0.0125, s = 0.95,
                         converge.eps = 0.001, correct.cor=TRUE )
{
  ## if asked for more than one q,r,s,..., combination, recurse and return
  ## a list of results

  parms <- cbind( q, r, s, converge.eps, correct.cor ) # so one can use
                                                      # different r/s values
                                                      # for each q

  if(dim(parms)[1] > 1 )
    {

      retval <- list();
      for(i in 1:length(q) )
        {
        retval[[i]] <- mcgibbsit( data, q=parms[i,"q"], r=parms[i,"r"],
                                 s=parms[1,"s"],
                                 converge.eps=parms[i,"converge.eps"],
                                 correct.cor=parms[i,"correct.cor"] )
      }
      return(retval)
    }

  if (is.mcmc.list(data))
    {
      nchains <- length(data)

      ## check that all of the chains are conformant #
      for( ch in 1:nchains)
        if( (start(data[[1]]) != start(data[[ch]] )) ||
           (end  (data[[1]]) != end  (data[[ch]] )) ||
           (thin (data[[1]]) != thin (data[[ch]] )) ||
           (nvar (data[[1]]) != nvar (data[[ch]] )) )
          stop(paste("All chains in mcmc.list must have same 'start',",
                     "'end', 'thin', and number of variables"));

      if(is.matrix(data[[1]]))
        combined <- mcmc(do.call("rbind",data))
      else
        combined <- mcmc(as.matrix(unlist(data)))
    }
  else
    {
      data <- mcmc.list(mcmc(as.matrix(data)))
      nchains <- 1
      combined <- mcmc(as.matrix(data))
    }


  resmatrix <- matrix(nrow = nvar(combined),
                      ncol = 6,
                      dimnames = list( varnames(data, allow.null = TRUE),
                                       c("M", "N", "Total", "Nmin",
                                         "I", "R" )) )

  # minimum number of iterations
  phi <- qnorm(0.5 * (1 + s))
  nmin <- as.integer(ceiling((q * (1 - q) * phi^2)/r^2))

  if (nmin > niter(combined))
    resmatrix <- c("Error", nmin)
  else
    for (i in 1:nvar(combined))
      {
        dichot <- list()

        if (is.matrix(data[[1]]))
          {
            quant <- quantile(combined[, i, drop = TRUE], probs = q)

            for(ch in 1:nchains)
              {
                dichot[[ch]] <- mcmc(data[[ch]][, i, drop = TRUE] <= quant,
                                     start = start(data[[ch]]),
                                     end   = end(data[[ch]]),
                                     thin  = thin(data[[ch]]))
              }
          }
        else
          {
            quant <- quantile(combined, probs = q)

            for(ch in 1:nchains)
              {
                dichot[[ch]] <- mcmc(data[[ch]] <= quant,
                                     start = start(data[[ch]]),
                                     end = end(data[[ch]]),
                                     thin = thin(data[[ch]]))
              }
          }
        kthin <- 0
        bic <- 1

        while (bic >= 0)
          {
            kthin <- kthin + thin(data[[1]])

            to.table <- function(dichot, kthin)
              {
                testres <- as.vector(window(dichot, thin = kthin))
                newdim <- length(testres)
                testtran <- table(testres[1:(newdim - 2)],
                                  testres[2:(newdim - 1)],
                                  testres[3:newdim])

                ## handle the case where one or more of the transition never
                ## happens, so that the testtran array has two few dimensions
                if( any(dim(testtran!=2)) )
                  {
                    tmp <- array( 0, dim=c(2,2,2),
                                dimnames=list( c("FALSE","TRUE"),
                                               c("FALSE","TRUE"),
                                               c("FALSE","TRUE") ) )
                    for(t1 in dimnames(testtran)[[1]] )
                      for(t2 in dimnames(testtran)[[2]])
                        for(t3 in dimnames(testtran)[[3]] )
                          tmp[t1,t2,t3] <- testtran[t1,t2,t3]

                    testtran <- tmp
                  }

                testtran <- array(as.double(testtran), dim = dim(testtran))
                return(testtran)
              }

            tmp <- sapply( dichot, to.table, kthin=kthin, simplify=FALSE )

            ## add all of the transition matrixes together
            testtran <- tmp[[1]]
            if(nchains>1)
              for(ch in 2:nchains)
                testtran <- testtran + tmp[[ch]]

            ## compute the likelihoood
            g2 <- 0
            for (i1 in 1:2)
              {
                for (i2 in 1:2)
                  {
                    for (i3 in 1:2)
                      {
                        if (testtran[i1, i2, i3] != 0) {
                          fitted <- (sum(testtran[i1, i2, 1:2]) *
                                     sum(testtran[1:2, i2, i3])) /
                                       (sum(testtran[1:2, i2, 1:2]))
                          g2 <- g2 + testtran[i1, i2, i3] *
                                     log(testtran[i1, i2, i3]/fitted) * 2
                        }
                      }
                  }
              }

            ## compute bic
            bic <- g2 - log( sum(testtran) - 2 ) * 2

          }

        ## estimate the parameters of the first-order markov chain
        alpha <- sum(testtran[1, 2, 1:2]) / (sum(testtran[1, 1, 1:2]) +
                                             sum(testtran[1, 2, 1:2]))
        beta  <- sum(testtran[2, 1, 1:2]) / (sum(testtran[2, 1, 1:2]) +
                                             sum(testtran[2, 2, 1:2]))

        ## compute burn in
        tempburn <- log((converge.eps * (alpha + beta)) /
                        max(alpha, beta))/(log(abs(1 - alpha - beta)))

        nburn <- as.integer(ceiling(tempburn) * kthin)

        ## compute iterations after burn in
        tempprec <- ((2 - alpha - beta) * alpha * beta * phi^2)/
                    (((alpha + beta)^3) * r^2)
        nkeep  <- ceiling(tempprec * kthin)

        ## compute the correlation
        if(nchains>1 && correct.cor)
          {
            dat <- do.call("cbind", dichot)
            varmat <- var(dat)
            denom <- mean(diag(varmat))  # overall variance
            diag(varmat) <- NA
            numer <- mean(c(varmat), na.rm=TRUE) # overall covariance
            rho <- numer / denom
          }
        else
          rho <- 1.0

        ## inflation factors
        iratio <- (nburn + nkeep)/nmin
        R      <- ( 1 + rho * (nchains - 1) )

        resmatrix[i, 1] <- M <- ceiling( nburn * nchains )  # M
        resmatrix[i, 2] <- N <- ceiling( nkeep * R       )  # N
        resmatrix[i, 3] <- M + N                            # Total
        resmatrix[i, 4] <- nmin                             # nmin
        resmatrix[i, 5] <- signif(iratio,  digits = 3)      # I
        if(nchains > 1 && correct.cor )
          resmatrix[i, 6] <- signif(R     ,  digits = 3)    # R
        else
          resmatrix[i, 6] <- NA                             # R

      }
  y <- list(params = c(r = r, s = s, q = q),
            resmatrix = resmatrix, call=match.call(),
            nchains = nchains, len=(end(data[[1]]) - start(data[[1]]) + 1) )
  class(y) <- "mcgibbsit"
  return(y)
}

"print.mcgibbsit" <-
  function (x, digits = 3, ...)
{
  cat("                  Multi-Chain Gibbsit \n")
  cat("                  ------------------- \n")
  cat("\n");

  cat("Call             = "); print(x$call)
  cat("\n");

  cat("Number of Chains =", x$nchains, "\n" )
  cat("Per-Chain Length =", x$len, "\n" )
  cat("Total Length     =", x$nchains * x$len, "\n")
  cat("\n");

  cat("Quantile (q)     =", x$params["q"] , "\n")
  cat("Accuracy (r)     = +/-", x$params["r"], "\n")
  cat("Probability (s)  =", x$params["s"], "\n")
  cat("\n")

  if (x$resmatrix[1] == "Error")
    cat("\nYou need a sample size of at least", x$resmatrix[2],
        "with these values of q, r and s\n")
  else {
    out <- x$resmatrix
    for (i in ncol(out)) out[, i] <- format(out[, i], digits = digits)

    maxM <- max(x$resmatrix[,1])
    maxN <- max(x$resmatrix[,2])

    out <- rbind(c("Burn-in ", "Estimation", "Total", "Lower bound ",
                   "Auto-Corr.", "Between-Chain"),
                 c("(M)",       "(N)",       "(M+N)", "(Nmin)",
                   "factor (I)", "Corr. factor (R)"),
                 rep('', ncol(out)),
                 out,
                 rep('-----', ncol(out) ),
                 c(maxM, maxN, maxM+maxN, "", "", "")
                 )


#    if (!is.null(rownames(x$resmatrix)) || all(rownames(x$resmatrix)==''))
#      out <- cbind(c("", "", rownames(x$resmatrix)), out)

    colnames(out) <- rep("", ncol(out))



    print.default(out, quote = FALSE, ...)
    cat("\n")

    cat("NOTE: The values for M, N, and Total are combined numbers",
        " of iterations \n")
    cat("      based on using", x$nchains, "chains.\n");
    cat("\n")

  }
  invisible(x)
}
