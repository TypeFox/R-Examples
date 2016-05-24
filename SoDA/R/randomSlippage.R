  randomSlippage  <- function(nRuns, expr1, expr2, slip = runif(1), check = FALSE) {
      slippage <- matrix(NA, nRuns, 2)
      if(!exists(".Random.seed", envir = .GlobalEnv))
        stop("The random number generator must be initiated before calling randomSlippage()")
      firstSeed <- .Random.seed
      e1 <- substitute(expr1); e2 <- substitute(expr2); eslip = substitute(slip)
      for(i in seq(length=nRuns)) {
           g1 = eval.parent(e1)
          saveSeed <- .Random.seed
          g21 = eval.parent(e2)
           .Random.seed <<- saveSeed
           u1 = eval.parent(eslip)
          g22 = eval.parent(e2)
          m <- match(g21, g22)
          seqn2 <- seq(along=g21)
          k <- seqn2[!is.na(m)]
          if(length(k) > 0) {
              k1 <- k[[1]]
              slippage[i,] <- c(k1, m[[k1]])
              n2 = length(g21)
              if(check && k1 <  n2  && ( any(diff(k) != 1) || any(diff(m[k]) != 1)))
                stop("Non-synchronized samples but with matching numbers!")
          }
        }
      attr(slippage, "seed") <- firstSeed
      attr(slippage, "expressions") <- list(e1, e2, eslip)
      slippage
  }
