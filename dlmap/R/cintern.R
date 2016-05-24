cintern <- function(cc, tau, vrb, sigma2){
   ccnams <- names(tau)
   nc <- length(tau)
   cc <- lapply(cc, function(el, ccnams)
                {
                  if(is.numeric(el[[1]])) {
                    if(max(el[[1]]) > length(ccnams))
                      stop("coefficient subscript out of of bounds")
                       names(el[[1]]) <- ccnams[el[[1]]]
                  }
                  else {
                    if(any(pmatch(el[[1]], ccnams, 0) == 0))
                      stop("Names of contrast do not match the names of coefficients of object")
                    temp <- pmatch(el[[1]], ccnams)
                    names(temp) <- el[[1]]
                    el[[1]] <- temp
                  }
                  el
                }
                , ccnams)	## split contrasts and other available tests
   cons <- cc[as.logical(pmatch(unlist(lapply(cc, function(el)
                                                 el[[2]])), "con", 0, duplicates.ok = TRUE))]
   zero <- cc[!as.logical(pmatch(unlist(lapply(cc, function(el)
                                                  el[[2]])), "con", 0, duplicates.ok = TRUE))]
   cse <- ctau <- zwtest <- cwtest <- zpval <- c()
   if(length(cons)) {
        CRows <- lapply(cons, function(el, nc)
                        {
                          if(length(el) < 3){
                            con <- contr.helmert(length(el[[1]]))[, (length(el[[1]]) - 1)]
                            names(con) <- cnam <- names(el[[1]])
                            cat("Warning: default contrast being taken for", cnam, "is", con, "\n")
                            row <- rep(0, nc)
                            row[el[[1]]] <- con
                            row
                          }
                          else {
                            if(is.matrix(el[[3]])) {
                              cons <- split(el[[3]], 1:nrow(el[[3]]))
                              rows <- lapply(cons, function(ell, first = el[[1]], nc)
                                             {
                                               row <- rep(0, nc)
                                               row[first] <- ell
                                               row
                                             }
                                             , first = el[[1]], nc)
                              rows <- unlist(rows, use.names = FALSE)
                              matrix(rows, nrow = nrow(el[[3]]), byrow = TRUE)
                            }
                            else {
                              row <- rep(0, nc)
                              row[el[[1]]] <- el[[3]]
                              row
                            }
                          }
                        }
                        , nc)
        numrow <- lapply(cons, function(el)
                         {
                           if(length(el) > 2 && is.matrix(el[[3]]))
                             nrow(el[[3]])
                           else 1
                         })
    
        ## next line is temp line to sort out matrix problems 
        if(is.matrix(CRows[[1]]))
          CRows[[1]] <- t(CRows[[1]])				
        Cmat <- matrix(unlist(CRows), nrow = sum(unlist(numrow)), byrow = TRUE)
        for(i in 1:nrow(Cmat)) {
          varmat <- sum(Cmat[i,  ]*crossprod(vrb, t(Cmat)[, i]))
          cse[i] <- sqrt(varmat * sigma2)
          ctau[i] <- sum(Cmat[i,  ]*tau)
         cwtest[i] <- (ctau[i]/cse[i])^2
        }
        cres <- list(Cmat = Cmat, cwald = cwtest, cse = cse, ctau = ctau, cpval = 1 - pchisq(cwtest, 1))
      }
      if(length(zero)) {
        ZRows <- lapply(zero, function(el, nc)
                        {
                          rows <- rep(rep(0, nc), length(el[[1]]))
                          dum <- seq(0, (length(el[[1]]) - 1) * nc, by = nc)
                          rows[el[[1]] + dum] <- 1
                          matrix(rows, nrow = length(el[[1]]), byrow = TRUE)
                        }
                        , nc)
        for(i in 1:length(ZRows)) {
          varmat <- ZRows[[i]] %*% crossprod(vrb, t(ZRows[[i]]))
          Ctau <- ZRows[[i]] %*% tau
          zwtest[i] <- sum(Ctau*crossprod(solve(varmat), Ctau))/sigma2
          zpval[i] <- 1 - pchisq(zwtest[i], nrow(ZRows[[i]]))
        }
        zres <- list(ZRows = ZRows, zwald = zwtest, zpval = zpval)
      }
      result <- list(tau = tau, vrb = vrb, sigma2 = sigma2, cc = cc, call = call, zres = if(exists("zres")) zres
      else NULL, cres = if(exists("cres")) cres else NULL)
     result
}
