
eisera <- function(used, available, scannf = TRUE, nf = 2)
{
    ## Verifications
    if (!all(dim(used)==dim(available)))
        stop("used and available should have the same dimension")
    ut <- as.matrix(used)
    av <- as.matrix(available)

    ## Computation of the table of selection ratios
    av <- av/apply(av,1,sum)
    wij <- ut/apply(ut,1,sum)/av - 1

    ## If 0 availability, therefore wij = 0
    wij[av<1e-07] <- 0

    ## the table to be analysed
    mT <- sqrt(av)*wij
    ## row weights
    D <- apply(ut,1,sum)

    ## The eigenanalysis of selection ratios
    o <- as.dudi(as.data.frame(mT), rep(1,ncol(ut)), D,
                 scannf, nf, call=match.call(), type="esr")

    ## Output

    ## Scores of the habitat types
    uuu <- wij
    uuv <- apply(uuu,2,function(x) x*o$lw)
    o$co <- t(as.matrix(uuv))%*%as.matrix(o$l1)
    o$c1 <- NULL

    ## original tables
    o$available <- as.data.frame(av)
    o$used <- as.data.frame(ut)

    ## selection ratios
    o$wij <- ut/apply(ut,1,sum)/av

    return(o)
  }


print.esr <- function (x, ...)
{
  cat("Factorial analysis of selection ratios\n")
  cat("\n$call: ")
  print(x$call)
  cat("\n$nf:", x$nf, "axis-components saved")
  cat("\n$rank: ")
  cat(x$rank)
  cat("\neigen values: ")
  l0 <- length(x$eig)
  cat(signif(x$eig, 4)[1:(min(5, l0))])
  if (l0 > 5)
    cat(" ...\n")
  else cat("\n")
  sumry <- array("", c(3, 4), list(1:3, c("vector", "length",
                                          "mode", "content")))
  sumry[1, ] <- c("$cw", length(x$cw), mode(x$cw), "column weights")
  sumry[2, ] <- c("$lw", length(x$lw), mode(x$lw), "row weights")
  sumry[3, ] <- c("$eig", length(x$eig), mode(x$eig), "eigen values")
  class(sumry) <- "table"
  print(sumry)
  cat("\n")
  sumry <- array("", c(6, 4), list(1:6, c("data.frame", "nrow",
                                          "ncol", "content")))
  sumry[1, ] <- c("$tab", nrow(x$tab), ncol(x$tab), "modified array")
  sumry[2, ] <- c("$li", nrow(x$li), ncol(x$li), "row coordinates")
  sumry[3, ] <- c("$co", nrow(x$co), ncol(x$co), "column coordinates")
  sumry[4, ] <- c("$available", nrow(x$available), ncol(x$available),
                  "available proportions")
  sumry[5, ] <- c("$used", nrow(x$used), ncol(x$used), "number of relocations")
  sumry[6, ] <- c("$wij", nrow(x$used), ncol(x$used), "selection ratios")

  class(sumry) <- "table"
  print(sumry)
}




scatter.esr <- function(x, xax = 1, yax = 2, csub = 1,
                        possub = "bottomleft", ...)
  {
      ## Verifications
      if (!inherits(x, "esr"))
          stop("x should be of class \"esr\"")

      opar <- par(mfrow=c(2,1), mar=c(0,0,0,0))
      s.label(x$co, xax = xax, yax = yax, ...)
      if (csub > 0)
          scatterutil.sub("Habitat types", csub, possub)
      s.arrow(x$li, xax = xax, yax = yax, ...)
      if (csub > 0)
          scatterutil.sub("Animals", csub, possub)
      par(opar)
  }

