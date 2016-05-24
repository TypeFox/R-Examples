"mlsolve" <-
function(ss, coef0, site.sel = "all", bruteforce = FALSE) {

  # Assume first column in ss is a site identifier
  names0 <- names(ss)
  site <- names0[1]
  imatch <- match(site, names0)
  names.all <- names(ss)[-imatch]
  ntaxa <- length(names.all)

  nvar <- length(coef0$xvar)
  tnames <- coef0$tnames
  csave <- coef0$csave

  m <- formtomat(as.character(coef0$form)[2], coef0$xvar)

  # if ss is based on otu, all taxon names in ss should be available
  # in tnames.  Check this first
  c <- matrix(NA, ncol = dim(coef0$csave)[[2]], nrow = ntaxa)
  for (i in 1:ntaxa) {
    imatch <- match(names.all[i], tnames)
    if (is.na(imatch)) {
      stop(paste(names.all[i],"not found."))
    }
    else {
      c[i,] <- csave[imatch,]
    }
  }
  
  if (site.sel == "all") {
    mat1 <- as.matrix(ss[, names.all])
    mat2 <- as.numeric(mat1 > 0)
    mat2[is.na(mat2)] <- 0
    dim(mat2) <- dim(mat1)
  }
  else {
    incvec <- ss[, site] == site.sel
    if (sum(incvec) > 0) {
      ss <- ss[incvec,]
      mat1 <- as.matrix(ss[, names.all])
      mat2 <- as.numeric(mat1 > 0)
      mat2[is.na(mat2)] <- 0
      dim(mat2) <- dim(mat1)
    }
    else{
      stop("Site selected was not found.")
    }
  }

  # define starting points for iterative solution search
  xs <- c(0.25,0.75)
  xsl <- as.list(rep(NA, times = nvar))
  for (i in 1:nvar) {
    xsl[[i]] <- xs
  }
  xstart.d <- expand.grid(xsl)

  ft <- function(x,ind) {
    xnew <- rep(NA, times = length(ind))
    xnew[1] <- 1
    for (i in 2:length(ind)) {
      xnew[i] <- prod(x[ind[[i]]])
    }
    return(xnew)
  }

  ftd <- function(x, ind.n, ind.p) {
    xnew <- rep(0, times = length(ind.n))
    for (i in 1:length(ind.n)) {
      xnew[i] <- prod(ind.n[i], x[ind.p[[i]]], na.rm = TRUE)
    }
    return(xnew)
  }

  loglik <- function(x,c, y) {
    xnew <- ft(x, m$ind)
    g <- c %*% xnew
    L <- sum(y*g + log(1/(1+exp(g))))
    return(L)
  }

  loglik.g <- function(x,c, y) {
    xnew <- ft(x, m$ind)
    g <- c %*% xnew
    dL <- rep(NA, times = nvar)
    dg <- matrix(NA, ncol = nvar, nrow = nrow(c))
    for (i in 1:nvar) {
      xnew2 <- ftd(x, m$ind.n[[i]], m$ind.p[[i]])
      dg[,i] <- c %*% xnew2
    }
    pg <- exp(g)/(1+exp(g))
    for (i in 1:nvar) {
      dL[i] <- sum(y*dg[,i]) - sum(pg*dg[,i])
    }
    return(dL)
  }
  
  infmat <- matrix(NA, nrow = nrow(mat2), ncol = nvar)
  flaginc <- rep(FALSE, times = nrow(mat2))

#  if (bruteforce) {
#    windows(height = 5, width = 5, pointsize = 10, restoreConsole = TRUE)
#    par(mfrow = c(1,1), pty = "s")
#  }

  for (j in 1:nrow(mat2)) {

    if (bruteforce) {
#      cat("Taxa observed at the site:\n")
#      print(names.all[as.logical(mat2[j,])])
#      cat("\n")
      xnew <- as.list(rep(NA, times = nvar))
      np <- 101
      for (i in 1:nvar) {
        xnew[[i]] <- seq(from = 0, to = 1,length = np)
      }
      new.data <- as.matrix(expand.grid(xnew))
      z <- rep(NA, times = nrow(new.data))
      for (i in 1:nrow(new.data)) {
        x <- drop(new.data[i,])
        z[i] <- loglik(x, c, mat2[j,])
      }
      if (nvar == 2) {
        dim(z) <- c(np,np)
      }
      if (is.factor(ss[, site])) {
        titlestr <- levels(ss[,site])[ss[,site]]
      }
      else {
        titlestr <- ss[,site]
      }
      if (nvar == 2) {
        par(mar = c(4,4,2,1))
        contour(xnew[[1]], xnew[[2]], z, nlevels = 10, main = paste(titlestr),
                axes = FALSE)
        box(bty = "o")
        at0 <- seq(from = 0, to = 1, length = 6)
        lab01 <- round(at0*diff(coef0$xlims[[1]]) + coef0$xlims[[1]][1],
                      digits = 1)
        lab02 <- round(at0*diff(coef0$xlims[[2]]) + coef0$xlims[[2]][1],
                      digits = 1)
        axis(1, at = at0, labels = lab01)
        axis(2, at= at0, labels = lab02)
        mtext(coef0$xvar[1], side = 1, line = 2.3)
        mtext(coef0$xvar[2], side = 2, line = 2.3)

      }
      else {
        par(mar = c(4,4,2,1))
        plot(xnew[[1]], z, type = "l", main = paste(titlestr), xlab = "",
             ylab = "", axes = FALSE)
        at0 <- seq(from = 0, to = 1, length = 6)
        lab0 <- round(at0*diff(coef0$xlims[[1]]) + coef0$xlims[[1]][1],
                      digits = 1)
        axis(1, at = at0, labels = lab0)
        axis(2)
        box(bty = "l")
        mtext(coef0$xvar[1], side = 1, line = 2.3)
        mtext("Log-likelihood", side = 2, line = 2.3)
      }
    }

    xgsave <- matrix(NA, ncol = nvar, nrow = nrow(xstart.d))
    liksav <- rep(NA, times = nrow(xstart.d))
    for (kk in 1:nrow(xstart.d)) {
      xg <- drop(xstart.d[kk,])
      o.out <- optim(xg, loglik, gr = loglik.g, method = "L-BFGS-B",
                     lower = rep(0, times = nvar),upper = rep(1, times = nvar),
                     c=c, y=mat2[j,],
                     control = list(fnscale = -1))
      xgsave[kk,] <- o.out$par
      liksav[kk] <- o.out$value
    }

    # Just check one dimension for consistency of results
    maxv <- apply(xgsave,2, max)
    minv <- apply(xgsave,2, min)

    incvec <- match(max(liksav), liksav)
    infmat[j,] <- xgsave[incvec,]

    difflike <- max(liksav) - min(liksav)
    if (sqrt(sum((maxv - minv)^2)) > 0.1) {
      if (difflike < 2) {
        if (bruteforce) {
          cat("\n")
          cat("Similar local maximums found.\n")
          cat("Likelihood values\n")
          print(liksav)
          cat("Point locations\n")
          print(xgsave)
        }
        flaginc[j] <- TRUE
      }
    }

    if (bruteforce) {
      if (nvar == 2) {
        points(infmat[j,1], infmat[j,2], pch = 16, cex = 1.2)
      }
      else {
        points(infmat[j,1], max(z), pch = 16, cex = 1.2)
      }
    }

  }


  inf2 <- infmat
  for (i in 1:nvar) {
    inf2[,i] <- infmat[,i]*diff(range(coef0$xlims[[i]])) + coef0$xlims[[i]][1]
  }

  dfout <- data.frame(ss[,site], inf2, flaginc)
  names(dfout) <- c(site, coef0$xvar, "Inconsistent")
  return(dfout)
}

