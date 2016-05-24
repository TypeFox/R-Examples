"plotKinBreakDown" <-
function (multimodel, multitheta, plotoptions) 
{
  if(dev.cur() != 1)
    dev.new()
  plotrow <- length(plotoptions@breakdown$plot)
  plotcol <- 1
  par(plotoptions@paropt)
  par(oma = c(0, 0, 4, 0))
  par(mfrow = c(plotrow, plotcol))
  m <- multimodel@modellist
  t <- multitheta
  
  superimpose <- plotoptions@breakdown$superimpose
  if(length(superimpose) == 0 || any(superimpose > length(m)))
    superimpose <- 1:length(m)
  indList <- list()
  divdrel <- plotoptions@divdrel
  m <- multimodel@modellist
  t <- multitheta
  res <- multimodel@fit@resultlist
  allx2 <- allx <- vector()
  pl <- plotoptions@breakdown$plot
  if (length(plotoptions@breakdown$tol) == 0) 
    tol <- abs(m[[1]]@x2[1] - m[[1]]@x2[2])/100
  else tol <- plotoptions@breakdown$tol
  groups <- multimodel@groups
  grtoplot <- vector("list", length(m))
  cnt <- 1
  ds <- vector()
  while (length(ds) < length(m)) {
    for (i in 1:length(groups[[cnt]])) {
      gds <- groups[[cnt]][[i]][2]
      if (!gds %in% ds) {
        ds <- append(ds, gds)
        grtoplot[[gds]] <- list(groups[[cnt]], i)
      }
      ds <- append(ds, gds)
    }
  }
  for (j in 1:length(pl)) {
    indList[[j]] <- vector()
    for (i in superimpose) {
      if (length(which(m[[i]]@x2 >= pl[j])) > 0) 
        x1 <- min(which(m[[i]]@x2 >= pl[j]))
      else x1 <- NA
      if (length(which(m[[i]]@x2 <= pl[j])) > 0) 
        x2 <- max(which(m[[i]]@x2 <= pl[j]))
      else x2 <- NA
      if (is.na(x1)) 
        x3 <- x2
      if (is.na(x2)) 
        x3 <- x1
      if (length(na.omit(c(m[[i]]@x2[x1], m[[i]]@x2[x2]))) == 2) {
        if (abs(m[[i]]@x2[x1] - pl[j]) <= abs(m[[i]]@x2[x2] - 
                 pl[j])) 
          x3 <- x1
        else x3 <- x2
      }
      if ((pl[j] <= m[[i]]@x2[x3] + tol) && (pl[j] >= m[[i]]@x2[x3] - tol)) {
        indList[[j]] <- append(indList[[j]], x3)
      }
      else indList[[j]] <- append(indList[[j]], NA)
    }
  }
  for (i in superimpose) {
    allx2 <- append(allx2, m[[i]]@x2)
    allx <- append(allx, m[[i]]@x)
  }
  
  allx2 <- sort(unique(allx2))
  allx <- sort(unique(allx))
  xmax <- max(allx)
  xmin <- min(allx)
  k <- t[[1]]@kinpar
  conmax <- list()
  spectraList <- getSpecList(multimodel, t)
  for (b in 1:length(pl)) {
    minr <- maxr <- 0
    resPList <- dataList <- fittedList <- list()
    newplot <- TRUE
    for (i in superimpose) {
      if (!is.na(indList[[b]][i])) {
        breakl <- indList[[b]][i]
        group <- grtoplot[[i]][[1]]
        place <- grtoplot[[i]][[2]]
        dset <- group[[place]][2]
        contoplot <- getKinConcen(group, 
                                  multimodel, t, oneDS = place, weight = FALSE)
        spec <- spectraList[[i]]
        resP <- matrix(0, nrow(contoplot), ncol(contoplot))
        for (xx in 1:ncol(contoplot)) resP[, xx] <- contoplot[, 
                                                      xx] * spec[breakl, xx]
        cohspec <- m[[i]]@cohspec
        cohcol <- m[[i]]@cohcol
        if (length(cohcol) >= 1) 
          resP[, cohcol[1]] <- rowSums(as.matrix(resP[, 
                                                      cohcol]))
        if (length(cohcol) >= 2) 
          resP <- resP[, -cohcol[2:length(cohcol)]]
        if (plotoptions@divdrel && length(t[[i]]@drel) != 0) 
          if (length(m[[i]]@dscalspec$perclp) != 0) 
            if (m[[i]]@dscalspec$perclp) {
              resP <- resP/t[[i]]@drel[breakl]
            }
            else {
              resP <- resP/t[[i]]@drel
            }
        data <- m[[i]]@psi.df[, breakl]
        fitted <- res[[i]]@fitted[[breakl]]
        if (divdrel && length(t[[i]]@drel) != 0) 
          if (length(m[[i]]@dscalspec$perclp) != 0) 
            if (m[[i]]@dscalspec$perclp) {
              data <- data/t[[i]]@drel[breakl]
              fitted <- fitted/t[[i]]@drel[breakl]
            }
            else {
              data <- data/t[[i]]@drel
              fitted <- fitted/t[[i]]@drel
            }
        if (m[[i]]@weight) 
          fitted <- fitted/m[[i]]@weightM[, breakl]
        resPList[[length(resPList) + 1]] <- resP
        fittedList[[length(fittedList) + 1]] <- fitted
        dataList[[length(dataList) + 1]] <- data
        for (j in 1:ncol(resPList[[i]])) {
          maxr <- max(na.omit(resP), maxr, fitted, data)
          minr <- min(na.omit(resP), minr, fitted, data)
        }
      }
      else {
        resPList[[length(resPList) + 1]] <- list()
        fittedList[[length(fittedList) + 1]] <- list()
        dataList[[length(dataList) + 1]] <- list()
      }
    }
    if (length(plotoptions@breakdown$ylim) == 0) 
      ylim <- c(minr, maxr)
    else ylim <- plotoptions@breakdown$ylim
    for (i in superimpose) {
      if (!is.na(indList[[b]][i])) {
        breakl <- indList[[b]][i]
        mu <- res[[i]]@irfvec[[breakl]][1]
        data <- resPList[[i]][, j] * spectraList[[i]][breakl, 
                                                      j]
        matlinlogplot(m[[i]]@x, resPList[[i]], mu = mu, 
                      alpha = plotoptions@linrange, add = !newplot, 
                      type = "l", lty = i, ylab = plotoptions@ylab, 
                      xlab = plotoptions@xlab, ylim = ylim, xlim = c(xmin, 
                                                              xmax),
                      main = m[[i]]@x2[breakl])
        newplot <- FALSE
        matlinlogplot(m[[i]]@x, dataList[[i]], mu = mu, 
                      alpha = plotoptions@linrange, add = !newplot, 
                      type = "l", lty = 1, col = i, lwd = 2)
        matlinlogplot(m[[i]]@x, fittedList[[i]], mu = mu, 
                      alpha = plotoptions@linrange, add = !newplot, 
                      type = "l", lty = 2, col = i, lwd = 2)
      }
    }
  }
  if (length(plotoptions@title) != 0) {
    tit <- plotoptions@title
    if (plotoptions@addfilename) 
      tit <- paste(tit, m[[i]]@datafile)
  }
  else {
    tit <- ""
    if (plotoptions@addfilename) 
      tit <- paste(tit, m[[i]]@datafile)
  }
  mtext(tit, side = 3, outer = TRUE, line = 1)
  par(las = 2)
  if (dev.interactive() && length(plotoptions@makeps) != 0) {
    if (plotoptions@output == "pdf") 
      pdev <- pdf
    else pdev <- postscript
    dev.print(device = pdev, file = paste(plotoptions@makeps, 
                               "_breakdown.", plotoptions@output, sep = ""))
  }
}
