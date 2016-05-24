# next release:
# make index a list
# make nbreaks a vector
# public functions: samur, summary.samur, plot.samur
# private functions: mdh, predict.mdh, generate.sample.wrapper, generate.sample.core

samur <- function(formula, data, matched.subset, nsmp = 100, use.quantile = TRUE, breaks = 10) {
  mycall <- match.call()
  # data checks:
  # 1) formula (factor response)
  
  # 2) matched.subset
  ndat <- nrow(data)
  matched.subset <- as.integer(matched.subset)
  if (length(matched.subset) > ndat) stop("matched subset cannot be larger than number of rows in data")
  if (any(matched.subset < 1 | matched.subset > ndat)) stop("out of range values in matched subset")
  if (length(unique(matched.subset)) < length(matched.subset)) stop("cannot have duplicate values in matched subset")
  
  # 3) nsmp: positive integer
  nsmp <- as.integer(nsmp)
  if (nsmp < 1) stop("nsmp must be a positive integer")
  
  # 4) breaks: positive integer larger than 1 and smaller than ??
  breaks <- as.integer(breaks)
  if (breaks < 1) stop("breaks must be a positive integer")
  
  my.mdh <- mdh(formula, data, breaks = breaks, use.quantile = use.quantile)
  myhist <- table(predict(my.mdh, data[matched.subset, ])$cell.assignment)
  mysmps <- generate.sample.wrapper(my.mdh, myhist, data, nsmp = nsmp)
  
  attr(mysmps, "call") <- mycall
  attr(mysmps, "formula") <- formula
  attr(mysmps, "mdg") <- my.mdh
  attr(mysmps, "mdh") <- myhist
  attr(mysmps, "data") <- data
  attr(mysmps, "matched.subset") <- matched.subset
  class(mysmps) <- c("samur", class(mysmps))
  
  return (mysmps)
}

print.summary.samur <- function(x, ...) {
  cat("minimum p-value of original matched subset:", x$min.pval.orig, "\n")
  cat("range of minimum p-values for augmented set:", range(x$min.pval.new), "\n")
  cat("coverage of original set:", x$coverage.orig, "\n")
  cat("coverage of augmented set:", x$coverage.new, "\n")
}

summary.samur <- function(object, ...) {
  min.pval.new <- sapply(1:ncol(object), function(n) {
    MatchBalance(attr(object, "formula"), attr(object, "data")[object[, n], ], print.level = 0, ...)$BMsmallest.p.value
    })
  min.pval.orig <- MatchBalance(attr(object, "formula"), attr(object, "data")[attr(object, "matched.subset"), ], print.level = 0, ...)$BMsmallest.p.value
  coverage.orig <- length(attr(object, "matched.subset")) / nrow(attr(object, "data"))
  coverage.new <- length(unique(as.vector(object))) / nrow(attr(object, "data"))
  ret <- list(min.pval.new = min.pval.new, min.pval.orig = min.pval.orig, coverage.orig = coverage.orig, coverage.new = coverage.new)
  class(ret) <- c("summary.samur", class(ret))
  return (ret)
}

print.samur <- function(x, ...) {
  cat("Call:\n")
  print(attr(x, "call"))
}

generate.sample.wrapper <- function(object, mytgt, newdata, nsmp = 100) { # remove formula
  sapply(1:nsmp, function(n) generate.sample.core(object, mytgt, newdata))
}

generate.sample.core <- function(object, mytgt, newdata) {
  nTreat <- nrow(mytgt)
  treatLevels <- rownames(mytgt)
  nCells <- ncol(mytgt)
  pred <- predict(object, newdata)
  idx <- pred$cell.assignment$cellno
  treatCol <- object$treatCol

  ret <- unlist(sapply(1:nTreat, function(m) {
    treat <- treatLevels[m]
    unlist(sapply(1:nCells, function(i) {
      #sample(which((idx==i) & (data[,treatCol]==treat)), size=tgtList[[m]][[i]][n], replace=F)
      if (mytgt[m,i] > 0) {
        candidate.set <- which((idx==i) & (newdata[,treatCol]==treat))
        if (length(candidate.set) < mytgt[m,i]) stop("requesting more samples than available data")
        if (length(candidate.set) == mytgt[m,i]) return (candidate.set) # combination of this and above line should take care of nasty bug where length(candidate.set)==1
        return (sample(candidate.set, size=mytgt[m,i], replace=F))
      } else {
        return (c())
      }
    }))
  }))
  return (as.vector(ret))
}

predict.mdh <- function(object, newdata, ...) {
  # we need more checks, e.g. that no new treatments appear in newdata, that treatment column type is compatible with training, etc
  lookup.obs <- function(x, part, i) {
    if (is.factor(x)) {
      return (x %in% part[[i]])
    } else {
      if (i==1) return (x >= part[i] & x <= part[i+1])
      return (x > part[i] & x <= part[i+1])
    }
  }
  if (missing(newdata)) {
    mf <- object$modelframe
  } else {
    mf <- model.frame(object$modelterms, droplevels(newdata), xlev = object$xlevels)
  }
  nobs <- nrow(mf)
  ncells <- prod(object$nparts)
  nmatch <- length(object$matchCols)
  
  idx2 <- rep(NA, nobs) # determining assignment of each observation to a cell
  tmpsink <- sapply(1:ncells, function(i) {
    sel <- rep(T,nobs)
    for (j in 1:nmatch) sel <- sel & lookup.obs(mf[,object$matchCols[j]], object$parts[[j]], object$cell.to.part.map[i,j])
    idx2[sel] <<- i # intentional "<<-", don't change to "<-"!!!
  })
  if (any(is.na(idx2))) warning("one or more cases were not assigned to histogram cells")
  treatCol <- object$treatCol
  cell.assignment <- data.frame(mf[,treatCol], idx2)
  colnames(cell.assignment) <- c(treatCol, "cellno") # rename to case assignment
  cell.assignment$cellno <- factor(cell.assignment$cellno, levels = 1:ncells)
  
  ret <- list(isFactor = object$isFactor, parts = object$parts, nparts = object$nparts
              , cell.to.part = object$cell.to.part.map, cell.assignment = cell.assignment)
  class(ret) <- c("predict.mdh", class(ret))
  return (ret)
}

# function for constructing multi-dimensional histograms
# TODO: determine if two-part formulas are needed and handle them if needed
mdh <- function(formula, data, breaks = 5, use.quantile = TRUE) {
  my.mf <- model.frame(formula, data, drop.unused.levels=TRUE, na.action = na.fail)
  my.mt <- attr(my.mf, "terms")
  my.levels <- .getXlevels(my.mt, my.mf)
  treatCol <- colnames(my.mf)[attr(my.mt, "response")]
  matchCols <- attr(my.mt, "term.labels")
  
  isFactor <- sapply(matchCols, function(x) is.factor(my.mf[,x]))
  
  nMatch <- length(matchCols)
  partition <- list()
  nparts <- rep(NA, nMatch)
  cellGrid.gen <- list()
  breaks.vec <- rep(NA, nMatch)
  for (i in 1:nMatch) {
    if (isFactor[i]) { # partitioning a factor variable
      partition[[i]] <- list()
      partition[[i]] <- my.levels[[matchCols[i]]]
      nparts[i] <- length(partition[[i]])
    } else { # partitioning a numeric variable
      breaks <- min(breaks, length(unique(my.mf[,matchCols[i]])))
      breaks.vec[i] <- breaks
      if (use.quantile) {
        partition[[i]] <- quantile(my.mf[,matchCols[i]], probs = seq(from=0.0, to=1.0, length.out=breaks+1))
      } else {
        partition[[i]] <- hist(my.mf[,matchCols[i]], plot = F, breaks = breaks)$breaks
      }
      nparts[i] <- length(partition[[i]])-1
    }
    cellGrid.gen[[i]] <- 1:nparts[i]
  }
  
  cellGrid <- expand.grid(as.list(cellGrid.gen)) # is "as.list" necessary?
  colnames(cellGrid) <- matchCols
  
  ret <- list(treatCol = treatCol, matchCols = matchCols, modelframe = my.mf, modelterms = my.mt, xlevels = my.levels
    , cell.to.part.map = cellGrid, parts = partition, nparts = nparts, ncells = prod(nparts), isFactor = isFactor
    , breaks.vec = breaks.vec)
  
  class(ret) <- "mdh"
  return (ret)
}

