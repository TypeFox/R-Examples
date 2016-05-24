
setMethod("ppa.normalize", signature(data="list"),
          function(data, ...) ppa.normalize.default(data, ...))

ppa.normalize.default <- function(data, prenormalize=FALSE) {

  isa.status("PPA normalization", "in")

  if (ncol(data[[1]]) != ncol(data[[2]])) {
    stop("The two matrices must have the same number of columns")
  }
  
  dimnames(data[[1]]) <- dimnames(data[[2]]) <- NULL

  Ecg <- t(scale(t(data[[1]])))
  if (prenormalize) {
    Egc <- t(scale(Ecg))
  } else {
    Egc <- t(scale(data[[1]]))
  }

  Ecd <- t(scale(t(data[[2]])))
  if (prenormalize) {
    Edc <- t(scale(Ecd))
  } else {
    Edc <- t(scale(data[[2]]))
  }

  res <- list(Egc=Egc, Ecd=Ecd, Edc=Edc, Ecg=Ecg)
  attr(res, "hasNA") <- c(any(is.na(Egc)) | any(is.na(Ecd)), 
                          any(is.na(Edc)) | any(is.na(Ecg)))
  attr(res, "prenormalize") <- prenormalize

  isa.status("PPA normalization", "out")
  
  res
}

setMethod("ppa.iterate", signature(normed.data="list"),
          function(normed.data, ...) ppa.iterate.default(normed.data, ...))

## TODO: imeplement the other seeds
ppa.iterate.default <- function(normed.data, row1.seeds, col1.seeds,
                                row2.seeds, col2.seeds,
                                thr.row1, thr.col=thr.row1, thr.row2=thr.row1,
                                direction="updown",
                                convergence=c("cor"), cor.limit=0.9,
                                oscillation=FALSE, maxiter=100) {

  isa.status("Doing PPA iteration", "in")

  if (missing(row1.seeds) && missing(col1.seeds) &&
      missing(row2.seeds) && missing(col2.seeds)) {
    stop("No seeds, nothing to do")
  }
  if (!missing(row1.seeds) && nrow(row1.seeds) != ncol(normed.data$Egc)) {
    stop("Invalid row1 seed length")
  }
  if (!missing(col1.seeds) && nrow(col1.seeds) != ncol(normed.data$Ecd)) {
    stop("Invalid col1 seed length")
  }
  if (!missing(row2.seeds) && nrow(row2.seeds) != ncol(normed.data$Edc)) {
    stop("Invalid row2 seed length")
  }
  if (!missing(col2.seeds) && nrow(col2.seeds) != ncol(normed.data$Ecg)) {
    stop("Invalid col2 seed length")
  }

  if (thr.row1 < 0 || thr.col < 0 || thr.row2 < 0) {
    warning("Negative threshold(s), are you sure about this?")
  }

  direction <- rep(direction, length=4)
  if (any(!direction %in% c("up", "down", "updown"))) {
    stop("Invalid `direction' argument, should be `down', `up' or `updown'.")
  }
  if (direction[1] != direction[3]) {
    warning("`direction[1]' and `direction[3]' refer to the common dimension and they differ")
  }

  convergence <- match.arg(convergence)
  if (convergence <= 0) {
    if (cor.limit <= 0) {
      warning("Non-positive correlation limit for convergence.")
    }
  }

  if (oscillation) {
    stop("Oscillating modules are not yet supported")
  }
  
  no.seeds <- 0
  if (!missing(row1.seeds)) {
    row1.seeds <- cbind(row1.seeds)
    no.seeds <- no.seeds + ncol(row1.seeds)
  }
  if (!missing(col1.seeds)) {
    col1.seeds <- cbind(col1.seeds)
    no.seeds <- no.seeds + ncol(col1.seeds)
  }
  if (!missing(row2.seeds)) {
    row2.seeds <- cbind(row2.seeds)
    no.seeds <- no.seeds + ncol(row2.seeds)
  }
  if (!missing(col2.seeds)) {
    col2.seeds <- cbind(col2.seeds)
    no.seeds <- no.seeds + ncol(col2.seeds)
  }
  
  orig.tg <- thr.row1 ; orig.tc <- thr.col ; orig.td <- thr.row2
  if (length(thr.row1) != 1 && length(thr.row1) != no.seeds) {
    stop("`thr.row1' does not have the right length")
  }
  if (length(thr.col) != 1 && length(thr.col) != no.seeds) {
    stop("`thr.col' does not have the right length")
  }
  if (length(thr.row2) != 1 && length(thr.row2) != no.seeds) {
    stop("`thr.row2' does not have the right length")
  }
  thr.row1 <- rep(thr.row1, length=no.seeds)
  thr.row2 <- rep(thr.row2, length=no.seeds)
  thr.col  <- rep(thr.col , length=no.seeds)

  ## Put all the seeds together
  all.seeds <- matrix(ncol=0, nrow=ncol(normed.data$Egc))
  if (!missing(row1.seeds)) {
    all.seeds <- cbind(all.seeds, row1.seeds)
  }
  if (!missing(col1.seeds)) {
    stop("`col1.seeds' is not implemented yet")
  }
  if (!missing(row2.seeds)) {
    stop("`row2.seeds' is not implemented yet")
  }
  if (!missing(col2.seeds)) {
    stop("`col2.seeds' is not implemented yet")
  }

  ## All the data about this PPA run
  rundata <- list(direction=direction, cor.limit=cor.limit,
                  maxiter=maxiter, N=no.seeds, convergence=convergence,
                  prenormalize=attr(normed.data, "prenormalize"),
                  hasNA=attr(normed.data, "hasNA"), unique=FALSE,
                  oscillation=oscillation)

  ## All the seed data, this will be updated at the end
  seeddata <- data.frame(iterations=NA, oscillation=0,
                         thr.row1=thr.row1, thr.row2=thr.row2,
                         thr.col=thr.col, freq=rep(1, no.seeds),
                         rob=rep(NA, no.seeds))

  if (length(all.seeds)==0) {
    return(list(rows1=matrix(ncol=0, nrow=ncol(normed.data$Egc)),
                rows2=matrix(ncol=0, nrow=nrow(normed.data$Ecd)),
                columns=matrix(ncol=0, nrow=nrow(normed.data$Egc)),
                rundata=rundata, seeddata=seeddata))
  }

  ## Choose convergence checking function
  if (convergence=="cor") {
    check.convergence <- function(row1.old, row1.new, row2.old, row2.new,
                                   col.old, col.new) {
      g.o <- scale(row1.old) ; g.n <- scale(row1.new)
      d.o <- scale(row2.old) ; d.n <- scale(row2.new)
      c.o <- scale(col.old)  ; c.n <- scale(col.new)
      res <- (colSums(g.o * g.n) / (nrow(g.o)-1) > cor.limit &
              colSums(d.o * d.n) / (nrow(d.o)-1) > cor.limit &
              colSums(c.o * c.n) / (nrow(c.o)-1) > cor.limit)
      res <- res & ! is.na(res)
    }
  }

  ## Initialize a couple of things
  iter <- 0
  index <- seq_len(ncol(all.seeds))
  row1.old <- all.seeds
  row2.old <- matrix(NA, nrow=nrow(normed.data$Ecd), ncol=no.seeds)
  col.old  <- matrix(NA, nrow=nrow(normed.data$Egc), ncol=no.seeds)

  row1.res <- matrix(NA, nrow=ncol(normed.data$Egc), ncol=no.seeds)
  row2.res <- matrix(NA, nrow=nrow(normed.data$Ecd), ncol=no.seeds)
  col.res  <- matrix(NA, nrow=nrow(normed.data$Egc), ncol=no.seeds)

  ## Main loop start here
  while (TRUE) {
    iter <- iter + 1
    one.step <- ppa.step(normed.data, row1.old,
                         thr.row1=thr.row1, thr.row2=thr.row2,
                         thr.col=thr.col, direction=direction)

    row1.new <- one.step$rows1
    row2.new <- one.step$rows2
    col.new  <- one.step$columns

    ## Mark converged seeds
    conv <- check.convergence(row1.old, row1.new, row2.old, row2.new,
                              col.old, col.new)

    ## Mar kall zero seeds
    zero <- apply(row1.new, 2, function(x) all(x==0))

    ## These are thrown out
    drop <- which(conv | zero)

    ## Drop the seeds to be dropped
    if (length(drop) != 0) {
      row1.res[,index[drop]] <- row1.new[,drop]
      row2.res[,index[drop]] <- row2.new[,drop]
      col.res[,index[drop]]  <- col.new[,drop]
      seeddata$iterations[index[drop]] <- iter
      row1.new <- row1.new[,-drop,drop=FALSE]
      row2.new <- row2.new[,-drop,drop=FALSE]
      col.new  <- col.new[,-drop,drop=FALSE]
      thr.row1 <- thr.row1[-drop]
      thr.row2 <- thr.row2[-drop]
      thr.col  <- thr.col [-drop]
      index <- index[-drop]
    }

    if (ncol(row1.new)==0 || iter > maxiter) { break; }

    row1.old <- row1.new
    row2.old <- row2.new
    col.old <- col.new

  } ## End of main looop

  isa.status("isa iteration", "out")

  list(rows1=row1.res, rows2=row2.res, columns=col.res,
       rundata=rundata, seeddata=seeddata)
}

ppa.step <- function(normed.data, rows1, thr.row1, thr.row2, thr.col,
                     direction) {

  direction <- rep(direction, length=4)
  if (any(!direction %in% c("up", "updown", "down"))) {
    stop("Invalid `direction' argument")
  }

  filter <- function(x, t, dir) {
    if (dir=="updown") {
      x <- .Call("beta_filter_updown_vart", x, as.double(t), PACKAGE="isa2")
    } else if (dir=="up") {
      x <- .Call("beta_filter_up_vart", x, as.double(t), PACKAGE="isa2")
    } else { ## dir=="down"
      x <- .Call("beta_filter_down_vart", x, as.double(t), PACKAGE="isa2")
    }
  }

  hasNA <- c(TRUE, TRUE)
  if ("hasNA" %in% names(attributes(normed.data))) {
    hasNA <- attr(normed.data, "hasNA")
  }

  ## Egc * g
  if (!hasNA[1]) {
    col.new <- filter(normed.data$Egc %*% rows1, thr.col, direction[1])
  } else {
    col.new <- filter(na.multiply(normed.data$Egc, rows1), thr.col,
                      direction[1])
  }
  
  ## Ecd * c
  if (!hasNA[2]) {
    row2.new <- filter(normed.data$Ecd %*% col.new, thr.row2, direction[2])
  } else {
    row2.new <- filter(na.multiply(normed.data$Ecd, col.new), thr.row2,
                       direction[2])
  }

  ## Edc * d
  if (!hasNA[2]) {
    col.new <- filter(normed.data$Edc %*% row2.new, thr.col, direction[3])
  } else {
    col.new <- filter(na.multiply(normed.data$Edc, row2.new), thr.col,
                      direction[3])
  }

  ## Ecg * c
  if (!hasNA[1]) {
    row1.new <- filter(normed.data$Ecg %*% col.new, thr.row1, direction[4])
  } else {
    row1.new <- filter(na.multiply(normed.data$Ecg, col.new), thr.row1,
                       direction[4])
  }

  list(columns=col.new, rows1=row1.new, rows2=row2.new)
}

setMethod("ppa.unique", signature(normed.data="list", pparesult="list"),
          function(normed.data, pparesult, ...)
          ppa.unique.default(normed.data, pparesult, ...))

ppa.unique.default <- function(normed.data, pparesult, method=c("cor"),
                               ignore.div=TRUE, cor.limit=0.9, neg.cor=TRUE,
                               drop.zero=TRUE) {

  method <- match.arg(method)
  
  if (ncol(pparesult$columns) == 0) { return(pparesult) }

  pparesult$rundata$unique <- TRUE
  
  drop <- rep(FALSE, ncol(pparesult$columns))
              
  if (ignore.div) { drop <- drop | is.na(pparesult$seeddata$iterations) }
  if (drop.zero)  { drop <- drop | (colSums(pparesult$rows1!=0) == 0)   }

  if (any(drop)) {
    pparesult$rows1 <- pparesult$rows1[,!drop,drop=FALSE]
    pparesult$rows2 <- pparesult$rows2[,!drop,drop=FALSE]
    pparesult$columns <- pparesult$columns[,!drop,drop=FALSE]
    pparesult$seeddata <- pparesult$seeddata[!drop,,drop=FALSE]
  }

  if (method=="cor") {
    if (cor.limit < 1 && ncol(pparesult$rows1) > 1) {
      cm <- cor(pparesult$rows1)
      if (neg.cor) { cm <- abs(cm) }
      cm[lower.tri(cm, diag=TRUE)] <- 0
      uni <- apply(cm < cor.limit, 2, all)
      freq <- sapply(seq_len(nrow(cm)), function(x) {
        sum(pparesult$seeddata$freq[cm[x,] >= cor.limit])
      }) + pparesult$seeddata$freq
      
      pparesult$rows1 <- pparesult$rows1[,uni,drop=FALSE]
      pparesult$rows2 <- pparesult$rows2[,uni,drop=FALSE]
      pparesult$columns <- pparesult$columns[,uni,drop=FALSE]
      pparesult$seeddata <- pparesult$seeddata[uni,,drop=FALSE]
      pparesult$seeddata$freq <- freq[uni]
    }
  }

  pparesult
}

setMethod("ppa", signature(data="list"),
          function(data, ...) ppa.default(data, ...))

ppa.default <- function(data,
                        thr.row1=seq(1,3,by=0.5),
                        thr.row2=seq(1,3,by=0.5),
                        thr.col=seq(1,3,by=0.5),
                        no.seeds=100, direction="updown") {

  isa.status("Performing complete PPA work flow", "in")

  if (!is.list(data) || length(data) != 2 ||
      !is.matrix(data[[1]]) || !is.matrix(data[[2]]) ||
      ncol(data[[1]]) != ncol(data[[2]])) {
    stop("`data' must be a list of two matrices, with matching number of columns")
  }

  ## Normalize the input
  normed.data <- ppa.normalize(data)

  ## Generate seeds
  row1.seeds <- generate.seeds(length=nrow(data[[1]]), count=no.seeds)

  ## Determine thresholds
  thr.list <- expand.grid(thr.row1=thr.row1, thr.row2=thr.row2,
                          thr.col=thr.col)
  thr.list <- unlist(apply(thr.list, 1, list), recursive=FALSE)

  ## Do the PPA, for all thresholds
  pparesults <- lapply(thr.list, function(x)
                       ppa.iterate(normed.data, row1.seeds=row1.seeds,
                                   thr.row1=x["thr.row1"],
                                   thr.row2=x["thr.row2"],
                                   thr.col=x["thr.col"],
                                   direction=direction))

  ## Make it unique for every threshold combination
  pparesults <- lapply(pparesults, function(x) ppa.unique(normed.data, x))

  ## Filter according to robustness
  pparesults <- lapply(pparesults, function(x)
                       ppa.filter.robust(data=data,
                                         normed.data=normed.data,
                                         ppares=x,
                                         row1.seeds=row1.seeds))
  
  ## Merge them
  result <- list()
  result$rows1 <- do.call(cbind, lapply(pparesults, "[[", "rows1"))
  result$rows2 <- do.call(cbind, lapply(pparesults, "[[", "rows2"))
  result$columns <- do.call(cbind, lapply(pparesults, "[[", "columns"))
  result$seeddata <- do.call(rbind, lapply(pparesults, "[[", "seeddata"))
  result$rundata <- pparesults[[1]]$rundata
  result$rundata$N <- sum(sapply(pparesults, function(x) x$rundata$N))

  ## Another filtering
  result <- ppa.unique(normed.data, result)
  
  isa.status("DONE", "out")

  result
}
