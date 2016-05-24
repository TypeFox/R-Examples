
setMethod("isa.normalize", signature(data="matrix"),
          function(data, ...) isa.normalize.default(data, ...))

isa.normalize.default <- function(data, prenormalize=FALSE) {

  isa.status("ISA normalization", "in")

  if (!is.matrix(data)) {
    stop("`data' must be a matrix")
  }
  
  ## Then normalize it
  Ec <- scale(t(data))
  if (prenormalize) {
    Er <- scale(t(Ec))
  } else {
    Er <- scale(data)
  }

  data <- list(Er=t(Er), Ec=t(Ec))
  
  attr(data, "prenormalize") <- prenormalize
  attr(data, "hasNA") <- (any(is.na(Er)) |
                          any(is.na(Ec)) )

  isa.status("DONE", "out")
  
  data
}

setMethod("isa.iterate", signature(normed.data="list"),
          function(normed.data, ...) isa.iterate.default(normed.data, ...))

isa.iterate.default <- function(normed.data, row.seeds, col.seeds,
                                thr.row, thr.col=thr.row,
                                direction=c("updown", "updown"),
                                convergence=c("corx", "cor", "eps"),
                                cor.limit=0.99, eps=1e-4,
                                corx=3,
                                oscillation=FALSE, maxiter=100) {

  isa.status("Doing ISA iteration", "in")

  if (( missing(row.seeds) &&  missing(col.seeds))) {
    stop("No seeds, nothing to do")
  }
  if (!missing(row.seeds) && nrow(row.seeds) != ncol(normed.data$Er)) {
    stop("Invalid row seed length")
  }
  if (!missing(col.seeds) && nrow(col.seeds) != ncol(normed.data$Ec)) {
    stop("Invalid column seed length")
  }
  
  if (thr.row < 0 || thr.col < 0) {
    warning("Negative thresholds, are you sure about this?")
  }
  
  direction <- rep(direction, length=2)
  if (any(!direction %in% c("up", "down", "updown"))) {
    stop("Invalid `direction' argument, should be `down', `up' or `updown'.")
  }

  convergence <- match.arg(convergence)
  if (convergence == "cor") {
    if (cor.limit <= 0) {
      warning("Non-positive correlation limit for convergence.")
    }
  }
  if (convergence == "eps") {
    if (eps >= 1) {
      warning("`eps' limit for convergence greater than one.")
    }
  }
  
  no.seeds <- 0
  if (!missing(row.seeds)) {
    no.seeds <- no.seeds + ncol(row.seeds)
  }
  if (!missing(col.seeds)) {
    no.seeds <- no.seeds + ncol(col.seeds)
  }
  
  orig.tg <- thr.row
  orig.tc <- thr.col
  if (length(thr.row) != 1 && length(thr.row) != no.seeds) {
    stop("`thr.row' does not have the right length")
  }
  if (length(thr.col) != 1 && length(thr.col) != no.seeds) {
    stop("`thr.col' does not have the right length")
  }
  thr.row <- rep(thr.row, length=no.seeds)
  thr.col <- rep(thr.col, length=no.seeds)

  ## Put the seeds together
  all.seeds <- matrix(ncol=0, nrow=nrow(normed.data$Ec))
  if (!missing(row.seeds)) {
    all.seeds <- cbind(all.seeds, row.seeds)
  }
  if (!missing(col.seeds)) {
    col.seeds <- isa.row.from.col(normed.data, col.seeds=col.seeds,
                                  thr.row=tail(thr.row, ncol(col.seeds)),
                                  direction=direction[2])
    all.seeds <- cbind(all.seeds, col.seeds)
  }

  ## All the data about this ISA run
  rundata <- list(direction=direction, eps=eps, cor.limit=cor.limit,
                  maxiter=maxiter, N=no.seeds, convergence=convergence,
                  prenormalize=attr(normed.data, "prenormalize"),
                  hasNA=attr(normed.data, "hasNA"), corx=corx,
                  unique=FALSE, oscillation=oscillation)

  ## All the seed data, this will be updated, of course
  seeddata <- data.frame(iterations=NA, oscillation=0,
                         thr.row=thr.row, thr.col=thr.col,
                         freq=rep(1, no.seeds), rob=rep(NA, no.seeds))
  
  if (length(all.seeds)==0) {
    return(list(rows=all.seeds, columns=matrix(ncol=0, nrow=ncol(normed.data$Ec)),
                rundata=rundata, seeddata=seeddata))
  }

  ## Choose convergence checking function
  if (convergence=="eps") {
    check.convergence <- function(row.old, row.new, col.old, col.new) {
      res <- (apply(row.old-row.new, 2, function(x) all(abs(x)<eps)) &
              apply(col.old-col.new, 2, function(x) all(abs(x)<eps)))
      res & !is.na(res)
    }
  } else if (convergence=="cor") {
    check.convergence <- function(row.old, row.new, col.old, col.new) {
      g.o <- scale(row.old)
      g.n <- scale(row.new)
      c.o <- scale(col.old)
      c.n <- scale(col.new)
      res <- (colSums(g.o * g.n) / (nrow(g.o)-1) > cor.limit &
              colSums(c.o * c.n) / (nrow(c.o)-1) > cor.limit)
      res & !is.na(res)
    }
  } else if (convergence=="corx") {
    if (corx < 2) { stop("Invalid `corx' value, shoudl be at least 2") }
    rows.old <- cols.old <- list()
    idx.old <- seq_len(corx)
    check.convergence <- function(row.old, row.new, col.old, col.new) {
      if (iter < corx+1) {
        rows.old <<- c(rows.old, list(row.new))
        cols.old <<- c(cols.old, list(col.new))
        return(rep(FALSE, ncol(row.old)))
      }
      row.new <- scale(row.new)
      col.new <- scale(col.new)
      res <- (colSums(rows.old[[idx.old[1]]]*row.new)/(nrow(row.new)-1) > cor.limit &
              colSums(cols.old[[idx.old[1]]]*col.new)/(nrow(col.new)-1) > cor.limit)
      rows.old[[ idx.old[1] ]] <<- row.new
      cols.old[[ idx.old[1] ]] <<- col.new
      idx.old <<- c(idx.old[-1], idx.old[1])
      res & !is.na(res)
    }
  }

  ## Initialize a couple of things
  iter <- 0
  index <- seq_len(ncol(all.seeds))
  if (oscillation) { fire <- character(no.seeds) }
  row.old <- all.seeds
  col.old <- matrix(NA, nrow=ncol(normed.data$Ec), ncol=no.seeds)
  row.res <- matrix(NA, nrow=nrow(normed.data$Ec), ncol=no.seeds)
  col.res <- matrix(NA, nrow=ncol(normed.data$Ec), ncol=no.seeds)
  
  ## Main loop starts here
  while (TRUE) {

    iter <- iter + 1
    one.step <- isa.step(normed.data, rows=row.old, thr.row=thr.row,
                         thr.col=thr.col, direction=direction)

    row.new <- one.step$rows
    col.new <- one.step$columns

    ## Mark converged seeds
    conv <- check.convergence(row.old=row.old, row.new=row.new,
                              col.old=col.old, col.new=col.new)

    ## Mark all zero seeds
    zero <- apply(row.new, 2, function(x) all(x==0))
    
    ## Mark oscillating ones, if requested
    if (oscillation && iter > 1) {
      new.fire <- apply(row.new, 2, function(x) sum(round(x, 4)))
      fire <- paste(sep=":", fire, new.fire)
      osc <- logical(ncol(row.new))
      osc[ (mat <- regexpr("(:.*:.*)\\1$", fire)) != -1] <- TRUE
      osc <- osc & !conv
      
      if (any(osc)) {
        mat <- cbind(mat[osc], attr(mat, "match.length")[[1]][osc])
        mat <- sapply(seq(length=nrow(mat)), function(x) substr(fire[osc][x],
                            mat[x,1], mat[x,1]+mat[x,2]))
        mat <- sapply(mat, function(x) sum(utf8ToInt(x) == 58), USE.NAMES=FALSE )
        seeddata$oscillation[index[osc]] <- mat/2
      }
    } else {
      osc <- FALSE
    }
    
    ## These are all to throw out
    drop <- which(conv | zero | osc)
    
    ## Drop the seeds to be dropped
    if (length(drop) != 0) {
      row.res[,index[drop]] <- row.new[,drop]
      col.res[,index[drop]] <- col.new[,drop]
      seeddata$iterations[index[drop]] <- iter
      row.new <- row.new[,-drop,drop=FALSE]
      col.new <- col.new[,-drop,drop=FALSE]
      if (oscillation) { fire <- fire[-drop] }
      thr.row <- thr.row[-drop]
      thr.col <- thr.col[-drop]
      if (convergence=="corx") {
        rows.old <- lapply(rows.old, function(x) x[,-drop,drop=FALSE])
        cols.old <- lapply(cols.old, function(x) x[,-drop,drop=FALSE])
      }
      index <- index[-drop]
    }
    
    if (ncol(row.new)==0 || iter>maxiter) { break; }

    row.old <- row.new
    col.old <- col.new
    
  } ## End of main loop
  
  isa.status("DONE", "out")
  
  list(rows=row.res, columns=col.res,
       rundata=rundata, seeddata=seeddata)
}

na.multiply <- function(A, B) {
  M <- !is.na(A)
  modA <- A
  modA[!M] <- 0
  v <- modA %*% B
  w <- sqrt(M %*% B^2)
  w2 <- sqrt(apply(B^2, 2, sum))
  ret <- v/w * rep(w2, each=nrow(v))
  ret[ w==0 ] <- 0
  ret
}

isa.step <- function(normed.data, rows, thr.row, thr.col, direction) {

  Ec <- normed.data$Ec
  Er <- normed.data$Er

  direction <- rep(direction, length=2)
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

  if ("hasNA" %in% names(attributes(normed.data)) && !attr(normed.data, "hasNA")) {
    col.new <- filter(Er %*% rows,    thr.col, direction[1])
    row.new <- filter(Ec %*% col.new, thr.row, direction[2])
  } else {
    col.new <- filter(na.multiply(Er, rows   ), thr.col, direction[1])
    row.new <- filter(na.multiply(Ec, col.new), thr.row, direction[2])
  }

  list(columns=col.new, rows=row.new)
}

setMethod("isa.unique", signature(normed.data="list", isaresult="list"),
          function(normed.data, isaresult, ...)
          isa.unique.default(normed.data, isaresult, ...))

isa.unique.default <- function(normed.data, isaresult, method=c("cor"),
                               ignore.div=TRUE, cor.limit=0.9, neg.cor=TRUE,
                               drop.zero=TRUE) {

  isa.status("Creating unique module set", "in")
  
  method <- match.arg(method)

  if (ncol(isaresult$rows) == 0) { return(isaresult) }

  ## drop divergent seeds
  if (ignore.div) {
    invalid <- is.na(isaresult$seeddata$iterations)
    if (any(invalid)) {
      valid <- !invalid
      isaresult$rows <- isaresult$rows[,valid,drop=FALSE]
      isaresult$columns <- isaresult$columns[,valid,drop=FALSE]
      isaresult$seeddata <- isaresult$seeddata[valid,,drop=FALSE]
    }
  }
  if (ncol(isaresult$rows) == 0) { return(isaresult) }
  
  ## drop all zero seeds
  if (drop.zero) {
    valid <- apply(isaresult$rows, 2, function(x) any(x != 0))
    if (!all(valid)) {
      isaresult$rows <- isaresult$rows[,valid,drop=FALSE]
      isaresult$columns <- isaresult$columns[,valid,drop=FALSE]
      isaresult$seeddata <- isaresult$seeddata[valid,,drop=FALSE]
    }
  }
  if (ncol(isaresult$rows) == 0) { return(isaresult) }

  if (method=="cor") {
    ## We reorder the results a bit, because we want to keep modules
    ## with higher (row,column) thresholds
    ord <- order(isaresult$seeddata$thr.row, isaresult$seeddata$thr.col,
                 decreasing=TRUE)
    isaresult$rows <- isaresult$rows[,ord,drop=FALSE]
    isaresult$columns <- isaresult$columns[,ord,drop=FALSE]
    isaresult$seeddata <- isaresult$seeddata[ord,,drop=FALSE]
    
    if (neg.cor) { ABS <- abs } else { ABS <- function(x) x }
    cm <- pmin(ABS(cor(isaresult$rows)), ABS(cor(isaresult$columns)))
    cm[ lower.tri(cm, diag=TRUE) ] <- 0
    uni <- apply(cm < cor.limit, 2, all)
    freq <- sapply(seq_len(nrow(cm)), function(x) {
      sum(isaresult$seeddata$freq[cm[x,] >= cor.limit])
    }) + isaresult$seeddata$freq
  } else if (method=="round") {
    ## TODO
    stop("The `round' method is currently not implemented")
  }

  isaresult$rows <- isaresult$rows[,uni,drop=FALSE]
  isaresult$columns <- isaresult$columns[,uni,drop=FALSE]

  isaresult$seeddata <- isaresult$seeddata[uni,,drop=FALSE]
  isaresult$seeddata$freq <- freq[uni]

  isaresult$rundata$unique <- TRUE

  isa.status("DONE", "out")
  
  isaresult      
}


isa.row.from.col <- function(normed.data, col.seeds, thr.row, direction) {


  Ec <- normed.data$Ec

  if (! direction %in% c("up", "updown", "down")) {
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

  if ("hasNA" %in% names(attributes(normed.data)) && !attr(normed.data, "hasNA")) {
    row.new <- filter(Ec %*% col.seeds, thr.row, direction)
  } else {
    row.new <- filter(na.multiply(Ec, col.seeds), thr.row, direction)
  }

  row.new
}  

generate.seeds <- function(length, count=100, method=c("uni"),
                           sparsity=2) {

  isa.status("Generating seeds", "in")

  if (method == "uni") {
    sparsity <- rep(sparsity, length=count)
    g <- matrix(0, nrow=length, ncol=count)
    for (i in 1:count) {
      g[sample(length, sparsity[i]),i] <- 1
    }
  }

  isa.status("DONE", "out")
  
  g
}

setMethod("isa.sweep", signature(data="matrix"),
          function(data, ...) isa.sweep.default(data, ...))

isa.sweep.default <- function(data, normed.data, isaresult, method=c("cor"),
                              neg.cor=TRUE, cor.limit=0.9) {
  
  isa.status("Performing an ISA sweep", "in")
  
  method <- match.arg(method)

  if (any(isaresult$seeddata$thr.row != isaresult$seeddata$thr.row[1]) &&
      any(isaresult$seeddata$thr.col != isaresult$seeddata$thr.col[1])) {
    stop(paste("Either the row or the column threshold",
               "must be constant in `isaresult'"))
  }

  if (all(isaresult$seeddata$thr.row==isaresult$seeddata$thr.row[1])) {
    key <- isaresult$seeddata$thr.col
  } else {
    key <- isaresult$seeddata$thr.row
  }

  isalist <- tapply(seq_along(key), key, list)
  isalist.value <- as.numeric(names(isalist))
  ord <- order(isalist.value, decreasing=TRUE)
  isalist <- isalist[ord]
  isalist.value <- isalist.value[ord]

  isalist.row <- lapply(isalist, function(x) isaresult$rows[,x,drop=FALSE])
  isalist.col <- lapply(isalist, function(x) isaresult$columns[,x,drop=FALSE])
  isalist.seed <- lapply(isalist, function(x) isaresult$seeddata[x,,drop=FALSE])

  thr.row <- sapply(isalist, function(x) isaresult$seeddata$thr.row[x[1]])
  thr.col <- sapply(isalist, function(x) isaresult$seeddata$thr.col[x[1]])

  NN <- isaresult$rundata$N
  
  conv.to <- list(numeric())
  for (i in seq_along(isalist)[-1]) {

    ## Are there any seeds at all?
    if (ncol(isalist.row[[i-1]]) == 0) {
      conv.to[[i]] <- numeric()
      next
    }

    ## Run ISA on the modules from step i-1 with threshold from step i
    NN <- NN + ncol(isalist.row[[i-1]])
    tmpres <- isa.iterate(normed.data, row.seeds=isalist.row[[i-1]],
                          thr.row=thr.row[i], thr.col=thr.col[i],
                          direction=isaresult$rundata$direction,
                          convergence=isaresult$rundata$convergence,
                          eps=isaresult$rundata$eps,
                          cor.limit=isaresult$rundata$cor.limit,
                          corx=isaresult$rundata$corx,
                          maxiter=isaresult$rundata$maxiter,
                          oscillation=isaresult$rundata$oscillation)

    tmpures <- isa.unique(normed.data, tmpres, method=method,
                         cor.limit=cor.limit, drop.zero=TRUE,
                         ignore.div=TRUE)

    ## instead of doing permutations, we just calculate the robustness
    ## and use the included robustsness limit
    if (!is.null(isaresult$seeddata$rob.limit)) {
      rob <- robustness(normed.data, tmpures$rows, tmpures$columns)
      keep <- rob > max(isalist.seed[[i]]$rob.limit)
      tmpures$rows <- tmpures$rows[,keep,drop=FALSE]
      tmpures$columns <- tmpures$columns[,keep,drop=FALSE]
      tmpures$seeddata <- tmpures$seeddata[keep,,drop=FALSE]
      tmpures$seeddata$rob <- rob[keep]
      tmpures$seeddata$rob.limit <- NA
      tmpures$seeddata$rob.limit[] <- max(isalist.seed[[i]]$rob.limit)
    }

    if (neg.cor) { ABS <- abs } else { ABS <- function(x) x }

    if (method=="cor") {
      ## Add the newly found modules
      if (ncol(tmpures$rows) != 0) {
        if (ncol(isalist.row[[i]]) != 0) {
          cm <- pmin(ABS(cor(isalist.row[[i]], tmpures$rows)),
                     ABS(cor(isalist.col[[i]], tmpures$columns)))
          first <- apply(cm < cor.limit, 2, all)
          isalist.row[[i]] <- cbind(isalist.row[[i]],
                                    tmpures$rows[,first,drop=FALSE])
          isalist.col[[i]] <- cbind(isalist.col[[i]],
                                    tmpures$columns[,first,drop=FALSE])
          isalist.seed[[i]] <- rbind(isalist.seed[[i]],
                                     tmpures$seeddata[first,,drop=FALSE])
        } else {
          isalist.row[[i]] <- tmpures$rows
          isalist.col[[i]] <- tmpures$columns
          isalist.seed[[i]] <- tmpures$seeddata
        }
      }
      
      ## Check what converged to what
      if (ncol(isalist.row[[i]]) != 0 && ncol(tmpres$rows) != 0) {
        cm <- pmin(ABS(cor(isalist.row[[i]], tmpres$rows)),
                   ABS(cor(isalist.col[[i]], tmpres$columns)))
        conv.to[[i]] <- apply(cm, 2, function(x) which(x >= cor.limit)[1])          
      } else {
        conv.to[[i]] <- numeric()
      }
    }
  } ## for

  v <- 0; for (i in seq_along(conv.to)) {
    v <- v + length(conv.to[[i]])
    conv.to[[i]] <- conv.to[[i]] + v
  }
  
  result <- list()
  result$rows <- do.call(cbind, isalist.row)
  result$columns <- do.call(cbind, isalist.col)
  result$seeddata <- do.call(rbind, isalist.seed)
  result$seeddata$father <- c(unlist(conv.to),
                              rep(NA, ncol(isalist.row[[i]])))
  result$seeddata$level <- rep(seq_along(isalist.row),
                               sapply(isalist.row, ncol))

  result$rundata <- isaresult$rundata
  result$rundata$N <- NN

  isa.status("DONE", "out")
  
  result
}

setMethod("sweep.graph", signature(sweep.result="list"),
          function(sweep.result, ...) sweep.graph.default(sweep.result, ...))

sweep.graph.default <- function(sweep.result) {

  if (is.null(sweep.result$seeddata$father)) {
    stop("Not a sweep result")
  }

  if (!requireNamespace("igraph")) {
    stop("The igraph package is required for sweeping")
  }

  nnodes <- nrow(sweep.result$seeddata)

  from <- seq_len(nnodes)
  to <- sweep.result$seeddata$father

  valid <- !is.na(to)
  from <- from[valid]
  to <- to[valid]

  G <- igraph::graph( rbind(from, to), n=nnodes )

  if (length(unique(sweep.result$seeddata$thr.row)) == 1) {
    igraph::V(G)$thr <- sweep.result$seeddata$thr.col
  } else {
    igraph::V(G)$thr <- sweep.result$seeddata$thr.row
  }

  igraph::V(G)$id <- seq(igraph::vcount(G))
  graphs <- igraph::decompose.graph(G)

  layouts <- lapply(graphs, function(g) {
    l <- igraph::layout.reingold.tilford(
      g, root=tail(igraph::topological.sort(g, mode="out"), 1))
    r <- sqrt(l[,1]^2 + l[,2]^2)
    phi <- atan2(l[,2], l[,1]) - pi/2
    l[,1] <- r * cos(phi)
    l[,2] <- r * sin(phi)
    labels <- sort(unlist(igraph::V(g)$thr))
    l <- igraph::layout.norm(l, labels[1]-2, labels[length(labels)]-2, NULL, NULL)
    l[,2] <- l[,2] - min(l[,2])
    l})

  offs <- 0
  for (i in 1:length(layouts)) {
    layouts[[i]][,2] <- layouts[[i]][,2] + offs
    r <- range(layouts[[i]][,2])
    r[2] <- r[2] + 1
    offs <- offs + r[2] - r[1]
  }
  offs <- offs-0.5

  G <- igraph::set.graph.attribute(G, "layout", do.call(rbind, layouts))
  lay <- igraph::get.graph.attribute(G, "layout")
  lay[unlist(sapply(graphs, igraph::get.vertex.attribute, "id")),] <- lay
  G <- igraph::set.graph.attribute(G, "layout", lay)
  G <- igraph::set.graph.attribute(G, "width", length(unique(G$layout[,1])) * 2)
  G <- igraph::set.graph.attribute(
    G, "height", (max(G$layout[,2])-min(G$layout[,2])) * 0.4)

  igraph::V(G)$color <- "lightgray"
  igraph::V(G)$shape <- "vrectangle"
  igraph::V(G)$size  <- 16
  igraph::V(G)$size2 <- offs * 1.5
  igraph::V(G)$label <- igraph::V(G)$id
  igraph::E(G)$arrow.size <- 0.5
  igraph::V(G)$rows <- colSums(sweep.result$rows != 0)
  igraph::V(G)$cols <- colSums(sweep.result$columns != 0)
  
  G
}

setMethod("isa", signature(data="matrix"),
          function(data, ...) isa.default(data, ...))

isa.default <- function(data, thr.row=seq(1,3,by=0.5),
                        thr.col=seq(1,3,by=0.5),
                        no.seeds=100, direction=c("updown", "updown")) {

  isa.status("Performing complete ISA work flow", "in")
  
  if (!is.matrix(data)) {
    stop("`data must be a matrix")
  }

  ## Normalize the matrix
  normed.data <- isa.normalize(data)
  
  ## Generate seeds
  row.seeds <- generate.seeds(length=nrow(data), count=no.seeds)

  ## Determine thresholds
  thr.list <- expand.grid(thr.row=thr.row, thr.col=thr.col)
  thr.list <- unlist(apply(thr.list, 1, list), recursive=FALSE)
  
  ## Do the ISA, for all thresholds
  isaresults <- lapply(thr.list, function(x)
                       isa.iterate(normed.data,
                                   row.seeds=row.seeds,
                                   thr.row=x["thr.row"],
                                   thr.col=x["thr.col"],
                                   direction=direction))
  
  ## Make it unique for every threshold combination
  isaresults <- lapply(isaresults, function(x)
                       isa.unique(normed.data, x))

  ## Filter according to robustness
  isaresults <- lapply(isaresults, function(x)
                       isa.filter.robust(data=data,
                                         normed.data=normed.data,
                                         isares=x,
                                         row.seeds=row.seeds))
  
  ## Merge them
  result <- list()
  result$rows <- do.call(cbind, lapply(isaresults, "[[", "rows"))
  result$columns <- do.call(cbind, lapply(isaresults, "[[", "columns"))
  result$seeddata <- do.call(rbind, lapply(isaresults, "[[", "seeddata"))
  result$rundata <- isaresults[[1]]$rundata
  result$rundata$N <- sum(sapply(isaresults, function(x) x$rundata$N))

  ## Another filtering
  result <- isa.unique(normed.data, result)

  isa.status("DONE", "out")
  
  ## We are done
  result
}
