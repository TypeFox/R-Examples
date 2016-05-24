pamr.makeclasses <- function(data,  sort.by.class = FALSE, ...) {
#  require(cluster)
  as.matrix.dist <- function (x)  {
    size <- attr(x, "Size")
    df <- matrix(0, size, size)
    df[row(df) > col(df)] <- x
    df <- df + t(df)
    labels <- attr(x, "Labels")
    dimnames(df) <- if (is.null(labels)) 
      list(1:size, 1:size)
    else list(labels, labels)
    df
  }
  as.dist <- function (m, diag = FALSE, upper = FALSE) {
    m <- as.matrix(m)
    retval <- m[row(m) > col(m)]
    attributes(retval) <- NULL
    if (!is.null(rownames(m))) 
      attr(retval, "Labels") <- rownames(m)
    else if (!is.null(colnames(m))) 
      attr(retval, "Labels") <- colnames(m)
    attr(retval, "Size") <- nrow(m)
    attr(retval, "Diag") <- diag
    attr(retval, "Upper") <- upper
    attr(retval, "call") <- match.call()
    class(retval) <- "dist"
    retval
  }
  
  if(!is.null(data$samplelabels)) {
    labs <- data$samplelabels
  }
  if(!is.null(data$samplelabels) & !is.null(data$y)) {
    labs <- paste(data$y, labs)
  }
  if(is.null(data$samplelabels)) {
    labs <- 1:ncol(data$x)
  }
  par(col = 1, cex = 1)
  d <- dist(t(data$x))
  dd <- as.matrix.dist(d)
  if(sort.by.class) {
    tt <- table(data$y)
    nc <- length(tt)
    for(i in 1:nc) {
      o <- data$y == names(tt[i])
      d1 <- max(dd[o, o])
      d2 <- min(dd[o, !o])
      fac <- ((0.2 + (0.7 * i)/nc) * d2)/d1
      dd[o, o] <- dd[o, o] * fac
    }
  }
  hc <- hclust(as.dist(dd), ...)
  plot(hc, labels = labs)
  aa <- vector("list", 100)
  go <- TRUE
  i <- 0
  while(go & i < 100) {
    go <- FALSE
    i <- i + 1
    print(c("Identify class", i))
    par(pch = as.character(i), col = 4)
    aa[[i]] <- locator(type = "p")
    if(!is.null(aa[[i]])) {
      go <- TRUE
    }
  }
  nclus <- i - 1
  res <- vector("list", nclus)
  for(i in 1:nclus) {
    res[i] <- aa[i]
  }
  hdelta <- 1
  clus <- vector("list", nclus)
  for(j in 1:nclus) {
    for(jj in 1:length(res[[j]]$x)) {
      r <- c(res[[j]]$x[jj], res[[j]]$y[jj])
      d <- abs(hc$hei - r[2])
      o <- rank(d)
      ncomp <- 5
      oo <- (1:length(o))[o < ncomp + 1 & d < hdelta]
      if(length(oo) == 0) {
        stop(
             "1 Ambigious selection; try pamr.makeclasses again"
             )
      }
      ncomp2 <- length(oo)
      good <- rep(FALSE, ncomp2)
      ordpos <- match(1:length(hc$ord), hc$ord)
      nodes <- vector("list", ncomp2)
      for(ii in 1:ncomp2) {
        ooo <- descendants(hc$mer, oo[ii])[[2]]
        o4 <- as.vector(hc$mer[ooo,  ])
        nodes[[ii]] <- -1 * o4[o4 < 0]
        op <- ordpos[nodes[[ii]]]
        if(r[1] > min(op) & r[1] < max(op)) {
          good[ii] <- TRUE
        }
      }
                                        #browser()
      if(sum(good) != 1) {
        stop(
             "2 Ambigious selection; try pamr.makeclasses again"
             )
      }
                                        #browser()
      ii2 <- (1:ncomp2)[good]
      clus[[j]] <- c(clus[[j]], nodes[[ii2]])
    }
  }
  newy <- rep(NA, ncol(data$x))
  temp <- NULL
  for(i in 1:nclus) {
    clus[[i]] <- unique(clus[[i]])
  }
  for(i in 1:nclus) {
    temp <- c(temp, clus[[i]])
  }
  if(length(unique(temp)) < length(temp)) {
    stop("Clusters overlap; try pamr.makeclasses again")
  }
  for(i in 1:nclus) {
    newy[clus[[i]]] <- i
  }
  labs2 <- as.character(newy)
  labs2[labs2 == "NA"] <- ""
  par(col = 1, cex = 1)
  plot(hc, labels = labs2)
  return(as.factor(newy))
}

