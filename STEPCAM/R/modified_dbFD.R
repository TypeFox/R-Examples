# Reduce the overload of the dbFD function
# Separate into two functions: one to calculate the ordination axes. 
# And the second function uses this axes to calculate FRic, FDiv, and FEve.
# Calculate the ordination axes only once and reuse them during STEPCAM_ABC
ordinationAxes <- function(x, corr = c("sqrt", "cailliez", "lingoes", "none"), 
                           ord = c("podani", "metric"), w, 
                           asym.bin = NULL, messages = FALSE, stand.x = FALSE){
  dist.bin <- NULL                        	
  tol <- .Machine$double.eps
  corr <- match.arg(corr)
  ord <- match.arg(ord)
  if (is.matrix(x) | is.data.frame(x)) {
    is.dist.x <- FALSE
    s.x <- nrow(x)
    t.x <- ncol(x)
    if (is.null(row.names(x))) 
      stop("'x' must have row names.", "\n")
    else x.rn <- row.names(x)
  }
  if (is.vector(x) | is.factor(x)) {
    is.dist.x <- FALSE
    s.x <- length(x)
    t.x <- 1
    if (is.null(names(x))) 
      stop("'x' must have names.", "\n")
    else x.rn <- names(x)
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    is.dist.x <- TRUE
    s.x <- attr(x, "Size")
    t.x <- 1
    if (is.null(attr(x, "Labels"))) 
      stop("'x' must have labels.", "\n")
    else x.rn <- attr(x, "Labels")
  }
  if (missing(w)) 
    w <- rep(1, t.x)/sum(rep(1, t.x))
  if (is.matrix(x) | is.data.frame(x)) {
    x <- data.frame(x)
    if (t.x >= 2) {
      x.class <- sapply(x, data.class)
      if (any(x.class == "character")) 
        x[, x.class == "character"] <- as.factor(x[, x.class == "character"])
      else x <- x
      if (all(x.class == "numeric") & all(!is.na(x))) {
        if (length(unique(w)) == 1) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        else {
          x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
        }
      }
      else {
        x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
      }
    }
    if (t.x == 1) {
      if (is.numeric(x[, 1])) {
        if (all(!is.na(x))) {
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
        }
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          x.s <- apply(x, 2, scale, center = TRUE, scale = stand.x)
          x.dist <- dist(x.s)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          if (messages) 
            cat("Warning: Species with missing trait values have been excluded.", 
                "\n")
        }
      }
      if (is.factor(x[, 1]) | is.character(x[, 1])) {
        if (is.ordered(x[, 1])) 
          x <- x
        else x[, 1] <- as.factor(x[, 1])
        if (any(is.na(x))) {
          pos.NA <- which(is.na(x), arr.ind = TRUE)
          x <- na.omit(x)
          row.excl.ab <- pos.NA[, 1]
          a <- a[, -row.excl.ab]
          x.rn <- x.rn[-pos.NA]
          if (messages) 
            cat("Warning: Species with missing trait values have been excluded.", 
                "\n")
        }
        if (is.ordered(x[, 1])) {
          x.s <- data.frame(rank(x[, 1]))
          names(x.s) <- x.rn
          x.dist <- dist(x.s)
        }
        else {
          x.f <- as.factor(x[, 1])
          x.dummy <- diag(nlevels(x.f))[x.f, ]
          x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
          sequence <- 1:10
          if (all(dist.bin != sequence[any(sequence)])) 
            stop("'dist.bin' must be an integer between 1 and 10.", 
                 "\n")
          x.dist <- dist.binary(x.dummy.df, method = dist.bin)
        }
      }
    }
  }
  if (is.vector(x) & is.numeric(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    x.s <- scale(x, center = T, scale = stand.x)
    x.dist <- dist(x.s)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (is.vector(x) & is.character(x)) {
    x <- as.factor(x)
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    dimnames(x) <- list(x.rn, "Trait")
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", 
           "\n")
    x <- data.frame(x)
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
  }
  if (is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      cat("Warning: Species with missing trait values have been excluded.", 
          "\n")
    }
    else x <- x
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
    x.dist <- gowdis(x, w = w, ord = ord, asym.bin = asym.bin)
  }
  if (is.factor(x) & !is.ordered(x)) {
    if (any(is.na(x))) {
      pos.NA <- which(is.na(x))
      x <- na.omit(x)
      a <- a[, -pos.NA]
      x.rn <- x.rn[-pos.NA]
      if (messages) 
        cat("Warning: Species with missing trait values have been excluded.", 
            "\n")
    }
    else x <- x
    x.dummy <- diag(nlevels(x))[x, ]
    x.dummy.df <- data.frame(x.dummy, row.names = x.rn)
    sequence <- 1:10
    if (all(dist.bin != sequence[any(sequence)])) 
      stop("'dist.bin' must be an integer between 1 and 10.", 
           "\n")
    x.dist <- dist.binary(x.dummy.df, method = dist.bin)
    x <- data.frame(x)
    dimnames(x) <- list(x.rn, "Trait")
  }
  if (class(x)[1] == "dist" | class(x)[1] == "dissimilarity") {
    if (any(is.na(x))) 
      stop("When 'x' is a distance matrix, it cannot have missing values (NA).", 
           "\n")
    x.dist <- x
  }
  if (any(is.na(x.dist))) 
    stop("NA's in the distance matrix.", "\n")
  if (!is.dist.x) {
    no.traits <- apply(x, 1, function(v) length(v[!is.na(v)]))
    if (any(no.traits == 0)) 
      stop("At least one species has no trait data.", "\n")
  } 
  attr(x.dist, "Labels") <- x.rn
  if (is.euclid(x.dist)) 
    x.dist2 <- x.dist
  if (!is.euclid(x.dist)) {
    if (corr == "lingoes") {
      x.dist2 <- lingoes(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Lingoes correction was applied.", 
            "\n")
    }
    if (corr == "cailliez") {
      x.dist2 <- cailliez(x.dist)
      if (messages) 
        cat("Species x species distance matrix was not Euclidean. Cailliez correction was applied.", 
            "\n")
    }
    if (corr == "sqrt") {
      x.dist2 <- sqrt(x.dist)
      if (!is.euclid(x.dist2)) 
        stop("Species x species distance matrix was still is not Euclidean after 'sqrt' correction. Use another correction method.", 
             "\n")
      if (is.euclid(x.dist2)) 
        if (messages) 
          cat("Species x species distance matrix was not Euclidean. 'sqrt' correction was applied.", 
              "\n")
    }
    if (corr == "none") {
      x.dist2 <- quasieuclid(x.dist)
      if (messages) 
        cat("Species x species distance was not Euclidean, but no correction was applied. Only the PCoA axes with positive eigenvalues were kept.", 
            "\n")
    }
  }
  x.pco <- dudi.pco(x.dist2, scannf = FALSE, full = TRUE)
  return(x.pco) # Best is to return the whole object, because dbFD needs sometimes more than the axes alone?!
}


# Strip away radically everything except FRic, FDiv, and FEve
# Because this function needs to be much much faster!
strippedDbFd <- function(x.pco, a, m = "max"){
  tol <- .Machine$double.eps
  a <- as.matrix(a)
  c <- nrow(a) # Number of communities
  traits <- x.pco$li 
  # Number of species per community where traits are present - do I need it later? Yes!
  nb.sp <- numeric(c)
  for (i in 1:c) {
    sp.pres <- which(a[i, ] > 0)
    traits.sp.pres <- traits[sp.pres, , drop = F]
    traits.sp.pres[traits.sp.pres != 0 & abs(traits.sp.pres) < tol] <- 0
    nb.sp[i] <- nrow(unique(traits.sp.pres))
  }
  min.nb.sp <- min(nb.sp)
  m.max <- min.nb.sp - 1
  Warning <- FALSE
  if (min.nb.sp < 3) {
    nb.sp2 <- nb.sp[nb.sp > 2]
    m.max <- min(nb.sp2) - 1    
  }
  else {
    m.max <- m.max
  }
  if(!is.numeric(m)){
    m <- m.max
  }
  
  if (m < x.pco$nf){ # If there is a community with less species than ordination axes
    traits.FRic <- x.pco$li[, 1:m] # Use only the first m axes
  }
  if (m >= x.pco$nf){
    traits.FRic <- x.pco$li
  }
  nbsp <- rep(NA, c)
  names(nbsp) <- row.names(a)
  FRic <- nbsp
  FEve <- nbsp
  FDiv <- nbsp
  
  for (i in 1:c) { # For each community
    sppres <- which(a[i, ] > 0) # Present species
    S <- length(sppres)
    nbsp[i] <- S
    tr <- data.frame(traits[sppres, ]) # Axes coordinates of the present species
    tr.FRic <- data.frame(traits.FRic[sppres, ])
    # Will I need relative abundances of the species?
    ab <- as.matrix(a[i, sppres])
    abundrel <- ab/sum(ab)
    
    if (ncol(tr.FRic) > 1 & nb.sp[i] >= 3) { # If there are more than 3 species present
      if (Warning) 
        thresh <- 4
      if (!Warning) 
        thresh <- 3
      if (nb.sp[i] >= thresh) {
        convhull <- convhulln(tr.FRic, "FA")
        FRic[i] <- convhull$vol
      }
    }
    if (ncol(tr.FRic) == 1) {
      tr.range <- range(tr.FRic[, 1])
      t.range <- tr.range[2] - tr.range[1]
      FRic[i] <- t.range
    }
    
    if (nb.sp[i] >= 3) {
      tr.dist <- dist(tr) # pair-wise distance of ordination coordinates
      linkmst <- mst(tr.dist)
      mstvect <- as.dist(linkmst)
      abund2 <- matrix(0, nrow = S, ncol = S)
      for (q in 1:S) for (r in 1:S) abund2[q, r] <- abundrel[q] + abundrel[r]
      abund2vect <- as.dist(abund2)
      EW <- rep(0, S - 1)
      flag <- 1
      # Loops should be avoided
      for (m in 1:((S - 1) * S/2)) {
        if (mstvect[m] != 0) {
          EW[flag] <- tr.dist[m]/(abund2vect[m])
          flag <- flag + 1
        }
      }
      minPEW <- rep(0, S - 1)
      OdSmO <- 1/(S - 1)
      for (l in 1:(S - 1)) minPEW[l] <- min((EW[l]/sum(EW)), OdSmO)
      FEve[i] <- ((sum(minPEW)) - OdSmO)/(1 - OdSmO)
    }
    if (ncol(tr.FRic) > 1 & nb.sp[i] >= 3) {
      vert0 <- convhulln(tr.FRic, "Fx TO 'vert.txt'")
      vert1 <- scan("vert.txt", quiet = T)
      vert2 <- vert1 + 1
      vertices <- vert2[-1]
      trvertices <- tr.FRic[vertices, ]
      baryv <- apply(trvertices, 2, mean)
      distbaryv <- rep(0, S)
      for (j in 1:S) distbaryv[j] <- (sum((tr.FRic[j, ] - baryv)^2))^0.5
      meandB <- mean(distbaryv)
      devdB <- distbaryv - meandB
      abdev2 <- abundrel * devdB
      ababsdev2 <- abundrel * abs(devdB)
      FDiv[i] <- (sum(abdev2) + meandB)/(sum(ababsdev2) + meandB)
    }    
  }
  res <- list()
  res$nbsp <- nbsp
  res$sing.sp <- nb.sp
  res$FRic <- FRic
  res$FEve <- FEve
  res$FDiv <- FDiv
  return(res)  
}


