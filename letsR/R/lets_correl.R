#' Compute correlogram based on the Moran's I index
#' 
#' @author Bruno Vilela, Fabricio Villalobos, Lucas Jardim & Jose Alexandre Diniz-Filho
#' 
#' @description Computes the Moran's I correlogram of a single or multiple variables.
#'
#' @param x A single numeric variable in vector format or multiple variables in matrix format 
#' (as columns). 
#' @param y A distance matrix of class \code{matrix} or \code{dist}.
#' @param z The number of distance classes to use in the correlogram.
#' @param equidistant Logical, if \code{TRUE} the classes will be equidistant.
#' If \code{FALSE} the classes will have equal number of observations.
#' @param plot Logical, if \code{TRUE} the correlogram will be ploted. 
#' 
#' @return Returns a matrix with the Moran's I Observed value, Confidence Interval (95%) 
#' and Expected value. Also the p value of the randomization test, the mean distance 
#' between classes, and the number of observations.   
#' quase tudo 
#' @references Sokal, R.R. & Oden, N.L. (1978) Spatial autocorrelation in biology.
#' 1. Methodology. Biological Journal of the Linnean Society, 10, 199-228.
#' @references Sokal, R.R. & Oden, N.L. (1978) Spatial autocorrelation in biology.
#' 2. Some biological implications and four applications of evolutionary and
#' ecological interest. Biological Journal of the Linnean Society, 10, 229-249.
#' 
#' @examples \dontrun{
#' data(PAM)
#' data(IUCN)
#' 
#' # Spatial autocorrelation in description year (species level)
#' midpoint <- lets.midpoint(PAM)
#' distan <- lets.distmat(midpoint[, 2:3])
#' moran <- lets.correl(IUCN$Description, distan, 12,
#'                      equidistant = FALSE, 
#'                      plot = TRUE)
#'                      
#' }
#' 
#' @export




lets.correl <- function(x, y, z, equidistant = FALSE,
                        plot = TRUE) {
  
  # Absent values in x (better error message)
  if (any(is.na(x))) {
    stop("Missing values in x argument")
  }
  
  # Allow dist classes
  if (class(y) == "dist") {
    y <- as.matrix(y)
  }
  
  # Absent values in y (better error message)
  if (any(is.na(y))) {
    stop("Missing values in y argument")
  }
  
  
  # Check if it has one
  if (is.matrix(x)) {
    if (ncol(x) == 1) {
      x <- as.vector(x)
    }
  }
  
  
  if (is.vector(x)) {
    return1 <- .br.correlogram(x, y, z, equidistant, plot)
    # Warning for removing some classes
    if (nrow(return1) < z) {
      warning(paste("Some of the distance classes", 
                    "were removed due to small number",
                    "of occurrences on it"))
    }
    return(return1)
  }
  
  if (!is.vector(x)) {
    n <- ncol(x)
    parcial <- list()
    
    for(i in 1:n) {
      parcial1 <- .br.correlogram(x[, i], y, z, plot = FALSE)
      parcial[[i]] <- parcial1[, c(1, 3, 5, 6)]
    }
    media <- apply(simplify2array(parcial), 1:2, mean)
    desvio <- apply(simplify2array(parcial), 1:2, sd) 
    resu <- cbind(media[, 1], desvio[, 1], media[, 2],
                  media[, 3], media[, 4] )
    colnames(resu) <- c("Observed", "Standard_Deviation",
                        "Expected_value", "Mean_Distance",
                        "Count")
    
    if (plot) {
      .plotcorrel(plot1 = resu[, 4],
                  plot2 = resu[, 1],
                  plot3 = resu[, 2],
                  plot4 = resu[, 3],
                  z)
    }
    
    if (nrow(resu) < z) {
      warning(paste("Some of the distance classes", 
                    "were removed due to few number",
                    "of occurrence on it"))
    }
    
    return(resu)
  }
}





############
.br.correlogram <- function(x, y, z, equidistant = FALSE,
                            plot = TRUE) {
  y3 <- y
  diag(y3) <- NA
  y2 <- y
  z2 <- 1 / z
  
  if (!equidistant) {
    quant <- quantile(y3, probs = seq(0, 1, z2),
                      na.rm = TRUE)
  }
  if (equidistant) {
    quant <- seq(min(y), max(y), ((max(y) - min(y)) / z))
  }
  
  quant <- as.vector(quant)
  n <- length(quant)
  ob <- rep(NA, (n - 1))
  CI <- rep(NA, (n - 1))
  ex <- rep(NA, (n - 1))
  dist_cl <- rep(NA, (n - 1))
  p <- rep(NA, (n - 1))
  count <- rep(NA, (n - 1))
  
  for(i in 1:(n - 1)) {
    if (i > 1) {
      pos <- (y > quant[i] & y <= quant[i + 1])
      diag(pos) <- FALSE
      count[i] <- sum(pos)
    }
    
    if (i == 1) {
      pos <- (y >= quant[i] & y <= quant[i + 1])
      diag(pos) <- FALSE
      count[i] <- sum(pos)
    }
    dist_cl [i] <- mean(c(quant[i], quant[i + 1]))
    
    if (count[i] > 0) {
      y2[pos] <- 1
      y2[!pos] <- 0
      m <- .br.moran(y2, x)
      ob[i] <- m$observed
      CI[i] <- m$ci
      ex[i] <- m$expected
      dist_cl[i] <- mean(c(quant[i], quant[i + 1]))
      p[i] <- m$p.value
    }
  }
  
  if (all(is.na(ob))) {
    stop(paste("None of the distance classes sets",
               "has enough sample size to calculate",
               "Moran's I. Set less classes or increase",
               "the sample size"))
    
  }
  
  resu <- cbind(ob, CI, ex, p, dist_cl, count)
  colnames(resu) <- c("Observed", "Confidence_Interval_(95%)",
                      "Expected_value", "p_value",
                      "Mean_Distance", "Count")
  resu <- resu[!is.na(resu[, 1]), , drop = FALSE]
  
  if (plot) { 
    
    .plotcorrel(plot1 = resu[, 5],
                plot2 = resu[, 1],
                plot3 = resu[, 2],
                plot4 = resu[, 3],
                z)
  }
  
  return(resu)
}


############

.br.moran <- function(w, y) {
  
  n <- sum(ifelse(rowSums(w) > 0, 1, 0))
  
  if (n > 3) {
    z <- y - mean(y)
    soma <- n * (sum(w * (z %o% z)))
    divi <- sum(w) * sum((z ^ 2))
    ob <- soma / divi
    es <- -1 / (n - 1)
    S1 <-  0.5 * sum((w + t(w)) ^ 2)
    S2 <- sum((apply(w, 1, sum) + apply(w, 2, sum)) ^ 2)
    k <- n * sum(z ^ 4) / ((sum(z ^ 2)) ^ 2)
    s.sq <- sum(w) ^ 2
    calc1 <- (n ^ 2 - 3 * n + 3) * S1 - n * S2 + 3 * s.sq
    calc2 <- n * (n - 1) * S1 - 2 * n * S2 + 6 * s.sq
    calc3 <- (n - 1) * (n - 2) * (n - 3) * s.sq
    sdi0 <- (n * calc1 - k * calc2) / calc3 - 1/((n - 1) ^ 2)
    if (sdi0 < 0) {
      sdi0 <- 0
    }
    sdi <- sqrt(sdi0)
    ci <- sdi * 1.96
    pv <- pnorm(ob, mean = es, sd = sdi)
    if (ob <= es) {
      pv <- 2 * pv
    } else {
      pv <- 2 * (1 - pv)
    }
  } else {
    ob <- NA
    es <- NA
    ci <- NA
    pv <- NA
  } 
  return(list("observed" = ob, "expected" = es, "ci" = ci,
              "p.value" = pv))
}


### Plot

.plotcorrel <- function(plot1, plot2, plot3, plot4, z) {
  epsilon <- max(plot1) / (14 * z)
  up <- plot2 + plot3
  low <- plot2 - plot3  
  plot(x = plot1, y = plot2, bty = "l", 
       ylab = "Moran's I", xlab = "Distance",
       type = "l", lty = 3, 
       ylim = c(min(low) - 0.2, max(up) + 0.2))
  abline(h = mean(plot4))
  points(x = plot1, y = plot2,
         pch = 20, cex = 1.5)
  segments(plot1, low,
           plot1, up)
  segments(plot1 - epsilon, up,
           plot1 + epsilon, up)
  segments(plot1 - epsilon, low,
           plot1 + epsilon, low)
  invisible(NULL)
}