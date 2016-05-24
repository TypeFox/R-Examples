# Z score section ---------------------------------


# Calculates the f score.
#
# Calculates the f score, defined as the difference between f max and f min.
#
# @param xValues the coordinates of the points describing the curve
# @param yValues the coordinates of the points describing the curve
#
# @return A vector contening the distance of each point to the diagonale
#
#
# @examples
# data(yer109c)
# cdf <- cumsum(yer109c$mean_type1) - cumsum(yer109c$mean_type2)
# plot(yer109c$pos, cdf)
# yCrypticRNAs:::fscore(yer109c$pos, cdf)
fscore <- function(xValues, yValues) {
  distances <- .perpendicular_distance(xValues, yValues)
  max(distances) - abs(min(distances))
}


# Calculates the perpendicular distance.
#
# Calculate the perpendicular distance from the curve to the diagonal
#
# @param xValues the coordinates of the points describing the curve
# @param yValues the coordinates of the points describing the curve
#
# @return A vector contening the distance of each point to the diagonale
#
# @examples
# data(yer109c)
# cdf <- cumsum(yer109c$mean_type1) - cumsum(yer109c$mean_type2)
#
# par(mfrow = c(1,2))
# plot(yer109c$pos, cdf, type = "l")
#
# d <- yCrypticRNAs:::.perpendicular_distance (yer109c$pos, cdf)
# plot(yer109c$original_pos, d, type = "l")
.perpendicular_distance <- function(xValues, yValues) {
  df <- data.frame(pos = xValues, cdf = yValues)

  x1 <- min(df$pos)
  x2 <- max(df$pos)
  y1 <- df$cdf[which(df$pos == x1)]
  y2 <- df$cdf[which(df$pos == x2)]

  A <- y1 - y2
  B <- x2 - x1
  C <- (x1 * y2) - (x2 * y1)

  equation <- paste(A, "x +", B, "y + ",C, " = 0")
  denominator <- sqrt(A ^ 2 + B ^ 2)
  (A * df$pos + B * df$cdf + C) / denominator
}


# ' calculate the CDF in the mutant - the CDF in the wild-type
.difference_of_cdf <- function(data) {
  cumsum(data[[1]]) - cumsum(data[[2]])
}


# A class to represent a cryptic score
#
# @param geneInfo A vector containing fields to represent the gene.
# @param cryptic_score The value of the cryptic score
#  comparing the mutant and the wild-type values
# @param controls The value of the cryptic score
# comparing replicates in the wild-type and in the mutant.
# @param method The method used to calculate the cryptic score
#
# @return An object of class 'CrypticScore' with the following components:
#        \item{geneAnnotation}{An object of class \code\link{annotationsSet}}
#                              containing the information on the gene.
#        \item{crypticScore}{Cryptic score obtained by comparing
#                            type2 data to type1 data}
#        \item{controls}{A list containing the scores obtained by
#                        comparing replicates of each type of data}
#        \item{method}{The method used.}
CrypticScore <- function(geneInfo, cryptic_score, controls, method) {
  cs <- list(
    geneAnnoation  = geneInfo,
    crypticScore = cryptic_score,
    controls = controls,
    method = method
  )

  ## Set the name for the class
  class(cs) <- append(class(cs),"CrypticScore")
  return(cs)
}


#' @title Z score
#'
#' @description Calculates the cryptic score (z score) using the probabilistic method.
#'
#' @param geneCoverage object of type geneCoverage containing the
#'        RNA-seq coverage values for one gene.
#' @param iterations A number indicating the number of iterations. Default = 10000.
#'
#' @return An object of class 'CrypticScore' with the following components:
#'       \item{geneAnnotation}{An object of class \code{\link{annotationsSet}}
#'                              containing the information on the gene.}
#'        \item{crypticScore}{Cryptic score obtained by comparing
#'                            type2 data to type1 data}
#'        \item{controls}{A list containing the scores obtained by
#'                        comparing replicates of each type of data}
#'        \item{method}{The method used.}
#'
#' @export
#' @import data.table
#'
#' @examples
#'
#' data(yer109c)
#' zscore_score(geneCoverage = yer109c, iterations = 1000)
zscore_score <- function(geneCoverage, iterations = 10000) {
  score <- .zscore(
    data = list(
      pos = geneCoverage$pos,
      wt = geneCoverage$mean_type1,
      mut = geneCoverage$mean_type2
    ),
    iterations = iterations
  )

  controls_names <- NULL
  if(is.na(score)){
    self_wt <- self_mut <- NA
  }else {
    type1 <- geneCoverage$samples$type1
    type2 <- geneCoverage$samples$type2

    if(length(type1) > 1){
      self_wt <- .zscore(
        data = list(
          pos = geneCoverage$pos,
          wt = geneCoverage[type1][[1]],
          mut = geneCoverage[type1][[2]]
        ),
        iterations = iterations
      )
      controls_names <- c(controls_names, paste0(type1[[2]],"/", type1[[1]]))
    }else{
      self_wt <- NA
      controls_names <- c(controls_names, type1)
    }

    if(length(type2) > 1){
      self_mut <- .zscore(
        data = list(
          pos = geneCoverage$pos,
          wt = geneCoverage[type2][[1]],
          mut = geneCoverage[type2][[2]]
        ),
        iterations = iterations
      )
      controls_names <- c(controls_names, paste0(type2[[2]], "/", type2[[1]]))
    }else{
      self_mut <- NA
      controls_names <- c(controls_names, type2)
    }
  }

  controls <- list(self_wt, self_mut)
  names(controls) <- controls_names
  CrypticScore(
    geneInfo = geneCoverage$gene_information,
    cryptic_score = score,
    controls = controls,
    method = "probabilistic method"
  )
}


## Simulations of values and calculation of the F score for each simulation
.simulate_fscore_values <- function(simulationIndex, values) {
  data <- lapply(values[c(2:3)], sample)
  cdf <- .difference_of_cdf(data)
  rm(data)
  res <- fscore(xValues = values[[1]], yValues = cdf)
  rm(cdf)
  return(res)
}

## calculate the z score between two samples.
#' @importFrom stats sd
.zscore <- function(data, iterations){

  if(sum(data$wt) == 0 || sum(data$mut)==0){
    return (NA)
  }

  observed_f <-
    .perpendicular_distance(xValues = data$pos, yValues <-
                             .difference_of_cdf(data[c(2:3)]))

  observed_f_value <- fscore(data$pos, yValues)

  if (length(which.max(observed_f)) > 1) {
    return (NA)
  }

  gc()
  simulated_values <-
    unlist(parallel::mclapply(seq_len(iterations), .simulate_fscore_values, values = data))

  (observed_f_value - mean(simulated_values)) / sd(simulated_values)
}


## Ratio score sectio -----------------------------------------------


#' @title 3'/5' ratio score
#'
#' @description Calculates the cryptic score (ratio score) using
#' the 3'/5' ratio method as described in Cheung et al., 2008.
#'
#' @param geneCoverage object of type geneCoverage containing the
#'        coverage values for all samples.
#' @param windowLength A integer indicating the length of the window to use at each end of the gene.
#'
#' @return An object of class 'CrypticScore' with the following components:
#'       \item{geneAnnotation}{An object of class \code{\link{annotationsSet}}
#'                              containing the information on the gene.}
#'        \item{crypticScore}{Cryptic score obtained by comparing
#'                            type2 data to type1 data}
#'        \item{controls}{A list containing the scores obtained by
#'                        comparing replicates of each type of data}
#'        \item{method}{The method used.}
#' @export
#' @examples
#' data(yer109c)
#' ratio_score(yer109c)
#'
ratio_score <- function(geneCoverage, windowLength = 100) {

  type1 <- geneCoverage$samples$type1
  type2 <- geneCoverage$samples$type2

  ratios <- lapply(X = geneCoverage[c(type1, type2)],
                   FUN = .ratio, windowLength = windowLength)

  score_f <- function(x,y){
    res <- NULL
    y <- unlist(y)
    x <- unlist(x)
    for (i in 1:length(y)){
      label <- paste0(names(x), "/", names(y[i]))
      z <- x/y[i]
      names(z) <- label
      res <- c(res, z)
    }
    res
  }

  self_f <- function(x){
    if(length(x) == 1){
      res <- NA
      names(res) <- names(x)
      return (res)
    }
    res <- NULL

    names <- names(x)
    for( i in 1:length(x)){
      z <- unlist(x[-i])
      label <- paste0(names[i], "/", names(z))
      z <- x[[i]]/z
      names(z) <- label
      res <- c(res, z)
    }
    res
  }

  self_wt <- self_f(ratios[type1])
  self_mut <- self_f(ratios[type2])
  scores <- score_f(ratios[type2], ratios[type1])

  score <- mean(unlist(ratios[type2])) / mean(unlist(ratios[type1]))

  CrypticScore(
    geneInfo = geneCoverage$gene_information,
    cryptic_score = list(mean = score, scores = scores),
    controls = c(self_wt, self_mut),
    method = "3'/5' ratio method"
  )
}

.ratio <- function(data, windowLength) {
  if (length(data) < windowLength * 2) {
    return (NA)
  }

  gene_length <- length(data)
  five_prime_region <- mean(data[seq_len(windowLength)])
  three_prime_region <-
    mean(data[(gene_length - windowLength):gene_length])
  (three_prime_region + 1) / (five_prime_region + 1)
}


# Enrichment score section -----------------------------------

#' @title 3' enrichment score
#'
#' @description Calculates the cryptic score (3' enrichment score) using the 3'/5' ratio method as described in DeGennaro et al., 2013.
#'
#' @param geneCoverage object of type geneCoverage containing the
#'        coverage values for all samples.
#'
#' @return An object of class 'CrypticScore' with the following components:
#'       \item{geneAnnotation}{An object of class \code{\link{annotationsSet}}
#'                              containing the information on the gene.}
#'        \item{crypticScore}{Cryptic score obtained by comparing
#'                            type2 data to type1 data}
#'        \item{controls}{A list containing the scores obtained by
#'                        comparing replicates of each type of data}
#'        \item{method}{The method used.}
#'
#' @export
#' @examples
#' data(yer109c)
#' enrichment_score (yer109c)
#'
enrichment_score <- function(geneCoverage){

  score <- .three_prime_enrichment(list(geneCoverage$mean_type1, geneCoverage$mean_type2))

  type1 <- geneCoverage$samples$type1
  type2 <- geneCoverage$samples$type2

  controls_names <- NULL
  if (length(type1) > 1) {
    score_wt <- .three_prime_enrichment(geneCoverage[type1][c(1:2)])
    controls_names <- c(controls_names, paste0(type1[[1]], "/", type1[[2]]))
  }else{
    score_wt = NA
    controls_names <- c(controls_names, type1)
  }
  if (length(type2) > 1) {
    score_mut <- .three_prime_enrichment(geneCoverage[type2][c(1:2)])
    controls_names <- c(controls_names, paste0(type2[[1]], "/", type2[[2]]))
  }else{
    score_mut = NA
    controls_names <- c(controls_names, type2)
  }

  controls <- list(score_wt, score_mut)
  names(controls) <- controls_names
  CrypticScore(
    geneInfo = geneCoverage$gene_information,
    cryptic_score = score,
    controls = controls,
    method = " 3' enrichment method"
  )
}

.area_under_diagonale <- function(x, data) {
  l <- sapply(data, length)
  vect1 <- data[[1]]
  vect2 <- data[[2]]

  x1 <- vect1[1]
  x2 <- vect1[l[1]]
  y1 <- vect2[1]
  y2 <- vect2[l[2]]

  A <- y1 - y2
  B <- x2 - x1
  C <- (x1 * y2) - (x2 * y1)

  equation <- paste(A, "x +", B, "y + ",C, " = 0")
  - (A * x + C) / B
}

.three_prime_enrichment <- function(data) {
  cdf <- lapply(data, function(values) {
    rev(cumsum(rev(values)))
  })

  if (sum(cdf[[1]]) == 0 || sum(cdf[[2]]) == 0) {
    score <- NA
  }else if (mean(cdf[[1]]) == cdf[[1]][1]) {
    score <- NA
  }else {
    areaUnderDiagonale <- stats::integrate(f = .area_under_diagonale,
                                           lower = min(cdf[[1]]),
                                           upper = max(cdf[[1]]),
                                           data = cdf)$value
    areaUnderCurve <- MESS::auc(cdf[[1]], cdf[[2]])
    score <- areaUnderCurve / areaUnderDiagonale
  }
  score
}

# Genome-wide section -----------------
#' Genome-wide cryptic scores.
#'
#' Genome-wide calculation of cryptic scores with a specified method.
#'
#' @param coverageDataSet A coverageDataSet containing the coverage values for all genes.
#' @param method A caracter vector indicating the method to be used.
#' @param outfile A character vector indicating the output file name.
#' @param introns An objet of type \code{\link{annotationsSet}}
#'        containing the annotations of the intronic regions.
#'        Note: The introns must have same name as the gene they
#'        are associated with.
#' @param windowLength If the \code{method} is "ratio", specify the length of the window to use
#'        at each end of the gene.
#' @param iterations If the \code{method} is "probabilistic", specify the number of iterations.
#'
#' @export
#' @examples
#' data("rna_seq_signals")
#' #genome_wide_scores(rna_seq_signals, "ratio", "ratio.txt")
#' genome_wide_scores(rna_seq_signals, "enrichment", "enrichment.txt")
#' #genome_wide_scores(rna_seq_signals, "probabilistic", "probabilistic.txt")
#' unlink(c("ratio.txt", "enrichment.txt", "probabilistic.txt"))
genome_wide_scores <- function (coverageDataSet, method = c("ratio", "enrichment",  "probabilistic"), outfile,
                                introns = NULL, windowLength = 100, iterations = 10000){

  if (!is.element("coverageDataSet", class(coverageDataSet))) {
    stop("You must provide an object of class coverageDataSet!")
  }

  method <- match.arg(method)
  all_genes <- unique(coverageDataSet$data$name)


  analyse <- function(gene, first = FALSE){
    cov <- gene_coverage(coverageDataSet, gene, introns)
    col.names = ifelse(first, TRUE, FALSE)
    append = ifelse(first, FALSE, TRUE)

    if(method == "ratio"){
      cs <- ratio_score(cov, windowLength)
      .writeFile(
        c(cs$geneAnnoation, cryptic_score = cs$crypticScore$mean,
          cs$crypticScore$scores, controls = cs$controls),
        outfile, append = append, col.names
      )
    }else{
      if( method == "enrichment"){
        cs <- enrichment_score(cov)
      }else{
        cs <- zscore_score(cov, iterations)
      }
      .writeFile(
        c(cs$geneAnnoation, cryptic_score = cs$crypticScore,
          controls = cs$controls),
        outfile, append = append, col.names
      )
    }
  }

  first_gene <- all_genes[1]
  if(method == "ratio" || method == "enrichment"){
    analyse(first_gene, TRUE)
    invisible(parallel::mclapply(all_genes[2:length(all_genes)], analyse))
  } else{
    analyse(first_gene, TRUE)
    invisible(lapply(all_genes[2:length(all_genes)], analyse))
  }
}


#   print generic functions ---------------------------------------------------------------------

#' @export
print.CrypticScore <- function(x, ...) {
  cat(paste0(
    "Cryptic score for ", x$gene$name,
    " using the ", x$method, ". \n\n"
  ))

  if (length(x$crypticScore) == 1) {
    cat(paste("cryptic score = ", x$crypticScore, "\n"))
  }else{
    cat(paste("cryptic score = ", x$crypticScore$mean, " [ "))
    values <- x$crypticScore[[2]]
    for (i in 1:length(values)) {
      cat (paste0(round(values[[i]], digits = 2), " (", names(values)[i], ") "))
    }
    cat("] \n")
  }

  cat("controls = ")
  for (i in 1:length(x$controls)) {
    cat (paste0(
      round(x$controls[[i]], digits = 4),
      " (", names(x$controls)[i], ") "
    ))
  }
}
