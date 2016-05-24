##' Calculate one of several measures of mating compatibility.
##'
##' @title Make potentials object--mating type compatibility
##' @param scene a matingScene object
##' @param method either "si_echinacea" or "dioecious" see details for
##' further description
##' @param subject whether you want pair, individual, population, or all.
##' Specifying more than one is allowed.
##' @param averageType whether to calculate individual and population proximity
##' using the mean or median
##' @return A potentials object containing one more more of the following, depending the
##' input for \code{subject}: \cr
##' If \code{subject} is "population" the return list will contain a numeric
##' value that has a range depending on the \code{method}. If
##' \code{subject} is "pair" the return list will contain a matrix
##' with all pairwise compatibilities. If \code{subject} is "individual"
##' the return list will contain a dataframe with a column containing IDs and
##' a column containing compatibility averages. If \code{subject} is "all"
##' the return list will contain all three of the items above.
##' @details When \code{method} is "si_echinacea" compatibility will be
##' calculated as sporophytic self incompatible (si) in the same manner as Echinacea
##' (and many other plants). For two individuals, they are incompatible if
##' they share any S alleles (columns s1 and s2) and they compatible otherwise.
##' When \code{method} is "dioecious" it is assumed that the column s1 will
##' contain either a 1 or 2 depending on the individual's sex. Thus, when
##' comparing two individuals, they are compatible if s1 of the first != s1
##' of the second, and s2 is ignored.
##' @author Danny Hanson
##' @examples
##' pop <- simulateScene()
##' compatibility(pop, "si_echinacea")
compatibility <- function(scene, method, subject = "all",
                          averageType = "mean"){

  method <- match.arg(method, c("si_echinacea", "dioecious"))
  subject <- match.arg(subject, c("population", "pairwise",
                                  "individual", "all"), several.ok = T)
  averageType <- match.arg(averageType, c("mean", "median"))

  if (averageType == "mean") {
    average <- mean
  } else if (averageType == "median") {
    average <- median
  }

  if (is.list(scene) & !is.data.frame(scene)) {
    potential <- lapply(scene, compatibility, method, subject, averageType)
  } else {

    if (method == "si_echinacea") {
      pairCompat <- pair_si_ech(scene$s1, scene$s2)
      attr(pairCompat, "idOrder") <- scene$id

      indCompat <- data.frame(id = scene$id, compatibility = -1)
      if (averageType == "mean") {
        indCompat$compatibility <- rowMeans(pairCompat)
      } else if (averageType == "median") {
        indCompat$compatibility <- row_medians(pairCompat)
      }

      popCompat <- average(indCompat$compatibility)
    } else if (method == "dioecious") {
      pairCompat <- pair_dioecious(scene$s1)
      attr(pairCompat, "idOrder") <- scene$id
      
      indCompat <- data.frame(id = scene$id, compatibility = -1)
      if (averageType == "mean") {
        indCompat$compatibility <- rowMeans(pairCompat)
      } else if (averageType == "median") {
        indCompat$compatibility <- row_medians(pairCompat)
      }
      
      popCompat <- average(indCompat$compatibility)
    }

    # return
    potential <- list()
    if ("population" %in% subject) {
      potential$pop <- popCompat
    }
    if ("individual" %in% subject) {
      potential$ind <- indCompat
    }
    if ("pairwise" %in% subject) {
      potential$pair <- pairCompat
    }
    if ("all" %in% subject) {
      potential$pop <- popCompat
      potential$ind <- indCompat
      potential$pair <- pairCompat

    }
    attr(potential, "t") <- FALSE
    attr(potential, "s") <- FALSE
    attr(potential, "c") <- TRUE
    potential

  }
}
