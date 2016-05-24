#' @rdname popStructStat
#' @export
#' 
statGst <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  if(is.null(strata.mat)) g <- g[, , strataNames(g)]
  if(ploidy(g) < 2 | nStrata(g) == 1) {
    return(list(
      stat.name = "Gst", 
      result = c(estimate = NA, p.val = NA),
      null.dist = NULL
    ))
  }
  
  strata.mat <- .checkStrataMat(strata.mat, g, nrep)
  
  result <- statGst_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    strata.mat, ploidy(g)
  )
  
  return(.formatResult(result, "Gst", keep.null))
  
#   hets <- Hstats(g, strata)
#   Hs <- mean(hets["Hs", ], na.rm = TRUE)
#   Ht <- mean(hets["Ht", ], na.rm = TRUE)
#   est <- 1 - (Hs / Ht) 
#   if(is.nan(est)) est <- NA
#   
#   names(est) <- "Gst"
#   est
}


#' @rdname popStructStat
#' @export
#' 
statGstPrime <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE,
                         prime.type = c("nei", "hedrick"), ...) { 
  if(is.null(strata.mat)) g <- g[, , strataNames(g)]
  if(ploidy(g) < 2 | nStrata(g) == 1) {
    return(list(
      stat.name = "G'st", 
      result = c(estimate = NA, p.val = NA),
      null.dist = NULL
    ))
  }
  
  strata.mat <- .checkStrataMat(strata.mat, g, nrep)
  
  result <- statGstPrime_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    strata.mat, ploidy(g),
    switch( match.arg(prime.type), nei = 0, hedrick = 1)
  )
  
  return(.formatResult(result, "G'st", keep.null))
  
#   hets <- Hstats(g, strata)
#   Hs <- mean(hets["Hs", ], na.rm = TRUE)
#   Ht <- mean(hets["Ht", ], na.rm = TRUE)
#   n <- if(is.null(strata)) nStrata(g) else length(unique(strata))
#   est <- if(prime.type == "nei") {
#     (n * (Ht - Hs)) / ((n * Ht) - Hs) 
#   } else {
#     gst.max <- ((n - 1) * (1 - Hs)) / (n - 1 + Hs)
#     (1 - (Hs / Ht)) / gst.max
#   }
#   if(is.nan(est)) est <- NA
  
#   names(est) <- "G'st"
#   est
}


#' @rdname popStructStat
#' @export
#' 
statGstDblPrime <- function(g, nrep = NULL, strata.mat = NULL, keep.null = FALSE, ...) {
  if(is.null(strata.mat)) g <- g[, , strataNames(g)]
  if(ploidy(g) < 2 | nStrata(g) == 1) {
    return(list(
      stat.name = "G''st", 
      result = c(estimate = NA, p.val = NA),
      null.dist = NULL
    ))
  }
  
  strata.mat <- .checkStrataMat(strata.mat, g, nrep)
  
  result <- statGstDblPrime_C(
    sapply(loci(g), function(x) as.numeric(x) - 1), 
    strata.mat, ploidy(g)
  )
  
  return(.formatResult(result, "G''st", keep.null))
  
#   hets <- Hstats(g, strata)  
#   Hs <- mean(hets["Hs", ], na.rm = TRUE)
#   Ht <- mean(hets["Ht", ], na.rm = TRUE)
#   n <- if(is.null(strata)) nStrata(g) else length(unique(strata))
#   est <- (n * (Ht - Hs)) / ((n * Ht) - Hs) / (1 - Hs)
#   if(is.nan(est)) est <- NA
#   
#   names(est) <- "G''st"
#   est
}