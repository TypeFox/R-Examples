#' Calculate variances of paramater contrasts
#' 
#' Calculate variances of paramater contrasts
#' 
#' See the vignette of this package for further details.
#' 
#' @param design Cross-over design.
#' @param model Model - one of the following numbers or Strings: 1 = "Standard
#' additive model", 2 = "Self-adjacency model", 3 = "Proportionality model", 4
#' = "Placebo model", 5 = "No carry-over into self model", 6 = "Treatment decay
#' model", 7 = "Full set of interactions", 8 = "Second-order carry-over
#' effects"
#' @param placebos The number of placebo treatments in the placebo model.
#' @param ppp The proportionality parameter for the proportionality model.
#' @param contrasts Optionally a contrast matrix or a list of contrast matrix.
#' If missing pairwise differences for treatment and carry-over parameters
#' are calculated.
#' @param v Number of treatments
#' @return A list with the variances of the pairwise differences or specified contrasts.
#' If contrasts are not estimable, NA is returned for variances.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @references Jones, B., & Kenward, M. G. (2003). Design and analysis of
#' cross-over trials (Vol. 98). Chapman & Hall.
#' @keywords misc
#' @examples
#' 
#' general.carryover(getDesign("fletcher1"), model=1)
#' general.carryover(getDesign("fletcher1"), model=2)
#' general.carryover(getDesign("fletcher1"), model=3)
#' general.carryover(getDesign("switchback4t"), model=7)
#' 
#' @export general.carryover
general.carryover <-function(design, v=length(table(design)), model, ppp=0.5, placebos=1, contrasts) {    
  if (class(design)=="CrossoverSearchResult") {
    if(missing(model)) {
      model <- design@model
    } else {
      if (model!=design@model) warning("Model from object does not equal specified model")
    }
    design <- design@design      
  }
  if (class(design)=="CrossoverDesign") {
    if(missing(model)) {
      model <- design@model
    } else {
      if (model!=design@model) warning("Model from object does not equal specified model")
    }
    design <- design@design
  }
  model <- getModelNr(model)    
  if(missing(contrasts)) {
    contrasts <- getPairwiseContrasts(model, v)
  }
  if (is.matrix(contrasts)) {
    contrasts <- list(contrasts)
  }
  
  H <- linkMatrix(model, v, ppp, placebos)  
  p <- dim(design)[1]
  s <- dim(design)[2]
  # Row-Column-Design:
  rcDesign <- rcd(design, v=v, model=model)
  Xr <- rcdMatrix(rcDesign, v=v, model=model)
  #X <- Xr %*% H
  Ar <- t(Xr) %*% (diag(s*p)-getPZ(s,p)) %*% Xr
  A <- t(H) %*% Ar %*% H
  
  result <- list()
  for (C in contrasts) {
    n <- sum(colSums(C==0)<dim(C)[1])
    m <- matrix(0, n, n)
    if (!estimable(design, v, model, C)) {
      # TODO We could check which rows of the contrast matrix are not estimable.
      m[row(m)!=col(m)] <- NA
    } else {
      m[lower.tri(m)] <- diag(C %*% ginv(A) %*% t(C))
      m[upper.tri(m)] <- t(m)[upper.tri(m)]
    }
    result <- c(result, list(m))
  }
  names(result)[1] <- "Var.trt.pair"
  if (length(result)==2) {
    names(result)[2] <- "Var.car.pair"
  }
  if (length(result)==3) {
    names(result)[2:3] <- c("Var.car.pair.1", "Var.car.pair.2")
  }
  result$model <- model
  return(result)
}
