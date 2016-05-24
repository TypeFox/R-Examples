transPrep <- function (x, trans = "ordinals", spline.intKnots = 4, spline.degree = 2,
                       missing = "none"){
  # Prepares for transformation 
  # and gives an initial form of the data
  #    trans:
  #       none     no transformation 
  #       linear   linear transformation
  #       interval   linear transformation
  #       nominal  nominal transformation
  #       ordinalp ordinal primary approach to ties (untie ties)
  #       ordinals secondary approach to ties (keep ties tied)
  #       ordinalt tertiary approach to ties (keep ties tied)
  #       spline   I-spline transformation
  #       mspline  monotone I-spline transformation
  #    missing: 
  #       none     missing values (NA) are deleted, that is, their 
  #                weights w are considered to be 0
  #       single   missing values (NA) are considered a single category that does not need to 
  #                follow the restrictions of trans
  #       multiple each missing value (NA) receives its own category and does follow the 
  #                restrictions of trans
  #            
  #
  n  <- length(x)
  knotSeq <- NULL
  base <- NULL
  #xSorted <- sort(as.vector(x), index.return = TRUE, na.last = TRUE)
  iord <- order(as.vector(x),  na.last = TRUE)  # Missing values are ordered last
  y    <- as.vector(x)[iord]
  y[is.na(y)] <- Inf       # Replace NA by Inf for the computation of tie blocks
  # Find tie blocks
  indTieBlock <-c(1,(2:n)[!y[-n]==y[-1]])
  ties <- c(indTieBlock[-1],n+1)-indTieBlock
  # Determine number of nonmissing data
  n_nonmis <- ifelse(is.infinite(y[n]), n - ties[length(ties)], n)
  iord_nonmis <- iord[1:n_nonmis]             # Order permutation for nonmissings
  if (n_nonmis < n){                          # If there are missings
    iord_mis <- iord[(n_nonmis+1):n]
    nties_nonmis <- length(ties) - 1
    if (missing=="multiple"){                 # add tieblocks of 1 for each missing value
      ties <- c(ties[-length(ties)],rep(1,n-n_nonmis))
    }
  } else {
    nties_nonmis <- length(ties)
    iord_mis <- NULL
  }
  x_unique <- x[iord[cumsum(ties[1:nties_nonmis])]]
  
  # Set xInit initial 
  if (trans %in% c("none","linear","interval")) {
    base <- cbind(rep(1,nties_nonmis),x_unique-x_unique[1])
    xInit <- rep(0,n)
    xInit[iord_nonmis] <- rep(x_unique,ties[1:nties_nonmis])
  } else if (trans %in% c("ordinalp","ordinals","ordinalt","ordinal","nominal","nominals","nominalp")){
    i <- 1:nties_nonmis
    xInit <- rep(0,n)
    xInit[iord_nonmis] <- rep(i,ties[1:nties_nonmis])
  } else if (trans %in% c("spline","mspline")){    
    if (ties[1]!=n){
      res <- splineSetUp(x_unique, spline.intKnots + 2, spline.degree)
      #base <- matrix(0, n, ncol(res$base))
      #base[iord[1:n_nonmis],] <- res$base
      base    <- res$base
      knotSeq <- res$knotSeq
    } else {
      base <- matrix(1, n_nonmis, 1)
    }
    xInit <- rowSums(base)
  }
  return(list(x = x, 
              x_unique = x_unique,
              n = n,
              n_nonmis = n_nonmis,
              trans = trans, 
              spline.allKnots = spline.intKnots + 2, 
              spline.degree = spline.degree, 
              spline.knotSeq = knotSeq,
              xInit = xInit, 
              iord = iord, 
              ties = ties, 
              nties_nonmis = nties_nonmis,
              base = base, 
              missing = missing,
              iord_nonmis = iord_nonmis,
              iord_mis = iord_mis,
              #factor = fact,
              class = "optScal"))
}