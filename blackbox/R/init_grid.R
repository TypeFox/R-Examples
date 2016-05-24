init_grid <- function(lower=c(par=0),upper=c(par=1),steps=NULL,nUnique=NULL,nRepl=min(10L,nUnique),jitterFac=0.5) {
  np <- length(lower)
  parnames <- names(lower)
  if ( length(parnames)!=np) {
    stop("Some names missing in 'lower'. Check input.")
  }
  if ( length(names(upper))!=np) {
    stop("names(lower) differ from names(upper). Check input.")
  }
  if ( ! all(parnames==names(upper))) {
    stop("names(lower) differ from names(upper). Check input.")
  }
  whichvar <- which(lower != upper)
  nvarp <- length(whichvar)
  if (is.null(nUnique)) nUnique <- floor(50^((nvarp/3)^(1/3)))
  if (is.null(steps)) steps <- ceiling(nUnique^(1/nvarp)) ## -> 2 in large dim
  if (length(steps)==1L) {
    steps <- rep(steps,nvarp)
    names(steps) <- names(whichvar)
  }
  if (length(steps)!=nvarp) stop("'steps' has incorrect length. Check input.")
  dx <- (upper[whichvar]-lower[whichvar])/(steps-1+2*jitterFac)
  arglist <- list()
  for (it in seq(np)) {
    if ( is.character(lower[it]) || lower[it]==upper[it] ) {
      arglist[[parnames[it]]] <- lower[[it]]
    } else {
      dxl <- dx[parnames[it]]
      arglist[[parnames[it]]] <- seq(lower[it]+dxl*jitterFac,upper[it]-dxl*jitterFac, dxl)
    }
  }
  ## regular grid
  grille <- do.call(expand.grid,arglist)
  ## add noise
  ng <- nrow(grille)
  for (vit in whichvar) grille[,vit] <- grille[,vit]+2*(runif(ng)-0.5)*dx[parnames[vit]]*jitterFac
  ## reduce to nbUnique
  rownames(grille) <- as.character(seq(nrow(grille)))
  if (ng>nUnique) {
    rownams <- greedyMAXMINwithFixed(as.matrix(grille[sample(nrow(grille)),whichvar,drop=FALSE]), nUnique, dx^2, fixedNbr=1) ## facon tordue de randomiser les fixed points...
    grille <- grille[rownams,,drop=FALSE]
  }
  ## add replicates
  ng <- nrow(grille)
  grille <- rbind(grille,grille[sample(ng,nRepl),,drop=FALSE])
  attr(grille,"LOWER") <- lower
  attr(grille,"UPPER") <- upper
  return(grille)
}
