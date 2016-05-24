trim <- function(...)
  UseMethod("trim")

"trim.pdb" <- function(pdb, ..., inds=NULL, sse=TRUE) {

  if(!is.pdb(pdb))
    stop("input 'pdb' must be a PDB list object as returned from 'read.pdb'")

  cl <- match.call()
  extra.args <- list(...)

  if(length(extra.args)>0) {
    if(!is.null(inds))
      warning("Multiple atom selection terms provided. Using only argument 'inds'")
    else if(is.select(extra.args[[1]]))
      # to be back-compatible with the habit calling trim.pdb(pdb, inds)
      inds = extra.args[[1]]
    else
      inds = atom.select(pdb, ...)
  }

  if(is.null(inds))
    stop("no selection indices provided")

  if(!is.list(inds))
    stop("selection indices must be provided i.e. from 'atom.select'")

  if(is.null(inds$atom) || is.null(inds$xyz))
    stop("selection indices must be provided i.e. from 'atom.select'")

  ## Trim atom components
  atom <- pdb$atom[inds$atom,]

  ## Add calpha indices if non-existing
  if(is.null(pdb$calpha)) {
    ca.inds <- atom.select(pdb, "calpha")
    pdb$calpha <- rep(FALSE, nrow(pdb$atom))
    pdb$calpha[ca.inds$atom] <- TRUE
  }

  ## Trim calpha indices
  calpha <- pdb$calpha[inds$atom]

  ## Trim xyz components
  xyz <- trim.xyz(pdb$xyz, col.inds = inds$xyz)

  ## Trim SSE components
  helix <- NULL; sheet <- NULL;
  if(sse) {
    ss <- pdb2sse(pdb, verbose = FALSE)

    ##- Trim sse vector
    calpha2 <- which(pdb$calpha) %in% inds$atom
    ss <- ss[calpha2]

    ##- New sse
    new.sse <- bounds.sse(ss)

    helix <- new.sse$helix
    if(length(helix$start) > 0) {
       ##- add back other components
       add <- pdb$helix[!names(pdb$helix) %in% names(new.sse$helix)]
       ##- match sse number in case some sse are completely removed
       add <- lapply(add, function(x) x[new.sse$helix$id])
       helix <- c(helix, add)
    }

    sheet <- new.sse$sheet
    if(length(sheet$start) > 0) {
       ##- add back other components
       add <- pdb$sheet[!names(pdb$sheet) %in% names(new.sse$sheet)]
       ##- match sse number in case some sse are completely removed
       add <- lapply(add, function(x) x[new.sse$sheet$id])
       sheet <- c(sheet, add)
    }

    ##- remove 'id'; Maybe we don't need it?
    helix$id <- NULL
    sheet$id <- NULL
  }

  output <- list(atom   = atom,
                 helix  = helix,
                 sheet  = sheet,
                 seqres = pdb$seqres, ## return unmodified
                 xyz    = xyz,
                 calpha = calpha, 
                 call   = cl)

  class(output) <- class(pdb)
  return(output)
}
