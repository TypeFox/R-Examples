atom2ele <- function(...)
  UseMethod("atom2ele")

atom2ele.default <- function(x, elety.custom=NULL, rescue=TRUE, ...){
  if(!is.null(elety.custom)) {
    if(!all(c("name","symb") %in% names(elety.custom)))
      stop("'elety.custom' must contains 'name' and 'symb' components")
    inds <- unlist(lapply(elety.custom, is.factor))
    elety.custom[inds] <- lapply(elety.custom[inds], as.character)
  }
  atom.index <- rbind(elety.custom[,c("name","symb")], atom.index[,c("name","symb")])
  # Why atom names starting by "H" are directly converted to "H" as follow?
  # x[substr(x,1,1) == "H"] <- "H"
  symb <- atom.index[match(x, atom.index[,"name"]), "symb"]
  is.unknown <- is.na(symb)
  if(any(is.unknown)) {
    if(rescue) {
      symb[is.unknown] <- substr(x[is.unknown],1,1)
      warning(paste("\n\tunknown element: mapped ", x[is.unknown], " to ", symb[is.unknown], sep=""))
    }
    else {
      stop(paste("\n\tatom2symb: element of '", x[is.unknown], "' unknown", sep=""))
    }
  }
  symb <- unlist(symb)
  return(symb)
}

atom2ele.pdb <- function(pdb, inds, elety.custom=NULL, rescue=TRUE, ...){
  if(!is.null(inds))
    pdb <- trim.pdb(pdb, inds)
  atom.names <- pdb$atom[,"elety"]
  symb <- atom2ele.default(atom.names, elety.custom, rescue, ...)
  return(symb)
}
