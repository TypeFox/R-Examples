atom2mass <- function(...)
  UseMethod("atom2mass")

atom2mass.default <- function(x, mass.custom=NULL, elety.custom=NULL,
                              grpby=NULL, rescue=TRUE, ...){
  if(!is.null(mass.custom)) {
    if(!all(c("symb","mass") %in% names(mass.custom)))
      stop("'mass.custom' must contains 'symb' and 'mass' components")
    inds <- unlist(lapply(mass.custom, is.factor))
    mass.custom[inds] <- lapply(mass.custom[inds], as.character)
  }
  elements <- rbind(mass.custom[,c("symb","mass")], elements[,c("symb","mass")])
  symb <- atom2ele.default(x, elety.custom, rescue, ...)
  M <- elements[match(symb, elements[,"symb"]), "mass"]
  
  if(any(is.na(M)))
    stop(paste("\n\tatom2mass: mass of element '", symb[is.na(M)], "' unknown", sep=""))
  
  if(!is.null(grpby)) {
    if(length(grpby) != length(M))
      warning("'grpby' as been recycled")
    M <- unlist(lapply(split(M, grpby), sum))
  }
  return(M)
}

atom2mass.pdb <- function(pdb, inds=NULL, mass.custom=NULL, elety.custom=NULL,
                          grpby=NULL, rescue=TRUE, ...){
  if(!is.null(inds))
    pdb <- trim.pdb(pdb, inds)
  
  atom.names <- pdb$atom[,"elety"]
  M <- atom2mass.default(atom.names, mass.custom, elety.custom, grpby, rescue, ...)
  return(M)
}
