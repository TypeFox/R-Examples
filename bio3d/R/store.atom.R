"store.atom" <-
function(pdb) {

  colpaste <- function(x, col.names = colnames(x)) {
    apply(x, 1, function(row) paste(row[col.names], collapse = "."))
  }
  getinds <- function(atoms, ref = atom.names) {
    sort(atom2xyz(charmatch(atoms, ref)))
  }
  repadd <- function(num, nrep = nres, toadd = nxyz) {
    c(num, rep(num, (nrep - 1)) + rep(cumsum(rep(toadd, (nrep - 1))), each = length(num)))
  }
  atom.data <- colpaste(pdb$atom, c("elety", "resno", "chain"))
  atom.list <- matrix(unlist(strsplit(atom.data, "\\.")), ncol = 3, byrow = TRUE)
  res.data <- colpaste(pdb$atom, c("resno", "chain"))
  res.list <- unique(res.data)
  atom.names <- c("N", "CA", "C", "O", "CB", "*G", "*G1", "*G2",
                  "*D", "*D1", "*D2", "*E", "*E1", "*E2", "*Z", "NH1",
                  "NH2", "OH") ##
  atom.greek <- c("N", "CA", "C", "O", "CB", "G", "G1", "G2",
                  "D", "D1", "D2", "E", "E1", "E2", "Z", "*", "*", "*")

  coords <- NULL
# Changed for PDB format v3.3  
#  blank <- matrix(NA, nrow = 13, ncol = length(atom.names))
  blank <- matrix(NA, nrow = ncol(pdb$atom), ncol = length(atom.names))
  for (i in 1:length(res.list)) {
    res.blank <- blank
    res.ind <- which(res.list[i] == res.data)
    
    blank.ind <- charmatch(atom.list[res.ind, 1], atom.names,
                 nomatch = 0) + charmatch(substr(atom.list[res.ind,1], 2, 4),
                   atom.greek, nomatch = 0)

    res.blank[, blank.ind[blank.ind != 0]] <-
      t( pdb$atom[(res.ind[blank.ind != 0]),] )

    coords <- cbind(coords, res.blank)
  }
  natm <- length(atom.names)
  # PDB format v3.3
  nxyz <- ncol(pdb$atom) * natm
  #  nxyz <- 13 * natm
  nres <- length(coords)/(nxyz)
  dim(coords) <- c(ncol(pdb$atom), natm, nres)
#  dim(coords) <- c(13, natm, nres)
  
  dimnames(coords) = list(atom = colnames(pdb$atom), type = atom.names,
            res = res.list)
  return(coords)
}

