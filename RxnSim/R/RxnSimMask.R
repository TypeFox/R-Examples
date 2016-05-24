.mask <- function (mol, map, mask) {
  count <- .jcall(mol, 'I', 'getAtomCount')
  atoms <- c(0:(count-1))
  atoms <- atoms[!atoms %in% map$mapping[[1]]]
  
  if (!is.null(mask)) {
    .jcall(mol, 'V', 'add', mask)
    newAtom <- .jcall(mol, 'Lorg/openscience/cdk/interfaces/IAtom;', 'getLastAtom')
  }
  
  mapping <- sort(map$mapping[[1]], TRUE)
  for (atomNo in mapping) {
    atom <- .jcall(mol, 'Lorg/openscience/cdk/interfaces/IAtom;', 'getAtom', atomNo)
    connectedAtoms <- as.list(.jcall(mol, 'Ljava/util/List;', 'getConnectedAtomsList', atom))
    for (nbdAtom in connectedAtoms) {
      nbdAtom <- .jcast(nbdAtom, 'org/openscience/cdk/interfaces/IAtom')
      nbdAtomNo <- .jcall(mol, 'I', 'getAtomNumber', nbdAtom)
      bond <- .jcall(mol, 'Lorg/openscience/cdk/interfaces/IBond;', 'removeBond', atom, nbdAtom)
      if (!nbdAtomNo %in% map$mapping[[1]]) {
        if (!is.null(mask)) {
          .jcall(bond, 'V', 'setAtoms', .jarray(c(nbdAtom, newAtom), contents.class = 'org/openscience/cdk/interfaces/IAtom'))
          .jcall(mol, 'V', 'addBond', bond)
        }
      }
    }
  }
  for(atomNo in mapping) {
    .jcall(mol, 'V', 'removeAtom', atomNo)
  }
}

.meta.mask <- function(substructure, mask, mol, recursive) {
  smi <- tryCatch({
    while(1) {
      if (mask != '') {
        maskX <- .smilesParser(mask, F, F)
      } else {
        maskX <- NULL
      }
      
      tryCatch({
        map <- rcdk::matches(substructure, mol, return.matches = T)
      }, error = function(err) {
        stop('Unable to find matches.', call. = F)
        #stop(err)
      })
      if (map[[1]]$match == TRUE) {
        .mask(mol, map[[1]], maskX)
      } else {
        break
      }
      if(recursive == FALSE) {
        break
      }
    }
    
    smi <- get.smiles(mol)
  }, error = function (err) {
    stop (err)
  })
  return(smi)
}

ms.mask <- function (substructure, mask, molecule, format = 'smiles', standardize = T, explicitH = F, recursive = F) {
  if(missing(substructure) || substructure == '') {
    stop('Enter a structure to mask in form of a SMILES or SMARTS.', call. = F)
  }
  if(missing(mask)) {
    stop('Mask not specified.', call. = F)
  }
  if(missing(molecule)) {
    stop('Input molecule missing.', call. = F)
  }
  
  format <- tolower(format)
  
  
  smi <- tryCatch({    
    if (format[[1]] == 'smiles') {
      mol <- .smilesParser(molecule, standardize = standardize, explicitH = explicitH)
    } else if (format[[1]] == 'mol') {
      mol <- .molParser(molecule, standardize = standardize, explicitH = explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    mask <- sub('\\s+', '', mask)
    .meta.mask(substructure, mask, mol, recursive)
  }, error = function (err) {
      stop (err)
  })
  return(smi)
}

.rct.mask <- function (substructure, mask, rct, recursive) {
  rct <- tryCatch({
    mask <- sub('\\s+', '', mask)
    
    for (mol in rct$Reactants) {
      .meta.mask(substructure, mask, mol, recursive)
    }
    for (mol in rct$Products) {
      .meta.mask(substructure, mask, mol, recursive)
    }
    
    react <- paste(lapply(rct$Reactants, get.smiles), collapse = '.')
    react <- sub('\\.+', '.', react)
    react <- sub('^\\.', '', react)
    react <- sub('\\.$', '', react)
    
    prod <- paste(lapply(rct$Products, get.smiles), collapse = '.')
    prod <- sub('\\.+', '.', prod)
    prod <- sub('^\\.', '', prod)
    prod <- sub('\\.$', '', prod)
    
    rct$RSMI <- paste(react, prod, sep=">>")
    rct
  }, error = function (err) {
    stop (err)
  })
  return(rct)
}

rs.mask <- function (substructure, mask, reaction, format = 'rsmi', standardize = T, explicitH = F, recursive = F) {
  if(missing(substructure) || substructure == '') {
    stop('Enter a structure to mask in form of a SMILES or SMARTS.', call. = F)
  }
  if(missing(mask)) {
    stop('Mask not specified.', call. = F)
  }
  if(missing(reaction)) {
    stop('Input reaction missing.', call. = F)
  }
  
  format <- tolower(format)
  
  rsmi <- tryCatch({
    if (format[[1]] == 'rsmi') {
      rct <- .rsmiParser(reaction, standardize = standardize, explicitH = explicitH)
    } else if (format[[1]] == 'rxn') {
      rct <- .mdlParser(reaction, standardize = standardize, explicitH = explicitH)
    } else {
      stop("Invalid input format.", call. = F)
    }
    
    rct <- .rct.mask(substructure, mask, rct, recursive)
    rct$RSMI
  }, error = function (err) {
    stop (err)
  })
  return(rsmi)
}