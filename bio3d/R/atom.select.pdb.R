".is.protein" <- function(pdb, byres=TRUE) {
  if(byres) {
    return(.is.protein1(pdb))
  }
  else {
    ## possible option to issue a warning when the two methods diverge
    sel1 <- .is.protein1(pdb)
    sel2 <- .is.protein2(pdb)
    
    if(!(identical(sel1, sel2))) {
      sel <- cbind(sel1, sel2)
      sums <- apply(sel, 1, sum)
      inds <- which(sums==1)
      unq <- paste(unique(pdb$atom$resid[inds]), collapse=",")
      warning(paste("possible protein residue(s) with non-standard residue name(s) \n   (", unq, ")"))
    }
    return(sel1)
  }
}

".is.protein1" <- function(pdb) {
  aa <- bio3d::aa.table$aa3
  return(pdb$atom$resid %in% aa)
}

".is.protein2" <- function(pdb) {
  resid <- paste(pdb$atom$chain, pdb$atom$insert, pdb$atom$resno, sep="-")

  at.ca <- resid[ pdb$atom$elety == "CA"]
  at.o  <- resid[ pdb$atom$elety == "O" ]
  at.c  <- resid[ pdb$atom$elety == "C" ]
  at.n  <- resid[ pdb$atom$elety == "N" ]

  common <- intersect(intersect(intersect(at.ca, at.o), at.n), at.c)
  return(resid %in% common)
}

".is.nucleic" <- function(pdb) {
  nuc.aa <- c("A",   "U",  "G",  "C",   "T",  "I",
              "DA", "DU", "DG", "DC",  "DT", "DI")
  return(pdb$atom$resid %in% nuc.aa)
}

".is.water" <- function(pdb) {
  hoh <- c("H2O", "OH2", "HOH", "HHO", "OHH", "SOL",
           "WAT", "TIP", "TIP2", "TIP3", "TIP4")
  return(pdb$atom$resid %in% hoh)
}

".is.hydrogen" <- function(pdb) {
  return(substr( gsub("^[123]", "",pdb$atom$elety) , 1, 1) %in% "H")
}

.match.type <- function(pdb, t) {
  if(!is.character(t))
    stop("'type' must be a character vector")
  pdb$atom$type %in% t
}

.match.eleno <- function(pdb, eleno) {
  if(!is.numeric(eleno))
    stop("'eleno' must be a numeric vector")
  pdb$atom$eleno %in% eleno
}

.match.elety <- function(pdb, elety) {
  if(!is.character(elety))
    stop("'elety' must be a character vector")
  pdb$atom$elety %in% elety
}

.match.resid <- function(pdb, resid) {
  if(!is.character(resid))
    stop("'resid' must be a character vector")
  pdb$atom$resid %in% resid
}

.match.chain <- function(pdb, chain) {
  if(!is.character(chain))
    stop("'chain' must be a character vector")
  pdb$atom$chain %in% chain
}

.match.resno <- function(pdb, resno) {
  if(!is.numeric(resno))
    stop("'resno' must be a numeric vector")
  pdb$atom$resno %in% resno
}

.match.segid <- function(pdb, segid) {
  if(!is.character(segid))
    stop("'segid' must be a character vector")
  pdb$atom$segid %in% segid
}

atom.select.pdb <- function(pdb, string = NULL,
                            type  = NULL, eleno = NULL, elety = NULL,
                            resid = NULL, chain = NULL, resno = NULL,
                            segid = NULL, operator = "AND", inverse = FALSE,
                            value = FALSE, verbose=FALSE,  ...) {

  if(!is.pdb(pdb))
    stop("'pdb' must be an object of class 'pdb'")

  ## check input operator
  op.tbl <- c(rep("AND",3), rep("OR",4))
  operator <- op.tbl[match(operator, c("AND","and","&","OR","or","|","+"))]
  if(!operator %in% c("AND", "OR"))
    stop("Allowed values for 'operator' are 'AND' or 'OR'")

  ## check input string
  if(!is.null(string)) {
    str.allowed <- c("all", "protein", "notprotein", "nucleic", "notnucleic", "water", "notwater",
                     "calpha", "cbeta", "backbone", "back", "ligand", "h", "noh")
    if(!(string %in% str.allowed))
      stop("Unknown 'string' keyword. See documentation for allowed values")
  }

  ## verbose message output
  if(verbose) cat("\n")
  .verboseout <- function(M, type) {
    cat(" .. ", sprintf("%08s", length(which(M))), " atom(s) from '", type, "' selection \n", sep="")
  }

  ## combine logical vectors
  .combinelv <- function(L, M, operator) {
    if(operator=="AND") M <- L & M
    if(operator=="OR") M <- L | M
    return(M)
  }
  
  cl <- match.call()
  if(operator=="AND")
    M <- rep(TRUE, nrow(pdb$atom))
  if(operator=="OR")
    M <- rep(FALSE, nrow(pdb$atom))
  
  if(!is.null(string)) {   
    M <- switch(string,
                all         =   M <- rep(TRUE, nrow(pdb$atom)),
                protein     =  .is.protein(pdb),
                notprotein  = !.is.protein(pdb),
                nucleic     =  .is.nucleic(pdb),
                notnucleic  = !.is.nucleic(pdb),
                water       =  .is.water(pdb),
                notwater    = !.is.water(pdb),
                calpha      =  .is.protein(pdb)  & .match.elety(pdb, "CA"),
                cbeta       =  .is.protein(pdb)  & .match.elety(pdb, c("CA", "N", "C", "O", "CB")),
                backbone    =  .is.protein(pdb)  & .match.elety(pdb, c("CA", "N", "C", "O")),
                back        =  .is.protein(pdb)  & .match.elety(pdb, c("CA", "N", "C", "O")),
                ligand      = !.is.protein(pdb)  & !.is.nucleic(pdb) & !.is.water(pdb),
                h           =  .is.hydrogen(pdb),
                noh         = !.is.hydrogen(pdb),
                NA
                )
    
    if(verbose) {
      .verboseout(M, 'string')
    }
  }
  
  if(!is.null(type)) {
    L <- .match.type(pdb, type)
    if(verbose) .verboseout(L, 'type')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(eleno)) {
    L <- .match.eleno(pdb, eleno)
    if(verbose) .verboseout(L, 'eleno')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(elety)) {
    L <- .match.elety(pdb, elety)
    if(verbose) .verboseout(L, 'elety')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(resid)) {
    L <- .match.resid(pdb, resid)
    if(verbose) .verboseout(L, 'resid')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(chain)) {
    L <- .match.chain(pdb, chain)
    if(verbose) .verboseout(L, 'chain')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(resno)) {
    L <- .match.resno(pdb, resno)
    if(verbose) .verboseout(L, 'resno')
    M <- .combinelv(L, M, operator)
  }
  if(!is.null(segid)) {
    L <- .match.segid(pdb, segid)
    if(verbose) .verboseout(L, 'segid')
    M <- .combinelv(L, M, operator)
  }

  if(verbose)
    cat(" ..", sprintf("%08s", length(which(M))), "atom(s) in final combined selection \n")

  if(inverse) {
    if(verbose) {
      cat(" ..", sprintf("%08s", length(which(!M))), "atom(s) in inversed selection \n")
    }
    sele <- as.select(which(!M))
  }
  else
    sele <- as.select(which(M))

  sele$call <- cl
  if(verbose) cat("\n")

  if(value)
    return(trim.pdb(pdb, sele))
  else
    return(sele)
}
