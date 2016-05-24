as.pdb <- function(...)
  UseMethod("as.pdb")

as.pdb.default <- function(pdb=NULL, xyz=NULL, type = NULL, resno = NULL,
                           resid = NULL, eleno = NULL,
                           elety = NULL, chain = NULL,
                           insert= NULL, alt = NULL,
                           o = NULL,  b = NULL,
                           segid = NULL, elesy = NULL,
                           charge = NULL, verbose=TRUE, ...) {
  cl <- match.call()

  ## which input argument to determine number of atoms from
  input <- list(pdb=pdb, xyz=xyz, eleno=eleno, resno=resno, resid=resid)
  nulls <- unlist(lapply(input, is.null))
  inds  <- which(!nulls)

  if(length(inds)==0)
    stop("insufficient arguments. provide 'pdb', 'xyz', 'eleno', 'resno', and/or 'resid'")

  ## check content of pdb
  if(!is.null(pdb)) {
    if(!is.pdb(pdb))
      stop("'pdb' must be of class 'pdb' as obtained from 'read.pdb'")
  }

  ## check content of xyz
  if(!is.null(xyz)) {
    if(!(is.numeric(xyz) & (is.matrix(xyz) | is.vector(xyz))))
      stop("'xyz' must be a numeric vector/matrix")
    
    xyz <- as.xyz(xyz)
  }

  ## if pdb is provided use it to determine natoms
  if (inds[1]==1) {
    natoms <- nrow(pdb$atom)
  }
  ## if xyz is provided use it to determine natoms
  else if (inds[1]==2) {
    natoms <- ncol(xyz)/3
  }
  ## else use eleno, resno, or resid
  else {
    natoms <- length(input[[inds[1]]])
  }

  if(verbose) {
    cat("\n")
    cat(" Summary of PDB generation:\n")
    cat(paste(" .. number of atoms in PDB determined by '", names(input)[inds[1]], "'\n", sep=""))
  }

  ## set value of 'xyz'
  if(!is.null(xyz)) {
    if((ncol(xyz)/3)!=natoms)
      stop("ncol(xyz)/3 != length(resno)")
  }
  else {
    if(!is.null(pdb))
      xyz <- as.xyz(pdb$xyz)
    else
      xyz <- as.xyz(rep(NA, natoms*3))
  }

  ## generic function to set the values of remaining columns of PDB
  .setval <- function(values=NULL, typ=NULL, default=NULL, class="character", repfirst=FALSE) {
    if(!is.null(values)) {
      if(class=="character") fun=is.character
      if(class=="numeric") fun=is.numeric
      if(!fun(values))
        stop(paste("'", typ, "' must be a ", class, " vector", sep=""))

      if(length(values)==1 & repfirst)
        values <- rep(values, natoms)

      if(length(values)!=natoms)
        stop(paste("length(", typ, ") != natoms", sep=""))
    }
    else {
      if(!is.null(pdb))
        values <- pdb$atom[[typ]]
      else
        values <- default
    }
    return(values)
  }

  type   <- .setval(type,   typ="type",   default=rep("ATOM", natoms), class="character", repfirst=TRUE)
  eleno  <- .setval(eleno,  typ="eleno",  default=seq(1, natoms),      class="numeric",   repfirst=FALSE)
  elety  <- .setval(elety,  typ="elety",  default=rep("CA", natoms),   class="character", repfirst=TRUE)
  resno  <- .setval(resno,  typ="resno",  default=seq(1, natoms),      class="numeric",   repfirst=FALSE)
  chain  <- .setval(chain,  typ="chain",  default=rep("A", natoms),    class="character", repfirst=TRUE)
  resid  <- .setval(resid,  typ="resid",  default=rep("ALA", natoms),  class="character", repfirst=TRUE)
  elesy  <- .setval(elesy,  typ="elesy",  default=rep(NA, natoms),     class="character", repfirst=TRUE)
  segid  <- .setval(segid,  typ="segid",  default=rep(NA, natoms),     class="character", repfirst=TRUE)
  o      <- .setval(o,      typ="o",      default=rep(NA, natoms),     class="numeric",   repfirst=TRUE)
  b      <- .setval(b,      typ="b",      default=rep(NA, natoms),     class="numeric",   repfirst=TRUE)
  alt    <- .setval(alt,    typ="alt",    default=rep(NA, natoms),     class="character", repfirst=FALSE)
  insert <- .setval(insert, typ="insert", default=rep(NA, natoms),     class="character", repfirst=FALSE)
  charge <- .setval(charge, typ="charge", default=rep(NA, natoms),     class="numeric",   repfirst=TRUE)

  ## make the data frame for the final PDB object
  atom        <- list()
  atom$type   <- type
  atom$eleno  <- eleno
  atom$elety  <- elety
  atom$alt    <- alt
  atom$resid  <- resid
  atom$chain  <- chain
  atom$resno  <- resno
  atom$insert <- insert
  atom$x      <- xyz[1, seq(1, natoms*3, by=3)]
  atom$y      <- xyz[1, seq(2, natoms*3, by=3)]
  atom$z      <- xyz[1, seq(3, natoms*3, by=3)]
  atom$o      <- o
  atom$b      <- b
  atom$segid  <- segid
  atom$elesy  <- elesy
  atom$charge <- charge
  atom        <- data.frame(atom, stringsAsFactors=FALSE)

  out        <- list()
  out$atom   <- atom
  out$xyz    <- xyz
  class(out) <- "pdb"

  ## should account for new resno and chain
  #if(!is.null(pdb)) {
    #out$helix  <- pdb$helix
    #out$sheet  <- pdb$sheet
    #out$seqres <- pdb$seqres
    #class(out) <- class(pdb)
  #}

  unwhich <- function(x, n) {
    out <- rep_len(FALSE, n)
    out[x] <- TRUE
    return(out)
  }

  ca.inds    <- atom.select(out, "calpha", verbose=verbose)
  out$calpha <- unwhich(ca.inds$atom, natoms)
  out$call   <- cl

  if(verbose) {
    resid <- unique(paste(atom$chain, atom$resno, sep="-"))
    cat(paste(" .. number of atoms in PDB: ", natoms, "\n"))
    cat(paste(" .. number of calphas in PDB:", sum(out$calpha), "\n"))
    cat(paste(" .. number of residues in PDB:", length(resid), "\n"))
    cat("\n")
  }

  return(out)
}
