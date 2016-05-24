#  Create an object of class 'atoms' containing the data related to ATOM and HETATM records
#  of a given molecular system stored in a PDB file.

atoms <- function(...)
  UseMethod("atoms")

atoms.default <- function(recname, eleid, elename, alt,
                                resname, chainid, resid, insert,
                                x1, x2, x3, occ, temp, segid, basis = "xyz", ...)
{
  
  recname <- as.character(recname)
  eleid   <- suppressWarnings(as.integer(eleid))
  elename <- as.character(elename)
  alt     <- as.character(alt)
  resname <- as.character(resname)
  chainid <- as.character(chainid)
  resid   <- suppressWarnings(as.integer(resid))
  insert  <- as.character(insert)
  x1      <- suppressWarnings(as.numeric(x1))
  x2      <- suppressWarnings(as.numeric(x2))
  x3      <- suppressWarnings(as.numeric(x3))
  occ     <- suppressWarnings(as.numeric(occ))
  temp    <- suppressWarnings(as.numeric(temp))
  segid   <- as.character(segid)
  
  if(any(is.na(eleid))) warning("In 'atoms': 'eleid' contains NA values")
  if(any(is.na(resid))) warning("In 'atoms': 'resid' contains NA values")
  if(any(is.na(   x1))) warning("In 'atoms':    'x1' contains NA values")
  if(any(is.na(   x2))) warning("In 'atoms':    'x2' contains NA values")
  if(any(is.na(   x3))) warning("In 'atoms':    'x3' contains NA values")
  if(any(is.na(  occ))) warning("In 'atoms':   'occ' contains NA values")
  if(any(is.na( temp))) warning("In 'atoms':  'temp' contains NA values")
  
  atoms <- data.frame(
    recname = recname,
    eleid   = eleid  ,
    elename = elename,
    alt     = alt    ,
    resname = resname,
    chainid = chainid,
    resid   = resid  ,
    insert  = insert ,
    x1     = x1      ,
    x2     = x2      ,
    x3     = x3      ,
    occ     = occ    ,
    temp    = temp   ,
    segid   = segid  ,
    stringsAsFactors = FALSE
  )
  attr(atoms, "basis") <- basis
  
  class(atoms) <- c("atoms","coords","data.frame")
  return(atoms)
}

is.atoms <- function(x)
{
  to.return <- any(attr(x,which="class") == "atoms")
  return(to.return)
}