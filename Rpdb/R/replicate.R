#  Duplicate atomic coordinates using periodical boundary conditions

replicate <- function(x, ...)
  UseMethod("replicate")

replicate.coords <- function(x, cryst1 = NULL, a.ind = 0, b.ind = 0, c.ind = 0, ...)
{
  if(!is.coords(x)) stop("'x' must be an object of class 'coords'")
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
  b <- basis(x)
  if(b == "xyz")
  {
    if( is.null(cryst1))   stop("Please specify a 'cryst1' object")
    if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
    x <- xyz2abc(x, cryst1) 
  }
  
  abc.ind <- expand.grid(a.ind, b.ind, c.ind)
  
  L <- apply(abc.ind, 1,
             function(abc, x)
             {
               x$x1 <- x$x1 + abc[1]
               x$x2 <- x$x2 + abc[2]
               x$x3 <- x$x3 + abc[3]
               return(x)
             }, x)
  
  x <- do.call(rbind, L)
  if(b == "xyz") x <- abc2xyz(x, cryst1)
  
  return(x)
}

replicate.atoms <- function(x, cryst1 = NULL, a.ind = 0, b.ind = 0, c.ind = 0, ...)
{
  if(!is.atoms(x)) stop("'x' must be an object of class 'atoms'")
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
  basis <- basis(x)
  if(basis == "xyz")
  {
    if( is.null(cryst1))   stop("Please specify a 'cryst1' object")
    if(!is.cryst1(cryst1)) stop("'cryst1' must be an object of class 'cryst1'")
    x <- xyz2abc(x, cryst1) 
  }
  
  abc.ind <- expand.grid(a.ind, b.ind, c.ind)
  
  L <- apply(abc.ind, 1,
             function(abc, x)
             {
               x$x1 <- x$x1 + abc[1]
               x$x2 <- x$x2 + abc[2]
               x$x3 <- x$x3 + abc[3]
               return(x)
             }, x)
  
  eleid <- rep(x$eleid, nrow(abc.ind)) + rep(1:nrow(abc.ind)-1,each=natom(x))*max(x$eleid)
  resid <- rep(x$resid, nrow(abc.ind)) + rep(1:nrow(abc.ind)-1,each=natom(x))*max(x$resid)
  
  x <- do.call(rbind, L)
  x$resid <- resid
  x$eleid <- eleid
  
  if(basis == "xyz") x <- abc2xyz(x, cryst1)
  
  return(x)
}

replicate.pdb <- function(x, a.ind = 0, b.ind = 0, c.ind = 0, cryst1 = NULL, ...)
{
  if(!is.pdb(x)) stop("'x' must be an object of class 'pdb'")
  
  if(is.null(cryst1))
    cryst1 <- x$cryst1
  
  a.ind <- unique(a.ind)
  b.ind <- unique(b.ind)
  c.ind <- unique(c.ind)
  
  na <- length(a.ind)
  nb <- length(b.ind)
  nc <- length(c.ind)
  
  ncell <- na*nb*nc
  
  eleid.1 <- x$conect$eleid.1
  eleid.2 <- x$conect$eleid.2
  eleid.1 <- rep(eleid.1, ncell) + rep(1:ncell-1,each=length(eleid.1))*max(x$atoms$eleid)
  eleid.2 <- rep(eleid.2, ncell) + rep(1:ncell-1,each=length(eleid.2))*max(x$atoms$eleid)
  conect <- conect.default(eleid.1, eleid.2)
  
  atoms <- replicate.atoms(x$atoms, cryst1, a.ind, b.ind, c.ind)
  
  abc <- x$cryst1$abc*c(diff(range(a.ind))+1, diff(range(b.ind))+1, diff(range(c.ind))+1)
  cryst1 <- cryst1(abc = abc, abg = x$cryst1$abg, sgroup = x$cryst1$sgroup)

  x <- pdb(atoms, cryst1, conect)
  
  return(x)
}