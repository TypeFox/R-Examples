
.dFvGet <- function()
{
   f.res <- .Fortran("zdfvals",io=to.integer(0),dfv=single(66))
   z <- f.res$dfv
   names(z) <- c(
     "tlo", "gma", "mxs", "mxt", "ntm", "tua", "tlu", "iop", "ix1", "iy1",
     "ic1", "ini", "isr", "itc", "icn", "alf", "ccc", "upr", "tli", "isq",
     "isg", "ite", "itw", "mxf", "mxn", "mxg", "iwg", "apr", "icv", "xfd",
     "ia1", "ia2", "fff", "ff1", "fu1", "fb1", "ial", "mxe", "fct", "ffc",
     "ial", "mxe", "tls", "tlr", "msx", "ik1", "ipt", "isd", "ich", "esp", 
     "ilc", "aa2", "bb2", "em" , "cr" , "enu", "tlv", "tlm", "ith", "ilm",
     "ddd", "ics", "mxx", "ilg", "ipo", "iug")
  as.list(z)
}

.dFvSet <- function(def)
{
  f.res <- .Fortran("zdfvals",io=to.integer(1),dfv=to.single(def))
  return()
}

.dFvPut <- function(vals,nams)
{
  alldef <-  c(
     "tlo", "gma", "mxs", "mxt", "ntm", "tua", "tlu", "iop", "ix1", "iy1",
     "ic1", "ini", "isr", "itc", "icn", "alf", "ccc", "upr", "tli", "isq",
     "isg", "ite", "itw", "mxf", "mxn", "mxg", "iwg", "apr", "icv", "xfd",
     "ia1", "ia2", "fff", "ff1", "fu1", "fb1", "ial", "mxe", "fct", "ffc",
     "ial", "mxe", "tls", "tlr", "msx", "ik1", "ipt", "isd", "ich", "esp", 
     "ilc", "aa2", "bb2", "em" , "cr" , "enu", "tlv", "tlm", "ith", "ilm",
     "ddd", "ics", "mxx", "ilg", "ipo", "iug")
  ll  <- length(nams)
  def <- .Fortran("zdfvals",io=to.integer(0),dfv=single(66))$dfv
  for (i in 1:ll) {
     zi <- match(nams[i],alldef,nomatch=0)
     if (zi==0) {cat("Unknown parameter",nams[i],"! \n"); next}
     def[zi] <- vals[i]
  }   
  f.res <- .Fortran("zdfvals",io=to.integer(1),dfv=to.single(def))
  return()
}

"dfvals" <- function() {
  f.res <- .Fortran("zdfvals",io=to.integer(-1),dfv=single(66))
  return()
}


