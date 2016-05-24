"ctsub"<-
  function(x, y, z)
{
  junk <- .Fortran("ctsub",
		   length(x),
		   as.double(x),
		   as.double(y),
		   as.double(z),
		   ans=double(length(x)),
                   PACKAGE = "bootstrap")
  return(junk$ans)
}
