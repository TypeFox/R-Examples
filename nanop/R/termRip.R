termRip <- function(pdf, Qmax, dr, maxR, maxRTermRip)
{
  lenRip <- min(maxRTermRip/dr, maxR/dr)
  rmax <- max(dr*length(pdf), maxR)
  .C("termRip",
     res = as.double(rep(0,length(pdf))),
     pdf = as.double(pdf), 
     len = as.integer(length(pdf)),
     qmax = as.double(Qmax),
     deltar = as.double(dr),
     rmax = as.double(maxR),
	 lenRip = as.integer(lenRip),
     PACKAGE="nanop")$res

}