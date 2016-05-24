calcPDFR <- function(nanop, calpha=1, dr=.01, rmin=1, rmax=20, p = 1,
                    foranalcs = FALSE) {

  r <- seq(rmin, rmax, by = dr)
  if(length(nanop) == 0)
    return(list(r=r,gr=rep(0,length(r))))
  if(!foranalcs) 
    ret <- list(r=r, gr=.C("calcPDFR",
                       res = as.double(rep(0,length(r))),
                       r = as.double(r), 
                       len = as.integer(length(r)),
                       np = as.double(as.vector(t(nanop))),
                       nrow = as.integer(nrow(nanop)),
                       calpha = as.double(calpha),
                       dr = as.double(dr),
                       rmin = as.double(rmin),
                       p = as.double(p),
                       PACKAGE="nanop")$res)
  else {
    rc <- attributes(nanop)$rowcore
    rsh <- attributes(nanop)$rowshell
    if(rc > 0)
      pdfC <- .C("calcPDFR",
               res = as.double(rep(0,length(r))),
               r = as.double(r), 
               len = as.integer(length(r)),
               np = as.double(as.vector(t(nanop[1:rc,]))),
               nrow = as.integer(rc),
               calpha = as.double(calpha),
               dr = as.double(dr),
               rmin = as.double(rmin),
               p = as.double(p),
               PACKAGE="nanop")$res
    else
      pdfC <- rep(0,length(r))
    if(rsh > 0)
    pdfS <- .C("calcPDFR",
               res = as.double(rep(0,length(r))),
               r = as.double(r), 
               len = as.integer(length(r)),
               np = as.double(as.vector(t(nanop[(rc+1):nrow(nanop),]))),
               nrow = as.integer(rsh),
               calpha = as.double(calpha),
               dr = as.double(dr),
               rmin = as.double(rmin),
               p = as.double(p),
               PACKAGE="nanop")$res
    else
      pdfS <- rep(0,length(r))
    if(rc==0||rsh==0)
      pdfCS <- rep(0,length(r))
    else
      pdfCS <- .C("calcPDF_CS",
                res = as.double(rep(0,length(r))),
                r = as.double(r), 
                len = as.integer(length(r)),
                npC = as.double(as.vector(t(nanop[1:rc,]))),
                nrowC = as.integer(rc),
                npS = as.double(as.vector(t(nanop[(rc+1):nrow(nanop),]))),
                nrowS = as.integer(rsh),
                calpha = as.double(calpha),
                dr = as.double(dr),
                rmin = as.double(rmin),
                p = as.double(p),
                PACKAGE="nanop")$res
    
    
    ret <- list(r=r,
                gr=(rc/nrow(nanop) * pdfC) + (rsh/nrow(nanop) * pdfS) + pdfCS,
                grC=rc/nrow(nanop) * pdfC, grS=rsh/nrow(nanop) * pdfS,
                grCS=pdfCS, rowcore=rc, rowshell=rsh)
    
  }
  
  ret
}
