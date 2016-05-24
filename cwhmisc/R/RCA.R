RCA <- function(dir, pkg, sw=0:6, echoonly=FALSE, verbose=TRUE)
{ # 2015-07-19, 17:24:10
  stopifnot(sw %in% 0:6)
  pp <- paste(dir,pkg,sep="/")
  wdold <- setwd(dir)
  oold <- options(warn=-1)
  x <- system2("grep", c( "-inr", 'version', paste(pp,"/DESCRIPTION",sep="")),,TRUE )
  version <- strsplit(x, ": ")[[1]][2]
            if (verbose)  cat("version = ", version,"\n")
  S <-paste("R CMD",c(
    " Rd2pdf --no-clean --force", ## look at mypkg.Rcheck/mypkg-manual.pdf for problems
    " build --force",
    " check",
    " check --as-cran",
    paste(" check --as-cran ",pp,"_",version,".tar.gz",sep=""),
    " install"
    ), rep( pp, 6) ) # paste
#            if (verbose) catn( "Package = ", paste(S[sw],"\n") )
  if ( (0 == sw)[1] ) for (ii in 1:6) print(paste(" RCA ( ",ii, " ): ", S[ii]))
            if (verbose) cat("old dir, current dir: ",wdold, ", ",getwd(),"\n")
  for (kk in setdiff(sw,0)) {
    if (echoonly) {  catn( kk, " ", S[kk])
    } else {
      if (verbose) print(paste(" ==== RC (", kk, ") : ",S[kk], " on: ",datetime()))
      system( S[kk] )
    } # !echoonly
 }
              if (verbose) print(datetime())
    setwd(wdold)
    options(warn=oold$warn)
  return( invisible() )
}  ##  RCA
