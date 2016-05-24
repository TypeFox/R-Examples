installPackages <- function(x,
                            minCount = sqrt(x[1, 'Count']),
                            ...){
##
## 1.  pkgs needed
##
  pkgsum <- attr(x, 'PackageSummary')
  xName <- substring(deparse(substitute(x)), 1, 25)
  if (is.null(pkgsum))
    stop('not a findFn object;  does not have ',
         'attribute PackageSummary:  ', xName)
  if (!is.data.frame(x))
    stop('not a findFn object:  is not a data.frame:  ', xName)
  if (!all(c('Count', 'Package') %in% names(x)))
    stop('Must have columns Count & Package:  ', xName)
  sel <- (pkgsum$Count >= minCount)
  toget <- pkgsum$Package[sel]
##
## 2.  Installed pkgs?
##
  instPkgs <- .packages(TRUE)
  notInst <- toget[!(toget %in% instPkgs)]
##
## 3.  get not installed
##
  if (length(notInst) > 0)
    install.packages(notInst)
}

