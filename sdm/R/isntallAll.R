# Author: Babak Naimi, naimi.b@gmail.com
# Date :  April 2016
# Version 1.3
# Licence GPL v3
#--------


.getPackageList <- function() {
  methodInfo <- NULL
  pkgs <- c()
  lst <- list.files(system.file("methods/sdm", package="sdm"),pattern='R$',full.names = TRUE)
  for (l in lst) {
    source(l,local=TRUE)
    p <- methodInfo$packages
    p <- p[!p == '.tmp']
    pkgs <- c(pkgs,p)
  }
  p <- c('shiny','rgdal','raster')
  unique(c(pkgs,p))
}


if (!isGeneric("installAll")) {
  setGeneric("installAll", function(pkgs,update,...)
    standardGeneric("installAll"))
}


setMethod('installAll', signature(pkgs='ANY'),
          function(pkgs,update=FALSE,...) {
            if (missing(update)) update <- FALSE
            pl <- .getPackageList()
            if (!update) {
              p <- pl[!.is.installed(pl)]
              if (length(p) > 0) {
                s <- rep(TRUE,length(p))
                for (i in seq_along(p)) {
                  pi <- try(install.packages(p[i],...),silent = TRUE)
                  if (inherits(pi, "try-error")) s[i] <- FALSE
                }
                if (any(!s)) {
                  if (any(s)) {
                    cat(paste('\n',length(p[s]),' packages are successfully installed...\n'))
                    cat(paste('The following packages could not be installed:\n.... ',paste(p[!s],collapse=', '),'\n'))
                  } 
                } else cat(paste('\n ',length(p[s]),' packages are successfully installed...\n'))
              } else cat(paste('\n All required packages have been already installed!\n'))
              
            } else {
              p <- pl[!pl %in% c('stats','utils','parallel','base','grDevice','tools','methods','graphics','compiler','datasets','profile','grid')]
              if (length(p) > 0) {
                .detachPackage(p)
                pi <- p[.is.installed(p)]
                if (length(pi) > 0) pi <- try(remove.packages(pi),silent = TRUE)
                
                s <- rep(TRUE,length(p))
                for (i in seq_along(p)) {
                  pi <- try(install.packages(p[i],...),silent = TRUE)
                  if (inherits(pi, "try-error")) s[i] <- FALSE
                }
                
                if (any(!s)) {
                  if (any(s)) {
                    cat(paste('\n',length(p[s]),' packages are successfully installed or updated...\n'))
                    cat(paste('The following packages could not be installed:\n.... ',paste(p[!s],collapse=', '),'\n'))
                  }
                } else cat(paste('\n ',length(p[s]),' packages are successfully installed or updated...\n'))
              } else cat(paste('\n There is no package to install!\n'))
            }
            .addMethods()
          }
)
