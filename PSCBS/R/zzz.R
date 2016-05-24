# This should ideally be added to 'R.utils' one day, e.g. as import().
# /HB 2014-06-08
.use <- function(name=NULL, package, ..., mode="function") {
  # Attach package?
  if (is.null(name)) {
    use(package, ...);
    return(invisible(NULL))
  }

  # Retrieve a particular function
  use(package, ..., how="load");
  ns <- getNamespace(package);
  res <- get(name, mode=mode, envir=ns);

  invisible(res);
} # .use()


.onLoad <- function(libname, pkgname) {
  ## covr: skip=5
  ns <- getNamespace(pkgname);
  pkg <- Package(pkgname);
  # Assign '.PSCBS' object [since 'PSCBS' is a constructor/Class].
  name <- sprintf(".%s", pkgname);
  assign(name, pkg, envir=ns, inherits=FALSE);
} # .onLoad()


.onAttach <- function(libname, pkgname) {
  # Copy some pre-memoized CBS-parameter calculations to the 'R.cache'
  # cache.  This speeds up the calculation for common CBS use cases.
  .prememoize();

  # Inform user if DNAcopy is missing
  if (!isPackageInstalled("DNAcopy")) {
    msg <- "The Bioconductor package 'DNAcopy' is not installed. Please see http://www.bioconductor.org/ on how to install it, or try calling PSCBS::installDNAcopy().";
    hrule <- paste(rep("*", times=getOption("width", 80L)-1L), collapse="");
    packageStartupMessage(sprintf("%s\nNOTE: %s\n%s", hrule, msg, hrule));
  }

  # Get '.PSCBS' object [since 'PSCBS' is a constructor/Class].
  name <- sprintf(".%s", pkgname);
  pkg <- get(name, envir=getNamespace(pkgname), inherits=FALSE);
  startupMessage(pkg);
}


############################################################################
# HISTORY:
# 2014-06-08
# o CLEANUP: Using R.utils::use() instead of require().
# 2014-02-04
# o Added .useDNAcopy() to simplify backward compatibility.
# 2013-10-13
# o Added .onLoad() that creates Package '.PSCBS' object, which is
#   used in .onAttach().  This is a workaround for not allocating a
#   local Package on in .onAttach(), which then will be garbage
#   collected and finalize():d, which in turn can generate cyclic
#   loading of namespaces in R.oo (< 1.16.0).
# 2013-09-27
# o Added .useAromaLight() to simplify backward compatibility.
# o Added .requirePkg() from the R.rsp package.
# 2011-07-23
# o Added a namespace to the package.
############################################################################
