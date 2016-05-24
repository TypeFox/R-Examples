
#' @importFrom methods reconcilePropertiesAndPrototype
NULL

install_quietly <- TRUE

with_wd <- function(dir, expr) {
  wd <- getwd()
  on.exit(setwd(wd))
  setwd(dir)
  eval(substitute(expr), envir = parent.frame())
}

#' @importFrom utils tar

build_pkg <- function(path, pkg_file = NULL) {
  if (!file.exists(path)) stop("path does not exist")
  pkg_name <- basename(path)
  if (is.null(pkg_file)) {
    pkg_file <- file.path(dirname(path), paste0(pkg_name, "_1.0.tar.gz"))
  }
  with_wd(dirname(path),
          tar(basename(pkg_file), pkg_name, compression = "gzip"))
  pkg_file
}

#' @importFrom utils package.skeleton install.packages

install_tmp_pkg <- function(..., pkg_name, lib_dir, imports = character()) {
  if (!file.exists(lib_dir)) stop("lib_dir does not exist")
  if (!is.character(pkg_name) || length(pkg_name) != 1) {
    stop("pkg_name is not a string")
  }

  ## Create a directory that will contain the source package
  src_dir <- tempfile()
  on.exit(try(unlink(src_dir, recursive = TRUE), silent = TRUE), add = TRUE)
  dir.create(src_dir)

  ## Create source package, need a non-empty environment,
  ## otherwise package.skeleton fails
  tmp_env <- new.env()
  assign("f", function(x) x, envir = tmp_env)
  suppressMessages(package.skeleton(pkg_name, path = src_dir,
                                    environment = tmp_env))
  pkg_dir <- file.path(src_dir, pkg_name)

  ## Make it installable: remove man, add imports
  unlink(file.path(pkg_dir, "man"), recursive = TRUE)
  if (length(imports) != 0) {
    cat("Imports: ", paste(imports, collapse = ", "), "\n",
        file = file.path(pkg_dir, "DESCRIPTION"), append = TRUE)
    cat(paste0("import(", imports, ")"), sep="\n",
        file = file.path(pkg_dir, "NAMESPACE"), append = TRUE)
  }

  ## Put the code in it, dput is noisy, so we need to redirect it to
  ## temporary file
  exprs <- list(...)
  unlink(file.path(pkg_dir, "R"), recursive = TRUE)
  dir.create(file.path(pkg_dir, "R"))
  code_file <- file.path(pkg_dir, "R", "code.R")
  tmp_file <- tempfile()
  on.exit(try(unlink(tmp_file), silent = TRUE), add = TRUE)
  sapply(exprs, function(x)
         cat(deparse(dput(x, file = tmp_file)),
             file = code_file, append = TRUE, "\n", sep="\n"))

  ## Build it
  pkg_file <- build_pkg(pkg_dir)

  ## Need to unset R_TESTS. This is set during R CMD check
  ## and it messes up things
  if ("R_TESTS" %in% names(Sys.getenv())) {
    R_TESTS <- Sys.getenv("R_TESTS")
    on.exit(Sys.setenv(R_TESTS = R_TESTS), add = TRUE)
    Sys.unsetenv("R_TESTS")
  }

  ## Install it into the supplied lib_dir
  install.packages(pkg_file, lib = lib_dir, repos = NULL, type = "source",
                   quiet = install_quietly)
}

with_libpath <- function(lib_path, ...) {
  cur_lib_path <- .libPaths()
  on.exit(.libPaths(cur_lib_path), add = TRUE)
  .libPaths(c(lib_path, cur_lib_path))
  exprs <- c(as.list(match.call(expand.dots = FALSE)$...))
  sapply(exprs, eval, envir = parent.frame())
}

#' Create, install, load and attach multiple disposable packages
#'
#' If a package with the same name as a disposable one, is
#' loaded, then it will be unloaded. If a package with same name
#' as a disposable on is installed in \code{lib_dir}, then
#' it will be overwritten. (\code{lib_dir} is usually a temporary
#' directory, so this is not a big problem.)
#'
#' @details
#' Note that if you specify \code{lib_dir} and it points to an
#' existing directory, \code{make_package} overwrites the packages
#' there. If an error happens during installation or loading of
#' the disposables packages, then it will \emph{not} restore the
#' original contents of \code{lib_dir}, but it will remove
#' all newly installed disposable packages, even the ones
#' that were installed cleanly.
#'
#' @param ... Named expressions.
#'   A separate package with the given name is created for each.
#' @param lib_dir Directory to install the package to.
#'   Defaults to a temporary directory that is
#'   deleted once the R session is over.
#' @param imports The 'Imports' field in the DESCRIPTION file,
#'   the packages to import in each disposable package.
#' @return A named list with entries: \itemize{
#'     \item \code{lib_dir} The directory in which the packages are
#'       installed.
#'     \item \code{package} The named of the packages.
#'   }
#'
#' @section Examples:
#' \preformatted{
#' pkg <- make_packages(
#'   foo1 = { f <- function() print("hello!") ; d <- 1:10 },
#'   foo2 = { f <- function() print("hello again!") ; d <- 11:20 }
#' )
#' foo1::f()
#' foo2::f()
#' foo1::d
#' foo2::d
#' dispose_packages(pkg)
#' }
#' 
#' @export
#' @seealso \code{\link{dispose_packages}}

make_packages <- function(..., lib_dir = tempfile(),
                          imports = character()) {

  remove_lib_dir <- !file.exists(lib_dir)
  if (remove_lib_dir) dir.create(lib_dir)
  exprs <- c(as.list(match.call(expand.dots = FALSE)$...))

  pkgs <- list(lib_dir = lib_dir, packages = names(exprs))

  ## Start clean
  dispose_packages(pkgs, delete_lib_dir = FALSE)

  ## Clean up on error
  on.exit(dispose_packages(pkgs, delete_lib_dir = remove_lib_dir))

  for (i in seq_along(exprs)) {
    expr <- exprs[[i]]
    name <- names(exprs)[i]
    install_tmp_pkg(expr, pkg_name = name,
                     lib_dir = lib_dir, imports = imports)
    with_libpath(lib_dir, suppressMessages(library(name, quietly = TRUE,
                                                   character.only = TRUE)))
  }

  on.exit()
  pkgs
}

#' Get rid of temporary packages
#'
#' @param packages A list returned by \code{\link{make_packages}}.
#' @param unattach Whether to unattach the packages.
#' @param unload Whether to unload the packages. It is not possible to
#'   unload without unattaching.
#' @param delete Whether to delete the installed packages from the
#'   \code{lib_dir}. If \code{delete_lib_dir} is \code{TRUE}, then
#'   this should be \code{TRUE} as well.
#' @param delete_lib_dir Whether to delete the the whole \code{lib_dir}.
#'
#' @section Examples:
#' \preformatted{
#' pkg <- make_packages(
#'   foo1 = { f <- function() print("hello!") ; d <- 1:10 },
#'   foo2 = { f <- function() print("hello again!") ; d <- 11:20 }
#' )
#'
#' foo1::f()
#' foo2::f()
#' foo1::d
#' foo2::d
#'
#' ## Unattach only
#' dispose_packages(pkg, unload = FALSE, delete = FALSE)
#' "package:foo1" \%in\% search()
#' "foo1" \%in\% loadedNamespaces()
#' dir(pkg$lib_dir)
#'
#' ## Unload
#' dispose_packages(pkg, delete = FALSE)
#' "package:foo1" \%in\% search()
#' "foo1" \%in\% loadedNamespaces()
#' dir(pkg$lib_dir)
#'
#' ## Delete completely
#' dispose_packages(pkg)
#' "package:foo1" \%in\% search()
#' "foo1" \%in\% loadedNamespaces()
#' file.exists(pkg$lib_dir)
#' }
#'
#' @export

dispose_packages <- function(packages, unattach = TRUE, unload = unattach,
                             delete = TRUE, delete_lib_dir = delete) {

  if (!unattach && unload) stop("Cannot unload without unattaching")
  if (!delete && delete_lib_dir) {
    stop("Cannot delete lib_dir without deleting packages")
  }

  if (unattach) {
    for (n in packages$packages) {
      pn <- paste0("package:", n)
      if (pn %in% search()) detach(pn, character.only = TRUE)
    }
  }

  if (unload) {
    for (n in packages$packages) {
      if (n %in% loadedNamespaces()) unloadNamespace(n)
    }
  }

  if (delete_lib_dir) {
    unlink(packages$lib_dir, recursive = TRUE)
  } else if (delete) {
    for (n in packages$packages) {
      unlink(file.path(packages$lib_dir, n), recursive = TRUE)
    }
  }

  invisible()
}
