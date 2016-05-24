###############################################################################
#' Manage CRAN-like repo on local network
#'
#' This package will allow easier maintainence of CRAN-like repos on local
#' networks (i.e. not on CRAN). This might be necessary where hosted packages
#' contain intellectual property owned by a corporation.
#'
#' @docType package
#' @name managelocalrepo
#' @aliases managelocalrepo package-managelocalrepo
NULL

###############################################################################
#' Determine the full path for local repo
#'
#' This will return the full file path of the terminal directory of a given
#' repo type (and R version number as is appropriate).
#'
#' @param repo_base the path of the base of the repository tree. This will have
#' the following child folders: {repo_base}/bin/ and {repo_base}/src/. A
#' character vector of length one.
#' @param version the version of R that the package should be made available for.
#' This is not relevant when \code{type} is \code{src}. A character vector
#' of length one.
#' @param type should be \code{"win"} (default), \code{"mac"} or \code{src}
#' for Windows binary, Mac binary and source package distributions respectively.
#' @param ... optional arguments to pass to \code{\link{file.path}}
#' @return Full path to the desired local repo's terminal directory.
#' @references
#' \href{http://cran.r-project.org/doc/manuals/R-admin.html}{Setting up a package repository}
#' @examples
#' \dontrun{
#' repo_base <- file.path(".")
#' version <- '3.0'
#' full_repo_dir(repo_base, version)
#' }

full_repo_dir <- function (repo_base, version, type="win", ...)
{
  if (type == 'src') {
    if (!missing(version)) {
      msg <- paste0('The version argument is not needed for type src. ',
                    'It has been ignored.\n')
      message(msg)
    }
    tree <- file.path("src", "contrib", ...)
  } else if (type == 'win') {
    tree <- file.path("bin", "windows", "contrib", version, ...)
  } else if (type == 'mac') {
    tree <- file.path("bin", "macosx")
    version_major_point <- as.numeric(stringr::str_sub(version, end=1))
    if (version_major_point < 3)
      tree <- file.path(tree, "leopard")
    tree <- file.path(tree, "contrib", version, ...)
  } else {
    stop ('Type should be one of: win, mac or src.\n')
  }
  file.path(repo_base, tree, ...)
}

###############################################################################
#' Create terminal repo directory
#'
#' This will create the required terminal directory in the repo.
#'
#' @inheritParams full_repo_dir
#' @return No return value. Will create directory if the directory does not
#' exist. Otherwise it will return an error message.
#' @param ... optional arguments to pass to \code{\link{file.path}}
#' @references
#' \href{http://cran.r-project.org/doc/manuals/R-admin.html}{Setting up a package repository}
#' @examples
#' \dontrun{
#' repo_base <- file.path(".")
#' version <- '3.0'
#' create_terminal_dir(repo_base, version)
#' }

create_terminal_dir <- function (repo_base, version, type="win", ...)
{
  full_dir <- full_repo_dir(repo_base, version, type, ...)
  dir.create(full_dir, recursive=TRUE)
}


###############################################################################
# Check whether package file name meets expectations:
# http://cran.r-project.org/doc/manuals/R-admin.html#Setting-up-a-package-repository

check_package_file_name <- function (package_location, type)
{
  pkg_file_ext <- tools::file_ext(package_location)
  assertthat::assert_that(any(
    all(type == 'src', pkg_file_ext %in% c('gz', 'bz2', 'xz')),
    all(type == 'win', pkg_file_ext == 'zip'),
    all(type == 'mac', pkg_file_ext == 'tgz')))
}

###############################################################################
#' Release a package to the local repo.
#'
#' This will take a package from a given location and populate it to the
#' relevant trees in the local repo. It will create terminal directories if
#' needed.
#'
#' @param package_location the full path to the package's file. This should
#' be the package's file name represented by a character vector of length one.
#' @param repo_base the base directory of the local repository represented by a
#' character vector of length one.
#' @param type the type of package file being released. Should be
#' \code{win} (default), \code{mac} or \code{src} for Windows binary,
#' Mac binary and source package distributions respectively.
#' @param version determines which terminal directory to release the package
#' file to, given a value for \code{type} that is either \code{wim} or
#' \code{mac}. Can also use \code{"all"} to release the package to all terminal
#' directories for either \code{wim} or \code{mac}. If \code{type} is
#' \code{src}, then \code{version} is ignored, as this is not needed.
#' @param r_versions a character vector of R versions. This is used to determine
#' the terminal directories when \code{type} is \code{"win"} or \code{"mac"} and
#' \code{version} is \code{"all"}. Should follow the specifications set out in
#' R Installation and Admin guide.
#' @param ... optional arguments to pass to \code{\link{file.path}}
#' @return No return value. Will release a package to the desired local repo's
#' terminal directory and update the relevant index file. Otherwise will return
#' a suitable error message.
#' @references
#' \href{http://cran.r-project.org/doc/manuals/R-admin.html}{Setting up a package repository}
#' @examples
#' \dontrun{
#' package_location <- file.path(".", "submissions")
#' package_location <- file.path(package_location, "test_package.zip")
#' repo_base <- file.path(".")
#' release_package(package_location, repo_base)
#' }
#' @export

release_package <- function (package_location, repo_base,
                             type='win', version='all',
                             r_versions = c('2.15', '3.0', '3.1', '3.2'), ...)
{
  # Check validity of pkg file ext against what is expected in R-admin
  check_package_file_name(package_location, type)

  # Determine repo trees to populate (i.e. the terminal directories)
  if (type == 'src')
    repo_trees <- full_repo_dir(repo_base, version, type, ...)
  else if (version == 'all')
    repo_trees <- sapply(X=r_versions, FUN=full_repo_dir,
                         repo_base=repo_base, type=type)
  else
    repo_trees <- full_repo_dir(repo_base, version, type, ...)

  # Check if terminal directory/directories exists
  trees_exist <- sapply(repo_trees, file.exists)

  # If it doesn't/they don't exist, create it/them
  lapply(X=repo_trees[!trees_exist], FUN=dir.create, recursive=TRUE)

  # If does exist, copy file across to terminal directories. Overwrites file
  # if it already exists
  lapply(X=repo_trees, FUN=file.copy,
         from=package_location, recursive=TRUE, copy.mode=TRUE)

  # Update the package indices
  type <- switch(type,
    win = "win.binary",
    mac = "mac.binary",
    src = "source")
  lapply(X=repo_trees, FUN=tools::write_PACKAGES, type = type)

  # Return NULL invisibly to avoid write_PACKAGES
  invisible(NULL)
}

###############################################################################
#' Quick release a package
#'
#' This is builds on top of \code{\link{release_package}} and makes it
#' quicker to release packages. The repo base directory can be specified in the
#' \code{managelocalrepo.base} option. The directory in which the package
#' file is located is assumed to be a "submissions" folder and can be specified
#' in the \code{managelocalrepo.submissions} option. These options can be set
#' using \code{\link{options}}, potentially by using \code{\link{.First}}
#'
#' @param file_name the package's file name (not full path)
#' @param ... optional arguments to pass to \code{\link{release_package}}
#' @return No return value. Will release a package to repo's
#' terminal directory and update the index. Otherwise will return a
#' suitable error message.
#' @examples
#' \dontrun{
#' quick_release_package("test_package.zip")
#' }
#' @export

quick_release_package <- function (file_name, ...)
{
  repo_base <- getOption('managelocalrepo.base')
  submissions <- getOption('managelocalrepo.submissions')
  package_location <- file.path(submissions, file_name)
  release_package(package_location, repo_base, ...)
}

