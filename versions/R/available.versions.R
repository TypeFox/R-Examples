#' available.versions
#'
#' @description List all of the past versions of the named packages ever
#'  uploaded to CRAN (and therefore in the CRAN source archives), their
#'  publication dates and whether they can be installed from MRAN via
#'  \code{\link{install.versions}} or \code{\link{install.dates}}.
#'
#' @param pkgs character vector of the names of packages for which to query
#' available versions
#'
#' @return a list of dataframes, each giving the versions and publication dates
#'  for the corresponding elements of \code{pkgs} as well as whether they can be
#'  installed from MRAN
#'
#' @export
#' @name available.versions
#' @examples
#'
#' \dontrun{
#'
#' # available versions of checkpoint
#' available.versions('checkpoint')
#'
#' # available versions of checkpoint and devtools
#' available.versions(c('checkpoint', 'devtools'))
#'
#' }
#'
available.versions <- function (pkgs) {

  # vectorise by recursion
  if (length(pkgs) > 1) {
    ans <- lapply(pkgs,
                  available.versions)

    # remove a level of listing
    ans <- lapply(ans, '[[', 1)

    names(ans) <- pkgs

    return (ans)
  }

  # get the current version
  current_df <- current.version(pkgs)

  # see if the package has been archived

  # get most recent MRAN image URL, Archive directory
  archive_url <- sprintf('%s/src/contrib/Archive',
                         latest.MRAN())

  # check for the package
  archived <- pkg.in.archive(archive_url, pkgs)

  # if it is archived, get the previous versions
  if (archived) {

    # get most recent MRAN image URL for the package archive
    # (inside the Archive directory)
    pkg_archive_url <- sprintf('%s/src/contrib/Archive/%s',
                           latest.MRAN(),
                           pkgs)

    # scrape the versions therein
    previous_df <- scrape.index.versions(pkg_archive_url,
                                         pkgs)

  } else {

    # otherwise, make it a blank row
    previous_df <- current_df[0, ]

  }

  # append previous versions to the current version
  df <- rbind(current_df,
              previous_df)

  # add whether they are on MRAN
  df$available <- as.Date(df$date) >= as.Date('2014-09-17')

  # wrap into a list
  ans <- list()

  ans[[pkgs]] <- df

  return(ans)

}
