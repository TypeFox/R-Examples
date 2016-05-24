# utility functions for versions package

# read lines from a url more quickly and with a clearer error
# message on failure than readLines
url.lines <- function (url) {

  # create a tempfile
  file <- tempfile()

  # stick the html in there
  suppressWarnings(success <- download.file(url, file,
                                            quiet = TRUE))

  # if it failed, issue a nice error
  if (success != 0) {
    stop(sprintf('URL does not appear to exist: %s',
                 url))
  }

  # get the lines, delete the file and return
  lines <- readLines(file)
  file.remove(file)
  return (lines)
}


# return the url for the latest date on an index page of dates
# (by default the MRAN snappshot index page)
latest.MRAN <- function(url = 'https://mran.revolutionanalytics.com/snapshot') {
  # get all full dates
  dates <- scrape.index.dates(url)

  # get latest
  max <- as.character(max(as.Date(dates)))

  # form the url and return
  ans <- paste(url, max, sep = '/')
  return (ans)

}


# list the dates in an index page file dates as subdirectories
scrape.index.dates <- function (url) {

  # get the lines
  lines <- url.lines(url)

  # keep only lines starting with hrefs
  lines <- lines[grep('^<a href="*', lines)]

  # take the sequence after the href that is between the quotes
  lines <- gsub('.*href=\"([^\"]+)\".*', '\\1', lines)

  # remove the trailing slash
  lines <- gsub('/$', '', lines)

  # remove any lines that aren't 10 characters long (a date only)
  lines <- lines[nchar(lines) == 10]

  # return list in reverse
  return (rev(lines))

}


# list the package versions in an index page
scrape.index.versions <- function (url, pkgs) {

  # get the lines
  lines <- url.lines(url)

  # keep only lines starting with hrefs
  lines <- lines[grep('^<a href="*', lines)]

  # take the sequence after the href that is between the quotes
  versions <- gsub('.*href=\"([^\"]+)\".*', '\\1', lines)

  # remove the leading package name
  versions <- gsub(sprintf('^%s_', pkgs),
                   '', versions)

  # remove the trailing tarball extension
  versions <- gsub('.tar.gz$', '', versions)

  # match the sequence in number-letter-number format
  dates <- gsub('.*  ([0-9]+-[a-zA-Z]+-[0-9]+) .*', '\\1', lines)

  # convert dates to standard format
  dates <- as.Date(dates, format = '%d-%b-%Y')

  # get them in date order
  o <- order(dates, decreasing = TRUE)

  # create dataframe, reversing both
  df <- data.frame(version = versions[o],
                   date = as.character(dates[o]),
                   stringsAsFactors = FALSE)

  return (df)

}


# given the url to an archive ('.../src/contrib/Archive'), a package name
# and version, see if the package is present and return a scalar logical
pkg.in.archive <- function (url, pkg) {

  # get the lines
  lines <- url.lines(url)

  # keep only lines starting with hrefs
  lines <- lines[grep('^<a href="*', lines)]

  # take the sequence after the href that is between the quotes
  items <- gsub('.*href=\"([^\"]+)\".*', '\\1', lines)

  # expected directory name
  dir <- paste0(pkg, '/')

  # search for the expected package directory
  archived <- dir %in% items

  return (archived)

}


# given packages name and required versions,
# return a date when it was live on CRAN
version2date <- function (pkgs, versions) {

  # vectorise by recursion
  if (length(pkgs) > 1) {
    ans <- mapply(version2date,
                  pkgs,
                  versions)

    return (ans)
  }

  # get available versions for the package
  df <- available.versions(pkgs)[[1]]

  # error if the version is not recognised
  if (!(versions %in% df$version)) {
    stop (sprintf('%s does not appear to be a valid version of %s.
                  Use available.versions("%s") to get valid versions\n\n',
                  versions,
                  pkgs,
                  pkgs))
  }

  # find the row corresponding to the version
  idx <- match(versions, df$version)

  # error if the version is recognised, but not available on MRAN
  if (!df$available[idx]) {
    stop (sprintf("%s is a valid version of %s, but was published before
                  2014-09-17 and can therefore not be downloaded from MRAN.
                  Try using devtools::install_version to install the package
                  from its source in the CRAN archives\n\n",
                  versions,
                  pkgs))
  }

  # get a middling date

  # append today's date (note idx is one off now)
  dates <- c(Sys.Date(),
             as.Date(df$date))

  # get the mean of the publication date and subsequent publication date
  # (or today) as the target date for version installation
  date <- as.character(mean(dates[idx + 0:1]))

  # return this
  return (date)

}


# get current version of package
current.version <- function (pkg) {

  # get all current contributed packages in latest MRAN
  current_url <- sprintf('%s/src/contrib',
                         latest.MRAN())

  # get the lines
  lines <- url.lines(current_url)

  # keep only lines starting with hrefs
  lines <- lines[grep('^<a href="*', lines)]

  # take the sequence after the href that is between the quotes
  tarballs <- gsub('.*href=\"([^\"]+)\".*', '\\1', lines)

  # match the sequence in number-letter-number format
  dates <- gsub('.*  ([0-9]+-[a-zA-Z]+-[0-9]+) .*', '\\1', lines)

  # convert dates to standard format
  dates <- as.Date(dates, format = '%d-%b-%Y')

  # get the ones matching the package
  idx <- grep(sprintf('^%s_.*.tar.gz$', pkg),
              tarballs)

  if (length(idx) == 1) {
    # if this provided exactly one match, it's the current package
    # so scrape the version and get the date

    versions <- tarballs[idx]

    # remove the leading package name
    versions <- gsub(sprintf('^%s_', pkg),
                     '', versions)

    # remove the trailing tarball extension
    versions <- gsub('.tar.gz$', '', versions)

    dates <- dates[idx]

  } else {
    # otherwise warn and return NAs
    warning (sprintf('The current version and publication date of %s could not
                     be detected',
                     pkg))
    versions <- dates <- NA
  }

  # create dataframe, reversing both
  df <- data.frame(version = versions,
                   date = as.character(dates),
                   stringsAsFactors = FALSE)

  return (df)

}
