#'@title List the version of cereal on github.
#'@return package_version. The vector of versions of cereal.
#'@details Use the github API to query the versions of cereal.
#'  The github page of cereal is \url{https://github.com/USCiLab/cereal}
#'  This function requires the package \code{httr}.
#'@export
list_version <- function() {
  response <- httr::GET("https://api.github.com/repos/USCiLab/cereal/tags")
  stopifnot(httr::status_code(response) == 200)
  v <- sapply(httr::content(response), function(obj) {
    substring(obj$name, 2, nchar(obj$name))
  })
  package_version(v)
}

#'@title Return the latest version of cereal on the github.
#'@return package_version. The latest version of cereal.
#'@details This function requires the package \code{httr}.
#'@export
last_version <- function() {
  max(list_version())
}

.package_file <- function(...) {
  system.file(..., package = .packageName)
}

.path <- "https://github.com/USCiLab/cereal.git"

#'@title Update the cereal to the Given Version
#'@param version package_version. The version of cereal
#'  you want to switch in Rcereal.
#'@details This function will change the version of cereal
#'  in Rcereal. The package \code{git2r} is required.
#'@export
update_version <- function(version = last_version()) {
  stopifnot(length(version) == 1)
  if (class(version)[1] == "character") {
    if (substring(version, 1, 1) == "v") {
      version <- substring(version, 2, nchar(version))
    }
    version <- package_version(version)
  }
  stopifnot(class(version)[1] == "package_version")
  .tmppath <- tempfile()
  dir.create(.tmppath)
  git2r::clone(.path, .tmppath)
  .repo <- git2r::repository(.tmppath)
  .commit <- git2r::tags(.repo)[[sprintf("v%s", version)]]
  git2r::checkout(.commit)
  .dst <- file.path(.package_file(""), "include")
  if (!file.exists(.dst)) dir.create(.dst)
  stopifnot(file.rename(.dst, .include <- file.path(.package_file(""), ".include")))
  tryCatch({
    .src <- file.path(.tmppath, "include")
    file.copy(.src, .package_file(""), overwrite = TRUE, recursive = TRUE)
    stopifnot(dir(.dst) == "cereal")
    unlink(.include, recursive = TRUE)
  }, error = function(e) {
    file.rename(.include, .dst)
    stop(conditionMessage(e))
  })
}
