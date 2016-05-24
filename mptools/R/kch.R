#' Extract Metapop KCH
#' 
#' Extract KCH time series for each RAMAS Metapop population.
#' 
#' @param meta The R object holding Metapop metadata returned by
#'   \code{\link{meta}}.
#' @param path A character string giving the full path to the directory
#'   containing .kch files referred to in object \code{meta}.
#' @return A \code{matrix} containing one column per population, giving the
#'   carrying capacity at each time step (i.e. each row).
#' @importFrom utils tail
#' @export
#' @examples
#' mp <- system.file('example.mp', package='mptools')
#' k <- kch(meta=meta(mp), path=dirname(mp))
kch <- function (meta, path) {
  f <- list.files(path, full.names=TRUE, pattern='\\.kch$', ignore.case=TRUE)
  kchs <- path.expand(file.path(path, meta$kch))
  if (any(!file.exists(kchs))) {
    kchs <- f[match(tolower(kchs), tolower(f))]
  } 
  if (any(!file.exists(kchs[!is.na(kchs)])))
    stop('Some file paths do not exist.')
  kch <- lapply(seq_along(kchs), function(x) {
    if (is.na(kchs[x])) meta[x, 'K'] else as.numeric(readLines(kchs[x]))
  })
  len <- sapply(kch, length)
  max.len <- max(len)
  kch <- sapply(kch, function(x) {
    l <- length(x)
    if(l < max.len) c(x, rep(utils::tail(x, 1), max.len - l)) else x
  })
  colnames(kch) <- meta$popName
  kch
}
