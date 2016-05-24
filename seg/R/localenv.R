# ------------------------------------------------------------------------------
# Function 'localenv'
#
# Author: Seong-Yun Hong <hong.seongyun@gmail.com>
# ------------------------------------------------------------------------------
localenv <- function(x, data, power = 2, useExp = TRUE, maxdist, sprel, 
  tol = .Machine$double.eps) {
  
  tmp <- suppressMessages(chksegdata(x, data))
  coords <- tmp$coords; data <- tmp$data; proj4string <- tmp$proj4string
  
  if (missing(maxdist))
    maxdist <- -1
  else if (!is.numeric(maxdist))
    stop("'maxdist' must be numeric", call. = FALSE)
  else if (maxdist < 0)
    stop("'maxdist' must be greater than or equal to 0", call. = FALSE)
  
  if (missing(sprel)) 
    sprel <- coords
  else if (class(sprel) != "nb" && class(sprel) != "dist")
    stop("invalid object 'sprel'", call. = FALSE)

  env <- localenv.get(sprel, data, power, useExp, maxdist, tol)

  SegLocal(coords, data, env, CRS(proj4string))
}
