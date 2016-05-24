## -----------------------------------------------------------------------------
## Vertical Volumetric Weighted Mean of Matter Concentrations in Water Bodies
## -----------------------------------------------------------------------------

vertmean <- function(depth, vari, level, top, bot, vol, total = FALSE) {
  idepth  <- which(depth >= top & depth <= bot)
  ndepth  <- c(depth[idepth])
  nvari   <- c(vari[idepth])
  mdepth  <- c(top, (ndepth[-1] + ndepth[-length(ndepth)])/2, bot)
  nvol    <- vol(mdepth, level)
  dvol    <- abs(diff(nvol))
  ret <- sum(nvari * dvol)
  if (total)
      return(ret)   # return sum instead of mean value
  else
      return(ret / sum(dvol))
}

