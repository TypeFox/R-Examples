s4class <- function(x)
{
  cx <- NULL
  for(j in c(".", "_")) {
    if(grepl(paste("_random", j, sep = ""), x, fixed = TRUE))
      cx <- "random.bayesx"
    if(grepl(paste("_spatialtotal", j, sep = ""), x, fixed = TRUE))
      cx <- "mrf.bayesx"
    if(grepl(paste("_pspline", j, sep = ""), x, fixed = TRUE))
      cx <- "sm.bayesx"
    if(grepl(paste("_season", j, sep = ""), x, fixed = TRUE))
      cx <- "sm.bayesx"
    if(grepl(paste("_rw", j, sep = ""), x, fixed = TRUE))
      cx <- "sm.bayesx"
    if(grepl(paste("_spatial", j, sep = ""), x, fixed = TRUE))
      cx <- "mrf.bayesx"
    if(grepl(paste("_geospline", j, sep = ""), x, fixed = TRUE))
      cx <- "geo.bayesx"
    if(grepl(paste("_geokriging", j, sep = ""), x, fixed = TRUE))
      cx <- "geo.bayesx"
    if(grepl(paste("_logbaseline", j, sep = ""), x, fixed = TRUE))
      cx <- "sm.bayesx"
    if(grepl(paste("_kriging", j, sep = ""), x, fixed = TRUE))
      cx <- "sm.bayesx"
    if(grepl(paste("_ridge", j, sep = ""), x, fixed = TRUE))
      cx <- "sm.bayesx"
  }
  if(is.null(cx)) {
    warning(paste("no appropriate class found for:", x))
    cx <- "sm.bayesx"
  }

  return(cx)
}

s4bs <- function(x)
{
  if(grepl("_random", x))
    bs <- "re"
  if(grepl("_pspline", x))
    bs <- "ps"
  if(grepl("_season", x))
    bs <- "season"
  if(grepl("_rw", x))
    bs <- "rw"
  if(grepl("_spatial", x))
    bs <- "mrf"
  if(grepl("_geospline", x))
    bs <- "gs"
  if(grepl("_geokriging", x))
    bs <- "gk"
  if(grepl("_logbaseline", x))
    bs <- "bl"
  if(grepl("_kriging", x))
    bs <- "kr"
  if(grepl("_ridge", x))
    bs <- "ridge"

  return(bs)
}
