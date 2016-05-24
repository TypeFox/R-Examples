time365fromDate <- function(x)
{
  myd <- month.day.year(x)
  cycle <- match(with(myd, paste(month, day)),
                 with(.time365md, paste(month, day)))
  return(myd$year + (cycle - 1)/365)
}
