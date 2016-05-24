specialmean <- function(x, type = "standard", tr = .1)
{
  mean(trim_or_win(x, type = type, tr = tr, na.rm = TRUE))
}