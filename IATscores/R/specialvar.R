specialvar <- function(x, type="standard", tr = .1)
{
  var(trim_or_win(x, type = type, tr = tr, na.rm = TRUE))
}