trim <- function(s)
{
  s <- sub(pattern="^ +", replacement="", x=s)
  s <- sub(pattern=" +$", replacement="", x=s)
  s
}