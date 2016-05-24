vec2spc <- function (x)
{
  tfl2spc(vec2tfl(x)) # could be optimized by calculating spectrum directly: table(table(x))
}
