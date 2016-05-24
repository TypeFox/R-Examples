`func.trim.RPD` <- function(data, trim = 0.25, deriv = c(0,1))
{
  depth.RPD(data, trim = trim, deriv = deriv)$mtrim
}

