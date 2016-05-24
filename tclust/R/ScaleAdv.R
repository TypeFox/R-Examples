
.ScaleAdv <- function (x, center, scale)
{
  x <- .Conv2Matrix (x)

  m <- .genapply (x, center, "center", 0)

  if (!is.null (attributes(m)$unknown))
    warning ("Unknown center specification. Centering will be omitted.\n")

  s <- .genapply (x, scale, "scale", 1)
  if (!is.null (attributes(s)$unknown))
    warning ("Unknown scale specification. Scaling will be omitted.\n")

  idx.s.z <- s == 0
  if (any (idx.s.z))
  {
    warning ("Scale vaules == 0 detected. Changing them to 1.")
    s [idx.s.z] <- 1
  }

  x <- t ((t (x) - m) / s)

  return (list (x = x, center = m, scale = s))
}

.genapply <- function (x, arg, arg.name, def)
{

  p <- ncol (x)
  if (missing (arg) || is.null (arg))
    return (rep (def, p))

  if (is.character (arg))
    arg <- eval (parse (text = arg))

  if (is.numeric (arg))
  {
    if (length (arg) == 1)
      return (rep (arg, p))
    if (length (arg) != p)
      stop (paste ("The argument \"", arg.name, "\" has wrong length\n",
                   "    (expected a numerical vector of length 1 or ", p, ")"),
                   sep = "")
    return (arg)
  }

  if (is.function (arg))
  {
    r <- arg (x)
    if (is.numeric (r))
    {
      if (length (r) == 1)
        return (apply (x, 2, arg))
      if (length (r) != p)
        stop (paste ("The function provided by argument \"", arg.name, 
                     "\" returns a vecor of wrong length\n",
                     "    (expected a numerical vector of length 1 or ", p, ")"),
                     sep = "")
      return (r)
    }

    if (is.list (r) &&
        !is.null (r$par))
        return (r$par)
  }

  ret <- rep (def, p)

  attributes (ret)$unknown <- TRUE

  return (ret)
}
