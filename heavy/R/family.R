heavy.family <-
function(object, ...)
UseMethod("family")

print.heavy.family <-
function (x, ...) 
cat(" Family:", deparse(x$call), "\n")

normal <-
function()
{
  cl <- match.call()
  pars <- NULL
  pnames <- NULL
  structure(list(family = "Normal",
                 call = cl,
                 pars = pars,
                 pnames = pnames,
                 npars = 0,
                 which = 0),
            class = "heavy.family")
}

Cauchy <-
function()
{
  cl <- match.call()
  pars <- NULL
  pnames <- NULL
  structure(list(family = "Cauchy",
                 call = cl,
                 pars = pars,
                 pnames = pnames,
                 npars = 0,
                 which = 1),
            class = "heavy.family")
}

Student <-
function(df = 4)
{
  cl <- match.call()
  if (df <= 0)
      stop("Deg. of freedom must be > 0")
  pars <- list(df = df)
  pnames <- list(df = "Deg. of Freedom")
  structure(list(family = "Student",
                 call = cl,
                 pars = pars,
                 pnames = pnames,
                 npars = 1,
                 which = 2),
            class = "heavy.family")
}

slash <-
function(df = 2)
{
  cl <- match.call()
  if (df <= 0)
      stop("Deg. of freedom must be > 0")
  pars <- list(df = df)
  pnames <- list(df = "Deg. of Freedom")
  structure(list(family = "slash",
                 call = cl,
                 pars = pars,
                 pnames = pnames,
                 npars = 1,
                 which = 3),
            class = "heavy.family")
}

contaminated <-
function(epsilon = 0.05, vif = 0.25)
{
  cl <- match.call()
  if ((epsilon < 0) || (epsilon > 1))
      stop("contamination percentage must be in [0,1]")
  if ((vif <= 0) || (vif >= 1))
      stop("variance inflation factor must be in (0,1)")
  pars <- list(epsilon = epsilon, vif = vif)
  pnames <- list(epsilon = "contamination percentage", vif = "variance inflation factor")
  structure(list(family = "contaminated",
                 call = cl,
                 pars = pars,
                 pnames = pnames,
                 npars = 2,
                 which = 4),
            class = "heavy.family")
}
