Student.family <-
function(object, ...)
UseMethod("family")

print.Student.family <-
function (x, ...) 
cat(" Family:", deparse(x$call), "\n")

Student <-
function(eta = 0.25)
{
  cl <- match.call()
  if ((eta < 0 ) || (eta > 1/2))
      stop("eta must be in [0,1/2)")
  pars <- list(eta = eta)
  pnames <- list(eta = "shape parameter")
  kind <- ifelse(eta == 0, 0, 1)
  structure(list(family = "Student",
                 call = cl,
                 pars = pars,
                 pnames = pnames,
                 npars = 1,
                 kind = kind),
            class = "Student.family")
}
