#  setClass for "bifd"

# setClass("bifd",  representation(coefs     = "array",
#                                  sbasis    = "basisfd",
#                                  tbasis    = "basisfd",
#                                  bifdnames = "list"))

#  Generator function of class bifd

bifd <- function (coef=matrix(0,2,1), sbasisobj=create.bspline.basis(),
                  tbasisobj=create.bspline.basis(), fdnames=defaultnames)
{
  #  This function creates a bivariate functional data object.
  #    A bivariate functional data object consists of two bases for expanding
  #      a bivariate functional observation and a set of coefficients defining
  #      this expansion.
  #    The bases are contained in "basisfd" objects; that is, a realization
  #    of the "basisfd" class.

  #  Arguments:
  #  COEF     ... a two-, three-, or four-dimensional array containing
  #               coefficient values for the expansion of each set of bivariate
  #               function values=terms of a set of basis function values
  #               If COEF is a two-way, it is assumed that there is only
  #                 one variable and only one replication, and then
  #                 the first and second dimensions correspond to
  #                 the basis functions for the first and second argument,
  #                 respectively.
  #               If COEF is a three-way, it is assumed that there is only
  #                 one variable per replication, and then
  #                 the first and second dimensions correspond to
  #                 the basis functions for the first and second argument,
  #                 respectively, and the third dimension corresponds to
  #                 replications.
  #               If COEF is a four-way array, then the fourth dimension
  #                 corresponds to variables
  #  SBASISOBJ ... a functional data basis object for the first  argument s
  #  TBASISOBJ ... a functional data basis object for the second argument t
  #  BIFDNAMES ... A list of length 3 with members containing
  #               1. a single name for the argument domain, such as 'Time'
  #               2. a name for the replications or cases
  #               3. a name for the function.

  #  Returns:
  #  BIFDOBJ ... a bivariate functional data object
#  last modified 2007 May 3 by Spencer Graves
  #  previously modified 20 September 2005

  #  check COEF and get its dimensions

    if(!is.numeric(coef)) stop(
		"coef must be numerical vector or matrix")
    else if (is.vector(coef)) stop(
		"Argument COEF is not at least 2 dimensional.")
    else if (is.matrix(coef)) {
            coefd <- dim(coef)
            ndim  <- length(coefd)
        }
    else if (is.array(coef)) {
            coefd <- dim(coef)
            ndim  <- length(coefd)
        }
    else stop("argument COEF is not correct")

    if (ndim > 4) stop(
        "First argument not of dimension 2, 3 or 4.")

    #  check SBASISOBJ

    if (!inherits(sbasisobj, "basisfd")) stop(
        "Argument SBASISOBJ is not of basis class")

    if (dim(coef)[1] != sbasisobj$nbasis) stop(
        paste("Number of coefficients does not match number of ",
              "basis functions for SBASISOBJ."))

    #  check TBASISOBJ

    if (!inherits(tbasisobj, "basisfd")) stop(
        "Argument TBASISOBJ is not of basis class.")

    if (dim(coef)[2] != tbasisobj$nbasis) stop(
        paste("Number of coefficients does not match number of ",
              "basis functions for TBASISOBJ."))

    #  setup number of replicates and number of variables

    if (ndim > 2) nrep <- coefd[3] else nrep <- 1
    if (ndim > 3) nvar <- coefd[4] else nvar <- 1

    #  set up default fdnames

    if (ndim == 2) defaultnames <- list("time", "reps", "values")
    if (ndim == 2) defaultnames <- list("time",
                                        paste("reps",as.character(1:nrep)),
                                        "values")
    if (ndim == 4) defaultnames <- list("time",
                                        paste("reps",as.character(1:nrep)),
                                        paste("values",as.character(1:nvar)) )

    names(defaultnames) <- c("args", "reps", "funs")

#  S4 definition
#	bifdobj <- new("bifd", coefs=coef, sbasis=sbasisobj, tbasis=tbasisobj,
#	               bifdnames=fdnames)

#  S3 definition

	bifdobj <- list(coefs=coef, sbasis=sbasisobj, tbasis=tbasisobj,
	               bifdnames=fdnames)
	oldClass(bifdobj) <- "bifd"
	
	bifdobj
}

#  "show" method for "bifd"

print.bifd <- function(x, ...)
{
  object <- x	
	cat("bifd:\n\n")
	
	cat("Dimensions of the data:\n")
	  cat(paste("  ",object$fdnames[[1]],"\n"))
	  cat(paste("  ",object$fdnames[[2]],"\n"))
	  cat(paste("  ",object$fdnames[[3]],"\n"))
	  cat("\n")
	
	print(object$sbasis)
	
	print(object$tbasis)
	
}

#  "summary" method for "bifd"

summary.bifd <- function(object, ...)
{
	
	cat("bifd:\n\n")
	
	cat("Dimensions of the data:\n")
	  cat(paste("  ",object$fdnames[[1]],"\n"))
	  cat(paste("  ",object$fdnames[[2]],"\n"))
	  cat(paste("  ",object$fdnames[[3]],"\n"))
	  cat("\n")
	
	print(object$sbasis)
	
	print(object$tbasis)
	
	cat("\nCoefficient matrix:\n\n")
	
	object$coefs
	
}

