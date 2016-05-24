# gamm/gamm4 support


# `gamm` <-
# function(formula, random = NULL, ..., lme4 = inherits(random, "formula")) {
	# pkg <- if(lme4) "gamm4" else "mgcv"
    # if (!require(pkg, character.only = TRUE)) stop("'gamm' requires package '",
												   # pkg, "' to be installed")
	# cl <- match.call()

	# if(lme4) {
		# pkg <- "gamm4"
		# funcname <- "gamm4"
	# } else {
		# pkg <- "mgcv"
		# funcname <- "gamm"
	# }

	# fun <- call("::", as.name(pkg), as.name(funcname))
	# cl2 <- match.call(eval(fun))
	# cl2$lme4 <- NULL
	# cl2[[1L]] <- fun
	# structure(c(eval(cl2, parent.frame()), list(call = cl)),
			  # class = c(if(lme4) "gamm4", "gamm", "list"))
# }

`update.gamm` <-
function(object, ...) {
# or, if call is as attribute: object$call <- attr(object, "call")
	update.default(object, ...)
}

`print.gamm` <-
function(x, ...) {
	cat("\nCall:\n", paste(asChar(x$call, nlines = -1L), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
	cat("--- \n")
	print(x[[if(inherits(x, "gamm4")) "mer" else "lme"]])
	cat("--- \n")
	print(x$gam)
	invisible(x)
}

`formula.gamm` <-
function (x, ...) formula(x$gam, ...)