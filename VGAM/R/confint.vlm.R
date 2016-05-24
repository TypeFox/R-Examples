# 20150827; confint.vlm.R

# Last modified: 20150827,
# 1. 20150827; Tingting Zhan prodded me to do this.
# 2. 
# 3.



confintvglm <- function(object, parm, level = 0.95, ...) {
  cf <- coef(object)
  pnames <- names(cf)
  if (missing(parm))
    parm <- pnames else
  if (is.numeric(parm))
    parm <- pnames[parm]
  a <- (1 - level)/2
  a <- c(a, 1 - a)
# ":::" may lead to warning in \pkg{devtools}:
# pct <- stats:::format.perc(a, 3)
  format.perc <- function(probs, digits)
  paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
        "%")
  pct <- format.perc(a, 3)
  fac <- qnorm(a)
  ci <- array(NA, dim = c(length(parm), 2L), dimnames = list(parm, pct))
  ses <- sqrt(diag(vcov(object)))[parm]
  ci[] <- cf[parm] + ses %o% fac
  ci
}



confintrrvglm <- function(object, parm, level = 0.95, ...) {
  stop("currently this function has not been written")


# 20150828; yettodo: write this function... it's not too difficult.
# Base it on vcovrrvglm() and summaryrrvglm() because they do all
# the work, which involves elements in A and C and B1.
# 

  
}


confintvgam <- function(object, parm, level = 0.95, ...) {
  stop("currently this function has not been written")
}






# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,


# if (FALSE) {
if (!isGeneric("confint"))
    setGeneric("confint",
               function(object, parm, level = 0.95, ...)
               standardGeneric("confint"),
           package = "VGAM")


setMethod("confint", "vglm",
          function(object, parm, level = 0.95, ...)
            confintvglm(object = object, parm = parm, level = level, ...))


# Stop other types of models from working... its dangerous:
setMethod("confint", "rrvglm",
          function(object, parm, level = 0.95, ...)
            confintrrvglm(object = object, parm = parm, level = level, ...))

setMethod("confint", "vgam",
          function(object, parm, level = 0.95, ...)
            confintvgam(object = object, parm = parm, level = level, ...))
# }


# ,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,,


