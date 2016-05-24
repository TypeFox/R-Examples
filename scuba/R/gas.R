#
# 	gas.R
#
#	$Revision: 1.3 $	$Date: 2008/06/14 10:43:09 $
#
###################################################################
#  Gases
#

is.gas <- function(x) { inherits(x, "gas") }

nitrox <- function(fO2) {
  if(fO2 <= 0 || fO2 > 1)
    stop("fO2 should be a fraction between 0 and 1")
  g <- list(fO2=fO2, fN2=1-fO2, fHe=0)
  class(g) <- c("gas", class(g))
  return(g)
}

trimix <- function(fO2, fHe, fN2) {
  nmiss <- missing(fO2) + missing(fN2) + missing(fHe)
  if(nmiss > 1)
    stop("Must specify at least two of the gas fractions")
  if(nmiss == 1) {
    if(missing(fO2)) fO2 <- 1 - fN2 - fHe
    else if(missing(fN2)) fN2 <- 1 - fO2 - fHe
    else fHe <- 1 - fO2 - fN2
  }
  stopifnot(fO2 >= 0 && fN2 >= 0 && fHe >= 0)
  if(abs(1 - (fO2+fN2+fHe)) > 0.00001)
    stop("gas fractions must add up to 1")
  g <- list(fO2=fO2, fN2=fN2, fHe=fHe)
  class(g) <- c("gas", class(g))
  return(g)
}

air <- nitrox(0.21)

is.nitrox <- function(g) { is.gas(g) && with(g, fHe == 0) }

is.air <- function(g) { is.nitrox(g) && with(g, fO2 == 0.21) }

as.character.gas <- function(x, ...) {
  with(x, gasnames(fO2, fN2, fHe, ...))
}
  
gasnames <- function(fO2, fN2=1-fO2, fHe=1 - fO2 - fN2, ..., full=FALSE) {
  # vectorised
  if(!missing(fO2) && !missing(fN2) && !missing(fHe)) {
    if(any(fO2 + fN2 + fHe != 1))
      warning("fO2+fN2+fHe is not equal to 1")
  } else {
    if(missing(fN2) && !missing(fHe) && !missing(fO2))
      fN2 = 1 - fO2 - fHe
    if(missing(fO2) && !missing(fN2) && !missing(fHe))
      fO2 = 1 - fN2 - fHe
  }
  stopifnot(length(fO2) == length(fN2) && length(fN2) == length(fHe))
  isnitrox <- (fHe == 0)
  isheliox <- (fN2 == 0)
  isair    <- isnitrox & (fO2 == 0.21)
  ispureOx <- isnitrox & (fO2 == 1)
  out <-
    ifelse(isair, "air",
         ifelse(ispureOx, "100% O2",
                ifelse(isnitrox, paste("EANx ", 100 * fO2, sep=""),
                       paste(ifelse(isheliox, "heliox", "trimix"),
                             paste(100 * fO2, "/", 100 * fHe, sep="")))))
  if(!full) return(out)
  details <- 
    ifelse(isair, " (21% oxygen, 79% nitrogen)",
         ifelse(ispureOx, "",
                ifelse(isnitrox,
                       paste("(", 100 * fO2, "% oxygen, ",
                                  100 * fN2, "% nitrogen)", sep=""),
                       paste("(", 100 * fO2, "% oxygen, ",
                                  100 * fN2, "% nitrogen, ", 
                                  100 * fHe, "% helium)", sep=""))))
  return(paste(out, details))
}

print.gas <- function(x, ...) {
  cat(paste(as.character(x), "\n"))
  invisible(NULL)
}

summary.gas <- function(object, ...) {
  mo <- mod(object, ppO2max=1.4)
  mo <- round(mo, 1)
  z <- list(name=as.character(object), mod=mo)
  z <- append(z, object[c("fO2", "fN2", "fHe")])
  class(z) <- c("summary.gas", class(z))
  return(z)
}

print.summary.gas <- function(x, ...) {
  blah <- with(x, gasnames(fO2, fN2, fHe, full=TRUE))
  cat(paste(blah, "\n"))
  cat(paste("Maximum operating depth", x$mod, "metres\n"))
  return(invisible(NULL))
}


