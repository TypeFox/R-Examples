### xtable package
###
### Produce LaTeX and HTML tables from R objects.
###
### Copyright 2000-2013 David B. Dahl <dahl@stat.byu.edu>
###
### This file is part of the `xtable' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 2, or at your option, any later version,
### incorporated herein by reference.
###
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
###
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA

xtable <- function(x, caption = NULL, label = NULL, align = NULL,
                   digits = NULL, display = NULL, auto = FALSE, ...) {
  UseMethod("xtable")
}


### data.frame and matrix objects

xtable.data.frame <- function(x, caption = NULL, label = NULL, align = NULL,
                              digits = NULL, display = NULL, auto = FALSE,
                              ...) {
  logicals <- unlist(lapply(x, is.logical))
  ##x[, logicals] <- lapply(x[, logicals], as.character)
  ## Patch for logicals bug, no 1911
  ## David Scott, <d.scott@auckland.ac.nz>, 2012-08-10
  x[, logicals] <- lapply(x[, logicals, drop = FALSE], as.character)
  characters <- unlist(lapply(x, is.character))
  factors <- unlist(lapply(x, is.factor))
  ints <- sapply(x, is.integer)
  class(x) <- c("xtable","data.frame")
  caption(x) <- caption
  label(x) <- label
  if(auto && is.null(align))   align   <- xalign(x)
  if(auto && is.null(digits))  digits  <- xdigits(x)
  if(auto && is.null(display)) display <- xdisplay(x)
  align(x) <- switch(1+is.null(align), align,
                     c("r",c("r","l")[(characters|factors)+1]))
  digits(x) <- switch(1+is.null(digits), digits, c(0,rep(2,ncol(x))))
  ## Patch from Seth Falcon <sfalcon@fhcrc.org>, 18-May-2007
  if (is.null(display)) {
      display <- rep("f", ncol(x))
      display[ints] <- "d"
      display[characters | factors] <- "s"
      display <- c("s", display)
  }
  display(x) <- display
  return(x)
}

xtable.matrix <- function(x, caption = NULL, label = NULL, align = NULL,
                          digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.data.frame(data.frame(x, check.names = FALSE),
                           caption = caption, label = label, align = align,
                           digits = digits, display = display, auto = auto,
                           ...))
}



### table objects (of 1 or 2 dimensions) by Guido Gay, 9 Feb 2007
### Fixed to pass R checks by DBD, 9 May 2007
xtable.table <- function(x, caption = NULL, label = NULL, align = NULL,
                       digits = NULL, display = NULL, auto = FALSE, ...) {
  if (length(dim(x)) == 1) {
    return(xtable.matrix(matrix(x,
                                dimnames = list(rownames(x),
                                                names(dimnames(x)))),
                         caption = caption, label = label, align = align,
                         digits = digits, display = display, auto = auto, ...))
  } else if (length(dim(x))==2) {
    return(xtable.matrix(matrix(x, ncol = dim(x)[2], nrow = dim(x)[1],
                                dimnames = list(rownames(x), colnames(x))),
                         caption = caption, label = label, align = align,
                         digits = digits, display = display, auto = auto, ...))
  } else {
    stop("xtable.table is not implemented for tables of > 2 dimensions")
  }
}


### anova objects
xtable.anova <- function(x, caption = NULL, label = NULL, align = NULL,
                         digits = NULL, display = NULL, auto = FALSE, ...) {
  suggested.digits <- c(0,rep(2, ncol(x)))
  suggested.digits[grep("Pr\\(>", names(x))+1] <- 4
  suggested.digits[grep("P\\(>", names(x))+1] <- 4
  suggested.digits[grep("Df", names(x))+1] <- 0

  class(x) <- c("xtable","data.frame")
  caption(x) <- caption
  label(x) <- label
  if(auto && is.null(align))   align   <- xalign(x)
  if(auto && is.null(digits))  digits  <- xdigits(x)
  if(auto && is.null(display)) display <- xdisplay(x)
  align(x) <- switch(1+is.null(align), align, c("l",rep("r", ncol(x))))
  digits(x) <- switch(1+is.null(digits), digits, suggested.digits)
  display(x) <- switch(1+is.null(display), display, c("s",rep("f", ncol(x))))
  return(x)
}


### aov objects
xtable.aov <- function(x, caption = NULL, label = NULL, align = NULL,
                       digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.anova(anova(x, ...), caption = caption, label = label,
                      align = align, digits = digits, display = display,
                      auto = auto, ...))
}

xtable.summary.aov <- function(x, caption = NULL, label = NULL, align = NULL,
                               digits = NULL, display = NULL, auto = FALSE,
                               ...) {
  return(xtable.anova(x[[1]], caption = caption, label = label, align = align,
                      digits = digits, display = display, auto = auto, ...))
}

xtable.summary.aovlist <- function(x, caption = NULL, label = NULL,
                                   align = NULL, digits = NULL, display = NULL,
                                   auto = FALSE, ...) {
    for (i in 1:length(x)) {
        if (i == 1) {
            result <- xtable.summary.aov(x[[i]], caption = caption,
                                         label = label,
                                         align = align, digits = digits,
                                         display = display, auto = auto, ...)
        } else {
            result <- rbind(result,
                            xtable.anova(x[[i]][[1]], caption = caption,
                                         label = label, align = align,
                                         digits = digits, display = display,
                                         auto = auto, ...))
        }
    }
    return(result)
}

xtable.aovlist <- function(x, caption = NULL, label = NULL, align = NULL,
                           digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.summary.aovlist(summary(x), caption = caption, label = label,
                                align = align, digits = digits,
                                display = display, auto = auto, ...))
}



### lm objects
xtable.lm <- function(x, caption = NULL, label = NULL, align = NULL,
                      digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.summary.lm(summary(x), caption = caption, label = label,
                           align = align, digits = digits, display = display,
                           auto = auto, ...))
}

xtable.summary.lm <- function(x, caption = NULL, label = NULL, align = NULL,
                              digits = NULL, display = NULL, auto = FALSE,
                              ...) {
  x <- data.frame(x$coef, check.names = FALSE)

  class(x) <- c("xtable","data.frame")
  caption(x) <- caption
  label(x) <- label
  if(auto && is.null(align))   align   <- xalign(x)
  if(auto && is.null(digits))  digits  <- xdigits(x)
  if(auto && is.null(display)) display <- xdisplay(x)
  align(x) <- switch(1+is.null(align), align, c("r","r","r","r","r"))
  digits(x) <- switch(1+is.null(digits), digits, c(0,4,4,2,4))
  display(x) <- switch(1+is.null(display), display, c("s","f","f","f","f"))
  return(x)
}


### glm objects
xtable.glm <- function(x, caption = NULL, label = NULL, align = NULL,
                       digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.summary.glm(summary(x), caption = caption,
                            label = label, align = align,
                            digits = digits, display = display,
                            auto = auto, ...))
}

xtable.summary.glm <- function(x, caption = NULL, label = NULL, align = NULL,
                               digits = NULL, display = NULL, auto = FALSE,
                               ...) {
  return(xtable.summary.lm(x, caption = caption, label = label, align = align,
                           digits = digits, display = display,
                           auto = auto, ...))
}


### prcomp objects
xtable.prcomp <- function(x, caption = NULL, label = NULL, align = NULL,
                          digits = NULL, display = NULL, auto = FALSE, ...) {
  x <- data.frame(x$rotation, check.names = FALSE)

  class(x) <- c("xtable","data.frame")
  caption(x) <- caption
  label(x) <- label
  if(auto && is.null(align))   align   <- xalign(x)
  if(auto && is.null(digits))  digits  <- xdigits(x)
  if(auto && is.null(display)) display <- xdisplay(x)
  align(x) <- switch(1+is.null(align), align, c("r",rep("r", ncol(x))))
  digits(x) <- switch(1+is.null(digits), digits, c(0,rep(4, ncol(x))))
  display(x) <- switch(1+is.null(display), display, c("s",rep("f", ncol(x))))
  return(x)
}

xtable.summary.prcomp <- function(x, caption = NULL, label = NULL, align = NULL,
                                  digits = NULL, display = NULL, auto = FALSE,
                                  ...) {
  x <- data.frame(x$importance, check.names = FALSE)

  class(x) <- c("xtable","data.frame")
  caption(x) <- caption
  label(x) <- label
  if(auto && is.null(align))   align   <- xalign(x)
  if(auto && is.null(digits))  digits  <- xdigits(x)
  if(auto && is.null(display)) display <- xdisplay(x)
  align(x) <- switch(1+is.null(align), align, c("r",rep("r", ncol(x))))
  digits(x) <- switch(1+is.null(digits), digits, c(0,rep(4, ncol(x))))
  display(x) <- switch(1+is.null(display), display, c("s",rep("f", ncol(x))))
  return(x)
}


### Slightly modified version of xtable.coxph contributed on r-help by
###   Date: Wed, 2 Oct 2002 17:47:56 -0500 (CDT)
###   From: Jun Yan <jyan@stat.wisc.edu>
###   Subject: Re: [R] xtable for Cox model output
xtable.coxph <- function (x, caption = NULL, label = NULL, align = NULL,
                          digits = NULL, display = NULL, auto = FALSE, ...)
{
  cox <- x
  beta <- cox$coef
  se <- sqrt(diag(cox$var))
  if (is.null(cox$naive.var)) {
    tmp <- cbind(beta, exp(beta), se, beta/se, 1 - pchisq((beta/se)^2, 1))
    dimnames(tmp) <- list(names(beta),
      c("coef", "exp(coef)", "se(coef)", "z", "p"))
  } else {
    tmp <- cbind( beta, exp(beta), se, beta/se,
      signif(1 - pchisq((beta/se)^2, 1), digits - 1))
    dimnames(tmp) <- list(names(beta),
      c("coef", "exp(coef)", "robust se", "z", "p"))
  }
  return(xtable(tmp, caption = caption, label = label, align = align,
                digits = digits, display = display, auto = auto, ...))
}

### Additional method: xtable.ts
### Contributed by David Mitchell (davidm@netspeed.com.au)
### Date: July 2003
xtable.ts <- function(x, caption = NULL, label = NULL, align = NULL,
                      digits = NULL, display = NULL, auto = FALSE, ...) {
  if (inherits(x, "ts") && !is.null(ncol(x))) {
    ## COLNAMES <- paste(colnames(x));
    tp.1 <- trunc(time(x))
    tp.2 <- trunc(cycle(x))
    day.abb <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
    ROWNAMES <- switch(frequency(x),
                       tp.1,
                       "Arg2", "Arg3",              # Dummy arguments
                       paste(tp.1, c("Q1", "Q2", "Q3", "Q4")[tp.2], sep = " "),
                       "Arg5", "Arg6",
                       paste("Wk.", tp.1, " ", day.abb[tp.2], sep = ""),
                       "Arg8", "Arg9", "Arg10", "Arg11",
                       paste(tp.1, month.abb[tp.2], sep = " "))
    tmp <- data.frame(x, row.names = ROWNAMES);
  } else if (inherits(x, "ts") && is.null(ncol(x))) {
    COLNAMES <- switch(frequency(x),
                       "Value",
                       "Arg2", "Arg3",              # Dummy arguments
                       c("Q1", "Q2", "Q3", "Q4"),
                       "Arg5", "Arg6",
                       day.abb,
                       "Arg8", "Arg9", "Arg10", "Arg11",
                       month.abb)
    ROWNAMES <- seq(from = start(x)[1], to = end(x)[1])
    tmp <- data.frame(matrix(c(rep(NA, start(x)[2] - 1), x,
                               rep(NA, frequency(x) - end(x)[2])),
                             ncol = frequency(x), byrow = TRUE),
                      row.names = ROWNAMES)
    names(tmp) <- COLNAMES
  }
  return(xtable(tmp, caption = caption, label = label, align = align,
                digits = digits, display = display, auto = auto, ...))
}

### Suggested by Ajay Narottam Shah <ajayshah@mayin.org> in e-mail 2006/07/22
xtable.zoo <- function(x, caption = NULL, label = NULL, align = NULL,
                       digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable(as.ts(x), caption = caption, label = label,
                align = align, digits = digits,
                display = display, auto = auto, ...))
}

### Date: Fri, 29 May 2015 11:41:04 +0200
### From: Martin G. <martin.gubri@framasoft.org>
### Subject: [xtable] Code for spdep, splm and sphet objects outputs
### package spdep
### sarlm objects
xtable.sarlm <- function(x, caption = NULL, label = NULL, align = NULL,
                         digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.summary.sarlm(summary(x), caption = caption, label = label,
                              align = align, digits = digits,
                              display = display, auto = auto, ...))
}

xtable.summary.sarlm <- function(x, caption = NULL, label = NULL, align = NULL,
                                digits = NULL, display = NULL, auto = FALSE,
                                ...) {
  x <- data.frame(x$Coef, check.names = FALSE)

  class(x) <- c("xtable","data.frame")
  caption(x) <- caption
  label(x) <- label
  if(auto && is.null(align))   align   <- xalign(x)
  if(auto && is.null(digits))  digits  <- xdigits(x)
  if(auto && is.null(display)) display <- xdisplay(x)
  align(x) <- switch(1+is.null(align), align, c("r","r","r","r","r"))
  digits(x) <- switch(1+is.null(digits), digits, c(0,4,4,2,4))
  display(x) <- switch(1+is.null(display), display, c("s","f","f","f","f"))
  return(x)
}

### spautolm objects: added by David Scott, 6/1/2016, after suggestion by
### Guido Schulz
### Date: Wed, 29 Apr 2015 10:45:16 +0200
### Guido Schulz <schulzgu@student.hu-berlin.de>
xtable.spautolm <- function(x, caption = NULL, label = NULL, align = NULL,
                            digits = NULL, display = NULL, auto = FALSE, ...) {
    return(xtable.summary.sarlm(summary(x), caption = caption, label = label,
                              align = align, digits = digits,
                              display = display, auto = auto, ...))
}

xtable.summary.spautolm <- function(x, caption = NULL, label = NULL,
                                    align = NULL, digits = NULL,
                                    display = NULL, auto = FALSE, ...) {
    return(xtable.summary.sarlm(summary(x), caption = caption, label = label,
                              align = align, digits = digits,
                              display = display, auto = auto, ...))
}


### gmsar objects
xtable.gmsar <- function(x, caption = NULL, label = NULL, align = NULL,
                         digits = NULL, display = NULL, auto = FALSE, ...) {
    return(xtable.summary.sarlm(summary(x), caption = caption, label = label,
                              align = align, digits = digits,
                              display = display, auto = auto, ...))
}

xtable.summary.gmsar <- function(x, caption = NULL, label = NULL, align = NULL,
                                 digits = NULL, display = NULL,
                                 auto = FALSE, ...) {
  return(xtable.summary.sarlm(x, caption = caption, label = label,
                              align = align, digits = digits,
                              display = display, auto = auto, ...))
}

### stsls objects
xtable.stsls <- function(x, caption = NULL, label = NULL, align = NULL,
                         digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.summary.sarlm(summary(x), caption = caption, label = label,
                              align = align, digits = digits,
                              display = display, auto = auto, ...))
}

xtable.summary.stsls <- function(x, caption = NULL, label = NULL, align = NULL,
                                 digits = NULL, display = NULL,
                                 auto = FALSE, ...) {
  return(xtable.summary.sarlm(x, caption = caption, label = label,
                              align = align, digits = digits,
                              display = display, auto = auto, ...))
}

### pred.sarlm objects
xtable.sarlm.pred <- function(x, caption = NULL, label = NULL, align = NULL,
                              digits = NULL, display = NULL,
                              auto = FALSE, ...) {
  return(xtable(as.data.frame(x), caption = caption, label = label,
                align = align, digits = digits,
                display = display, auto = auto, ...))
}


### This method removed because of the need to copy code to pass CRAN checks
### lagImpactMat is neither exported nor documented in spdep

###lagImpact objects
xtable.lagImpact <- function(x, caption = NULL, label = NULL, align = NULL,
                             digits = NULL, display = NULL,
                             auto = FALSE, ...) {
  xtable(lagImpactMat(x), caption = caption, label = label,
         align = align, digits = digits,
         display = display, auto = auto, ...)
}

### package splm
### splm objects
xtable.splm <- function(x, caption = NULL, label = NULL, align = NULL,
                        digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.summary.splm(summary(x), caption = caption, label = label,
                             align = align, digits = digits,
                             display = display, auto = auto, ...))
}

xtable.summary.splm <- function(x, caption = NULL, label = NULL, align = NULL,
                                digits = NULL, display = NULL, auto = FALSE,
                                ...) {
  x <- data.frame(x$CoefTable, check.names = FALSE)

  class(x) <- c("xtable","data.frame")
  caption(x) <- caption
  label(x) <- label
  if(auto && is.null(align))   align   <- xalign(x)
  if(auto && is.null(digits))  digits  <- xdigits(x)
  if(auto && is.null(display)) display <- xdisplay(x)
  align(x) <- switch(1+is.null(align), align, c("r","r","r","r","r"))
  digits(x) <- switch(1+is.null(digits), digits, c(0,4,4,2,4))
  display(x) <- switch(1+is.null(display), display, c("s","f","f","f","f"))
  return(x)
}

### package sphet
### sphet objects
xtable.sphet <- function(x, caption = NULL, label = NULL, align = NULL,
                         digits = NULL, display = NULL, auto = FALSE, ...) {
  return(xtable.summary.splm(summary(x), caption = caption, label = label,
                             align = align, digits = digits,
                             display = display, auto = auto, ...))
}

xtable.summary.sphet <- function(x, caption = NULL, label = NULL, align = NULL,
                                 digits = NULL, display = NULL,
                                 auto = FALSE, ...) {
  return(xtable.summary.splm(x, caption = caption, label = label,
                             align = align, digits = digits,
                             display = display, auto = auto, ...))
}
