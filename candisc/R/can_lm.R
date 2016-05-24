#' Transform a Multivariate Linear model mlm to a Canonical Representation

#' This function uses \code{\link{candisc}} to transform the responses in a
#' multivariate linear model to scores on canonical variables for a given term and then uses
#' those scores as responses in a linear (lm) or multivariate linear model (mlm).

#' @param mod A \code{mlm} object
#' @param term One term in that model
#' @param ... Arguments passed to \code{\link{candisc}}
#' @return A \code{lm} object if \code{term} is a rank 1 hypothesis, otherwise a \code{mlm} object

can_lm   <- function(mod, term, ...) {
#  require(candisc)
  if (!(inherits(mod, "mlm"))) stop("model must be a mlm")
  term.names <- gsub(" ", "", labels(terms(mod)))
  which.term <- which(term == term.names)
  if (length(which.term)==0) stop(term, " does not appear in the model")
  
  can <- candisc(mod, term, ...)
#  term <- mod$term                        # term for which candisc was done
  lm.terms <- can$terms                   # terms in original lm
  scores <- can$scores
  resp <- if (can$rank==1) "Can1" else
          paste("cbind(", paste("Can",1:can$rank,sep="", collapse = ","), ")")
  txt <- paste("lm(", resp, " ~ ",
              paste( lm.terms, collapse = "+"), ", data=scores)" )
  can.mod <- eval(parse(text=txt))
  can.mod
}


