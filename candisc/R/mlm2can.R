# transform an mlm to an equivalent canonical representation???
# Q: this seems wrong -- should have only term as RHS of the canonical model

mlm2can <- function(mod, term, ...) {
	require(candisc)
	can <- candisc(mod, term, ...)
#	term <- mod$term                        # term for which candisc was done
	lm.terms <- can$terms                   # terms in original lm
	scores <- can$scores
	resp <- if (can$rank==1) "Can1" else
	        paste("cbind(", paste("Can",1:can$rank,sep="", collapse = ","), ")")
	terms <- paste( lm.terms, collapse = "+")
  txt <- paste( "lm(", resp, " ~ ", terms, ", data=scores)" )
  can.mod <- eval(parse(text=txt))
#  can.mod <- lm(as.formula(paste(resp, " ~ ", terms)), data=scores)
	can.mod
}

