# describe terms in a manova object by index of unique error terms, size and type

termInfo <- function(x) {
	if (!inherits(x,"Anova.mlm")) stop("Only for Anova.mlm objects")
  terms <- x$terms
	n.terms <- length(terms)
	result <- matrix(0, n.terms, 3)
	SSPE <- x$SSPE

  if (!x$repeated) {       # wholly Between-S design
		result[,1] <-1
		result[,2] <- nrow(SSPE)
		result[,3] <- 0
  	}

	else {                   # some Within-S terms
		P <- x$P
		unique.P <- unique(P)
		for (term in 1:n.terms) {
			for (j in 1:length(unique.P)) {
				if (identical(P[[term]], unique.P[[j]])) {
					result[term,1] <- j
					result[term,2] <- nrow(SSPE[[term]])
					result[term,3] <- if (all(P[[term]]==1)) 0 else 1
					break
					}
				}
			}
		}
		rownames(result) <- terms
		colnames(result) <- c("index", "size", "within")
		as.data.frame(result)
}
