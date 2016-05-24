"marginal.distr" <- function(otree, with.errors=TRUE,
       edge.weights=if (with.errors) "estimated" else "observed"){

	edge.weights <- match.arg(edge.weights, c("observed","estimated"))
	if (with.errors & is.null(otree$eps)) stop("Need false positive and negative rates")
	
	add.children.probs <- function(event, event.prob){
		children <- which(otree$parent$parent.num==event)
		for (ch in children){
			res[ch] <<- event.prob * ifelse(edge.weights=="estimated", otree$parent$est.weight[ch],
			                                              otree$parent$obs.weight[ch])
			add.children.probs(ch, res[ch])
		}
	}
	
	res <- numeric(otree$nmut)
	res[1] <- 1
	add.children.probs(1,1)
	names(res) <- otree$parent$child
	if (with.errors)  res <- otree$eps["epos"]+(1-otree$eps["epos"]-otree$eps["eneg"])*res
	res["Root"] <- 1
	res	
}

