"error.rates<-" <-
function(x, value){
	calc.orig.edge.probs <- function(otree, epos, eneg){	
		orig.edge.children <- function(event){
			pj.star <- md[event]
			children <- which(otree$parent$parent.num==event)
			for (ch in children){
				pi.star <- md[ch]
				orig.edge[ch] <<- (otree$parent$obs.weight[ch] * pj.star - epos*(pi.star+pj.star) + epos^2)/
				                 ((pj.star-epos)*(1-eneg-epos))
				orig.edge.children(ch)
			}
		}
		
		md <- colMeans(otree$data)  #observed
		orig.edge <- numeric(otree$nmut)
		orig.edge[1] <- 1  #the root
		orig.edge.children(1);
		orig.edge
	}
  
	value <- as.numeric(value)	
	x$eps <- c(epos=value[1], eneg=value[2])
    oe <- calc.orig.edge.probs(x, value[1], value[2])
    if (any(oe[-1]>=1)) stop("Edge transition probability >= 1")
    x$parent$est.weight <- oe
    x
}

