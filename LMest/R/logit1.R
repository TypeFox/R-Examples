logit1 <- function(p,ref=1){
	
# compute logit with respect to category ref
	lp = log(p[-ref]/p[ref])
	out = list(lp=lp)
		
}