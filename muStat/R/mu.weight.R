`mu.weight` <- 
function(x, y=NULL, frml=NULL, dsgn=1, wght=TRUE) 
	mu.Sums(mu.AND(mu.GE(x, y), frml), dsgn=dsgn, wght=wght)$weight

