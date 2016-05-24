`mu.score` <- 
function(x, y=NULL, frml=NULL, dsgn=1, wght=FALSE) 
	mu.Sums(mu.AND(mu.GE(x, y), frml), dsgn=dsgn, wght=wght)$score

