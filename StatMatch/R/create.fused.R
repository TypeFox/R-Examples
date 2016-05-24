`create.fused` <-
function (data.rec, data.don, mtc.ids, z.vars, dup.x=FALSE, match.vars=NULL) 
{
	r.id <- mtc.ids[,1]
	d.id <- mtc.ids[,2]
	A <- data.rec[r.id,]
	if(dup.x) {
		B <- data.frame(data.don[d.id, c(match.vars, z.vars)])
		lab  <- c(paste(match.vars, "don", sep="."), z.vars)
		names(B) <- lab
	}	
	else {
		B <- data.frame(data.don[d.id, z.vars])
		names(B) <- z.vars
	}	
	AB <- cbind(A, B) 
	row.names(AB) <- row.names(A) 
	AB
}

