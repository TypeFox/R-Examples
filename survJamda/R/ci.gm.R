ci.gm <-
function(x){
	gm1 = mean(log(x), na.rm = TRUE)
	cil = exp(gm1-(1.96*(sd(log(x), na.rm = TRUE)/sqrt(length(x)))))
	ciupp = exp(gm1+(1.96*(sd(log(x), na.rm = TRUE)/sqrt(length(x)))))
	vec = c(round(cil,2), round(ciupp,2))
	return (vec)
}

