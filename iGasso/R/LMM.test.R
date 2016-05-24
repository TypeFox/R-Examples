LMM.test = function(g, y1, K, S0.1)
{
    p = mean(g, na.rm=T)/2
    if (is.na(p) | p<0.01 | p>0.99) stat = NA   # MAF > 0.01, the same criterion kinf.BN
    else {
    	s = !is.na(g) 
    	if (sum(!s) > 0) {
    		if (sum(!s) <2) S0.1 = S0.1[s, s] - sum(S0.1[!s,s]^2)/S0.1[!s,!s]
    		else S0.1 = S0.1[s, s] - crossprod(S0.1[!s,s], solve(S0.1[!s,!s],S0.1[!s,s])) 
    	}
    	gg = (g[s]-2*p)/sqrt(2*p*(1-p))
    	G = outer(gg, gg) - K[s, s]
    	S.y1 = crossprod(S0.1, y1[s])
    	S.G = crossprod(S0.1, G)
    	b = crossprod(S.y1, crossprod(G, S.y1)) - sum(diag(S.G))
    	stat = b^2/(2*sum(S.G * t(S.G)))
    }
    
	names(stat) = "statistic"
    structure(list(statistic = stat, p.value = (1-pchisq(stat, 1))/2,
                   method = paste("Score Statistic for the Linear Mixed Model"), 
                   data.name = deparse(substitute(LMM.test(g, y1, K)))), class = "htest") 
}

