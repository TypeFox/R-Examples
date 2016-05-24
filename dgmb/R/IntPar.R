IntPar <-
function(rs,rie,modo,bv)
	{
	rs[lower.tri(rs)] <- 0
	rie[lower.tri(rie)] <- 0
	
	nat <- rep(0,length(modo))
	for (i in 1:ncol(rs))
		{
		for (j in 1:ncol(rs))
			{
			if (rs[i,j]==1 & i<j)
				{
				nat[i] <- "ex" 
				nat[j] <- "en"
				}
			} 
		}

	bex <- 0; ben <- 0	
	for (i in 1:(length(modo)))
		{
		ifelse(nat[i]=="ex", bex <- bex + 1, ben <- ben + 1) 
		}

	ind.ex <- 0; ind.en <- 0
	for (i in 1:length(modo))
		{
		ifelse(modo[i]=="F" & nat[i]=="ex", ind.ex <- ind.ex + bv[i], ind.en <- ind.en + bv[i])
		}

	list(nat=nat,bex=bex,ben=ben,ind.ex=ind.ex,ind.en=ind.en,rie=rie)
	}

