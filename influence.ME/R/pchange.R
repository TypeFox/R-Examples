pchange <- 
function(estex, parameters=0, sort=FALSE, to.sort=NA, abs=FALSE)
	{

		n.groups <- dim(estex$alt.fixed)[1]
		n.parameters <- dim(estex$alt.fixed)[2]
		ifelse(parameters==0, sel <- 1:n.parameters, sel <- parameters)	
					
		a <- as.matrix(estex$or.fixed[,sel])
		b <- as.matrix(estex$alt.fixed[,sel])

		e <- t(apply(b,1,function(x) {((x-a) / x) * 100}))
		colnames(e) <- rownames(a)
		
		
		if(abs == TRUE)
			{
			e <- abs(e)
			}
		
		if(sort == TRUE)
			{
			if(dim(e)[2] == 1)
				{
				to.sort <- colnames(e)
				}
	
			if(dim(e)[2] > 1 & is.na(to.sort))
				{
				stop("Indicate on which variable to sort the percentile change: please specify to.sort")
				}
			e <- e[order(e[, to.sort]), ]
			}
		
		return(e)
		
	}
