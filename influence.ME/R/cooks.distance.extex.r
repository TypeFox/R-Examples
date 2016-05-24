cooks.distance.estex <- function(model, parameters=0, sort=FALSE, ...) 
	{
 		
		estex <- model
		
		n.groups <- dim(estex$alt.fixed)[1]
		n.parameters <- dim(estex$alt.fixed)[2]
		ifelse(parameters==0, sel <- 1:n.parameters, sel <- parameters)	
		
		a <- as.matrix(estex$or.fixed[,sel])
		b <- as.matrix(estex$alt.fixed[,sel])
		e <- NA
		
		if (n.groups==1)
			{
				d <- as.matrix(estex$alt.vcov[[1]][sel, sel])
				e <- t(a-b) %*% solve(d) %*% (a-b) / length(a)
			}
			
		if(n.groups > 1)
			{
			for (i in 1:n.groups)
				{	
					d <- as.matrix(estex$alt.vcov[[i]][sel, sel])
					e[i] <- t(a-b[i,]) %*% solve(d) %*% (a-b[i,]) / length(a)
				}
			
			e <- as.matrix(e) 
			rownames(e) <- rownames(b)
			}
		
		if(sort == TRUE)
			{
			e <- as.matrix(sort(e[,1]))

			}
			
		return(e)
	}