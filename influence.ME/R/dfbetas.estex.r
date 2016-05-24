dfbetas.estex <- function(model, parameters=0, sort=FALSE, to.sort=NA, abs=FALSE, ...) 
{
		estex <- model
		
		n.groups <- dim(estex$alt.fixed)[1]
		n.parameters <- dim(estex$alt.fixed)[2]
		ifelse(parameters==0, sel <- 1:n.parameters, sel <- parameters)	
					
		a <- as.matrix(estex$or.fixed[,sel])
		b <- as.matrix(estex$alt.fixed[,sel])
		c <- as.matrix(estex$alt.se[,sel])
		
		e <- matrix(data=NA, nrow=n.groups, ncol=length(sel))
		
		if(n.groups == 1)
			{
			e <- (a-b) / c
			return(e)
			}
		
		if(n.groups > 1)
			{
			for (i in 1:n.groups) 
				{
				e[i,] <- (a-b[i,]) / c[i,]
				}
			dimnames(e) <- dimnames(estex$alt.fixed[, sel])
			}	
		
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
				stop("Indicate on which variable to sort the DFBETAS: please specify to.sort")
				}
			e <- e[order(e[, to.sort]), ]
			}

                e
	}
