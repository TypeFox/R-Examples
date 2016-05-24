plot.estex <- function(x, which="dfbetas",
						sort=FALSE, to.sort=NA, abs=FALSE, cutoff=0,
						parameters=seq_len(ncol(estex$alt.fixed)),
                        groups=seq_len(nrow(estex$alt.fixed)), ...)
{
	estex <- x
	
	# Define a panel function only to be used in this plot.estex function
	estex.dotplot <- function(x,y,...,subscripts,col,pch) 
	{
        panel.dotplot(x,y,col=col[subscripts],pch=pch[subscripts], ...)
    }
	
	# Plot DFBETAS
	if(which=="dfbetas")
	{
		
		plot.matrix <- as.matrix(dfbetas(estex, sort=sort, to.sort=to.sort, abs=abs)[groups, parameters])
		
        if (!cutoff) return(dotplot(plot.matrix, groups=FALSE, ...))
    
    	col <- as.vector(ifelse(abs(plot.matrix) >= abs(cutoff), "red", "blue"))
    	pch <- as.vector(ifelse(abs(plot.matrix) >= abs(cutoff), 17, 19))
     
	    print(dotplot(plot.matrix, col=col, pch=pch, panel=estex.dotplot, groups=FALSE, ...))
        
	}
		
	# Plot Cooks' Distance	
	if(which=="cook")
	{
		plot.matrix <- cooks.distance(estex, parameters=parameters, sort=sort)[groups,]
    
    	if (!cutoff) return(dotplot(plot.matrix, ...))
    
	    col <- rep("blue", length(plot.matrix))
    	col[which(abs(plot.matrix) >= cutoff)] <- "red"
    
	    pch <- rep(19, length(plot.matrix))
    	pch[which(abs(plot.matrix) >= cutoff)] <- 17
    
	    print(dotplot(plot.matrix, col=col, pch=pch, ...))
	}	
		
	# Plot Percentage Change	
	if(which=="pchange")
	{
		
		plot.matrix <- as.matrix(pchange(estex, sort=sort, to.sort=to.sort, abs=abs)[groups, parameters])
		
        if (!cutoff) return(dotplot(plot.matrix, groups=FALSE, ...))
    
		col <- as.vector(ifelse(abs(plot.matrix) >= abs(cutoff), "red", "blue"))
    	pch <- as.vector(ifelse(abs(plot.matrix) >= abs(cutoff), 17, 19))
    
	    print(dotplot(plot.matrix, col=col, pch=pch, panel=estex.dotplot, groups=FALSE, ...))
	}
	
	# Plot Test of Changing Significance	
	if(which=="sigtest")
	{
		
		sigs <- sigtest(estex, test=cutoff, sort=sort, to.sort=to.sort, ...)
		test.stats 		<- lapply(sigs, function(x) return(x[,1]))
		test.changed 	<- lapply(sigs, function(x) return(x[,3]))
		
		plot.matrix 			<- matrix(data=unlist(test.stats), ncol=length(test.stats))
		colnames(plot.matrix) 	<- names(test.stats)
		rownames(plot.matrix)	<- rownames(estex$alt.fixed)
		changed.matrix 			<- matrix(data=unlist(test.changed), ncol=length(test.changed))
		
		col <- as.vector(ifelse(changed.matrix, "red", "blue"))
    	pch <- as.vector(ifelse(changed.matrix, 17, 19))
		
		print(dotplot(plot.matrix, col=col, pch=pch, panel=estex.dotplot, groups=FALSE, ...))
	}	
	
}