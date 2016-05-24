summary.clustconst<-function(object, mode="all", name=NULL, sort=TRUE, minconst=0.5, digits=3,...) {
	x <- as.data.frame(object)
	if(mode=="cluster") {
		clustIndex = which(names(x)==name)
		spnames = row.names(x)
		values = x[,clustIndex]
		if(sort) {
			o = order(x[, clustIndex],decreasing=TRUE)
			spnames = spnames[o]
			values = values[o]
		}
		spIndices = which(values>= minconst)
		for(i in 1:length(spIndices)) {
			cat(paste(spnames[spIndices[i]],format(round(values[spIndices[i]],digits=digits),nsmall=digits),"\n"))
		}
	}
	else if(mode=="species") {
		spIndex = which(row.names(x)==name)
		clnames = names(x)
		values = x[spIndex,]
		if(sort) {
			o = order(x[spIndex,],decreasing=TRUE)
			clnames = clnames[o]
			values = values[o]
		}	
		clIndices = which(values>=minconst)
		for(i in 1:length(clIndices)) {
			cat(paste(clnames[clIndices[i]],format(round(values[clIndices[i]],digits=digits),nsmall=digits),"\n"))
		}
	}
	else if(mode=="all") {
		x =x[apply(x,1,max)>minconst,]
		if(sort) {
			y = sweep(x,1,rowSums(x),"/")
			oc = order(colSums(y), decreasing=TRUE)
			y = y[,oc]
			x = x[,oc]
			o = integer(nrow(x))
			t = 0
			used=logical(nrow(x))
			for(k in 1:ncol(x)) {
				indices = which(apply(y,1,which.max)==k & !used)
				if(length(indices)>0) {
					o[(t+1):(t+length(indices))] = indices[order(x[indices,k], decreasing=TRUE)]
					t = t+length(indices)
					
					cat(paste("------------",names(x)[k],"-------------\n"))
					print(format(round(x[indices[order(x[indices,k], decreasing=TRUE)],],digits=digits), nsmall=digits))
					used[indices] = TRUE
				}
			}
			if(sum(!used)>0) {
				o[(t+1):(t+sum(!used))] = which(!used)
				cat(paste("------------ REMAINING -------------\n"))
				print(format(round(x[which(!used),],digits=digits), nsmall=digits))
			}
			x = x[o,]
		}
		invisible(x)
	}
}