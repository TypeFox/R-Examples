alias_fn =  function(genes, data, abs,filename, terms)  {indices = NULL; for (i in 1:dim(genes)[1]){indices = c(indices,which(genes[i,1] == data[,2]))  }
data_theme <- data[indices,]
alias = NULL;for (i in 1:dim(genes)[1]){ alias = c(alias, list(unlist(strsplit(as.character(data_theme[i,9]),"|",fixed=T)))) }
result_genes_alias = NULL; for (i in 1:dim(genes)[1]){if ( !is.na(alias[[i]][1]) ) {for (j in 1:length(alias[[i]])) result_genes_alias = Give_Sentences(alias[[i]][j],abs) ; if (length(result_genes_alias) != 0) {print(c(i,j));write(paste(">>",data_theme[i,2],alias[[i]][j],sep=" "), file = paste(filename,"alias.txt",sep=""), append=T); for (k in 1:length(result_genes_alias)){for(l in 1: length(result_genes_alias[[k]])) {temp = result_genes_alias[[k]][l];  for(s in 1:length(terms)){temp1 = regexpr(terms[s],temp);if (temp1 != -1) write(c(attr(result_genes_alias,"PMID")[k],temp), file = paste(filename,"alias.txt",sep=""),append=T)}   }}}    }}
}
