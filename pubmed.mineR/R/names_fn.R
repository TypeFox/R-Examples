names_fn =  function(genes, data, abs,filename, terms){indices = NULL; for (i in 1:dim(genes)[1]){indices = c(indices,which(genes[i,1] == data[,2]))  }
data_theme <- data[indices,]
name = NULL;for (i in 1:dim(genes)[1]){ name = c(name, list(unlist(strsplit(as.character(data_theme[i,3]),"|",fixed=T)))) }
result_genes_name = NULL; for (i in 1:dim(genes)[1]){if ( !is.na(name[[i]][1]) ) {for (j in 1:length(name[[i]])) result_genes_name = Give_Sentences(name[[i]][j],abs) ; if (length(result_genes_name) != 0) {print(c(i,j));write(paste(">>",data_theme[i,2],name[[i]][j],sep=" "), file = paste(filename,"names.txt",sep=""), append=T); for (k in 1:length(result_genes_name)){for(l in 1: length(result_genes_name[[k]])) {temp = result_genes_name[[k]][l];  for(s in 1:length(terms)){temp1 = regexpr(terms[s],temp);if (temp1 != -1) write(c(attr(result_genes_name,"PMID")[k],temp), file = paste(filename,"names.txt",sep=""),append=T)}   }}}    }}
}

