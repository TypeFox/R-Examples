genes_BWI = function(currentabs,previousabs,theme, genes){
check_BWI= lapply(genes, function(x){BWI(currentabs, previousabs, paste(" ",x," ",sep=""), theme)})
unlist(check_BWI);
temp = which(check_BWI != 0 );
new_check_BWI = check_BWI[temp];
new_genes = genes[temp]
new_genes_mat = cbind(new_genes,new_check_BWI )
new_genes_mat = as.data.frame(new_genes_mat)
names(new_genes_mat ) = c("Gene", "BWI")
return(new_genes_mat)}
