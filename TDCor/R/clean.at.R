clean.at <-
function(dataset,l_genes)

{rownames(dataset)[as.vector(na.omit(match(l_genes,rownames(dataset))))]}
