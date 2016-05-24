sp.ortholog.z <-
function(sp1_sp2_orthologs,sp_gene_expr,i){
    sp_ortholog_data<-t( sapply(sp1_sp2_orthologs[,i], FUN=function(x) sp_gene_expr[which(sp_gene_expr[,1]==as.character(x)),]) )
    z<-cbind(as.character(sp_ortholog_data[,1]), t(apply(sp_ortholog_data[,2:ncol(sp_ortholog_data)], 1, FUN=function(x) {
      x <- as.numeric(x)
      (x-mean(x))/sd(x)
    } )) )
    colnames(z) <- colnames(sp_ortholog_data)
    return(z)
  }
