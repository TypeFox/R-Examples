sp.associated.idx <-
function(sp_gene_expr,z_thre=1.5){
  sp.std <- function(sp_gene_expr){
    cbind(sp_gene_expr[,1], t(apply(sp_gene_expr[,2:ncol(sp_gene_expr)], 1, FUN=function(x) {
      x <- as.numeric(x)
      (x-mean(x))/sd(x)
    } )) )
  }
  sp_z=sp.std(sp_gene_expr)
  sapply(2:ncol(sp_z), FUN=function(i) 
  which(sp_z[,i]>z_thre & rowSums(sp_gene_expr[,2:ncol(sp_gene_expr)])>0), simplify=FALSE)  
}
