select.associated.orthologs <-
function(sp_gene_expr, sp1_sp2_orthologs, z_thre=1.5,  i, save = TRUE,
                                        plot_distribution = FALSE){
  
  sp_ortholog_data <- t( sapply(sp1_sp2_orthologs[,i], FUN=function(x) sp_gene_expr[which(sp_gene_expr[,1]==as.character(x)),]) )
  
  sp_genes <- as.character(sp_ortholog_data[,1])
  
  sp_ortholog_z<-sp.ortholog.z(sp1_sp2_orthologs,sp_gene_expr,1)
  Ind1 <- apply(sp_ortholog_data[,2:ncol(sp_ortholog_data)], 1, function(r){
    sum(unlist(r))>0
  })
  sp_specific_idx <- sapply(2:ncol(sp_ortholog_z), FUN=function(i) which(sp_ortholog_z[,i]>z_thre & Ind1), simplify=FALSE)
  
  sp_specific_genes_w_orth <- sapply( sp_specific_idx, FUN=function(x) unique(sp_genes[x]) )
  
  if (plot_distribution==TRUE) {
    sp_count_table <- table(unlist(sp_specific_genes_w_orth))
    
    sub.bar <- unique(quantile(sp_count_table))
    sp_gene_counts <- sapply(sp_specific_genes_w_orth, FUN=function(x) {
      return(sp_count_table[x])
    })
    num = matrix(nrow=length(sub.bar), ncol=length(sp_gene_counts))
    for (i in 1:nrow(num)){
      if (i==1){
        num[i,] <- sapply(sp_gene_counts, FUN=function(x) sum(x <= sub.bar[i]))
      }else {
        num[i,] <- sapply(sp_gene_counts, FUN=function(x) sum(x > sub.bar[i-1] & x<=sub.bar[i]))
      }
    }
    
    pdf("number of sample associated orthologous genes.pdf",width=15,height=8)
      legend <- paste("num <=",sub.bar[1])
      for (i in 2:nrow(num)){
        legend <- c(legend, paste(sub.bar[i-1],"< num <=",sub.bar[i]))
      }
      main=paste("Number of associated orthologous genes\n", "( z_thre=",z_thre,")")
    barplot(num, names.arg=colnames(sp_gene_expr)[2:ncol(sp_gene_expr)], las=2,
            cex.names=0.7, main=main, font.main=4, 
            legend.text=legend)
    dev.off()
  }
  
  sp_specific_genes_w_orth_raw <- sp_specific_genes_w_orth
  
  if (save == TRUE){
    listlen<-as.vector(unlist(lapply(sp_specific_genes_w_orth,length)))
    mlen<-max(listlen)
    sp_specific_genes_w_orth<-lapply(sp_specific_genes_w_orth, FUN=function(x,mlen){
      x<-list(c(as.vector(unlist(x)),rep("",mlen-length(x))))
    },mlen=mlen)
    sp_specific_genes_w_orth<- data.frame(matrix(unlist(sp_specific_genes_w_orth), ncol=ncol(sp_gene_expr)-1,
                                                 byrow=F))
    colnames(sp_specific_genes_w_orth)=colnames(sp_gene_expr)[2:ncol(sp_gene_expr)]
    write.xlsx(sp_specific_genes_w_orth,"associated genes within ortholog genes",colNames=TRUE)
  }
  
  return(sp_specific_genes_w_orth_raw)
}
