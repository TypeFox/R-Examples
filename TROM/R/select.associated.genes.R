select.associated.genes <-
function(sp_gene_expr,z_thre=1.5, save = TRUE, 
                                  plot_distribution = FALSE){
  sp_associated_idx <- sp.associated.idx(sp_gene_expr,z_thre)
  sp_associated_genes_all <- sapply( sp_associated_idx, FUN=function(x) unique(sp_gene_expr[x,1]) )
  
  if (plot_distribution==TRUE) {
    sp_count_table <- table(unlist(sp_associated_genes_all))
    
      sub.bar <- unique(quantile(sp_count_table))
 
    sp_gene_counts <- sapply(sp_associated_genes_all, FUN=function(x) {
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
    
    pdf("number of sample associated genes.pdf",width=16,height=8)
    legend <- paste("num <=",sub.bar[1])
      for (i in 2:nrow(num)){
        legend <- c(legend , paste(sub.bar[i-1],"< num <=",sub.bar[i]))
      }
    main=paste("Number of associated genes\n", "( z_thre=",z_thre, ")")
    barplot(num, names.arg=colnames(sp_gene_expr)[2:ncol(sp_gene_expr)], las=2,
            cex.names=0.7, 
            main=main,
            font.main=4, 
            legend.text=legend)
    dev.off()
  }
  
  sp_associated_genes_all_raw <- sp_associated_genes_all
  
  if (save == TRUE){
    listlen<-as.vector(unlist(lapply(sp_associated_genes_all,length)))
    mlen<-max(listlen)
    sp_associated_genes_all<-lapply(sp_associated_genes_all, FUN=function(x,mlen){   
      x<-list(c(as.vector(unlist(x)),rep("",mlen-length(x))))
    },mlen=mlen)
    sp_associated_genes_all<- data.frame(matrix(unlist(sp_associated_genes_all), ncol=ncol(sp_gene_expr)-1,
                                                byrow=F))
    colnames(sp_associated_genes_all)=colnames(sp_gene_expr)[2:ncol(sp_gene_expr)]
    write.xlsx(sp_associated_genes_all,"associated genes",colNames=TRUE)
  }
  
  return(sp_associated_genes_all_raw)
}
