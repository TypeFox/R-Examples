ws.trom <-
function(sp_gene_expr = data.frame(), single = TRUE, sp_gene_expr2 = NULL,
                  z_thre = 1.5,  
                  provide = FALSE, gene_lists = "", save_overlap_genes = FALSE){
  if (sum(dim(sp_gene_expr))==0 && gene_lists==""){
    print("Error: Please either input sp_gene_expr or gene_lists.")
    return (0)
  }
  
  if (single == TRUE){ # use single expr data
    if(provide==TRUE){
      # here sp_associated_idx is character
      sp_associated_idx <- read.xlsx(gene_lists,colNames=TRUE)
      names <- colnames(sp_associated_idx)  
      sp_associated_idx <- lapply(1:length(sp_associated_idx), FUN=function(x) 
        sp_associated_idx[[x]][!is.na(sp_associated_idx[[x]])])
      sp_associated_idx <- lapply(1:length(sp_associated_idx), FUN=function(x) 
        sp_associated_idx[[x]][!(sp_associated_idx[[x]]=="")])
      N <- length(unique(unlist(sp_associated_idx)))     
    }
    
    else {
      # here sp_associated_idx is number
      sp_associated_idx <- sp.associated.idx(sp_gene_expr,z_thre)
      N<-as.numeric(nrow(sp_gene_expr))
      names <- colnames(sp_gene_expr)[2:ncol(sp_gene_expr)]
    }
    
    #### calculate TROM scores
    sp_overlap<-sapply(1:length(sp_associated_idx), FUN=function(k) 
      sapply(1:length(sp_associated_idx), FUN=function(j) {
        num <- length(intersect(sp_associated_idx[[k]], sp_associated_idx[[j]]))
        m <- length(sp_associated_idx[[k]])
        n <- length(sp_associated_idx[[j]])
        sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=N))
        #m:associated genes in stage k; n:associated genes in stage j; N: totla#of genes
      }))
    rownames(sp_overlap) <- names
    colnames(sp_overlap) <- names
    sp_overlap_old <- sp_overlap
    sp_overlap[sp_overlap_old==0]=300
    sp_overlap[sp_overlap_old!=0] <- -log(sp_overlap_old[sp_overlap_old!=0]*nrow(sp_overlap_old)*ncol(sp_overlap_old), base=10)
    sp_overlap[sp_overlap<=0] <- 0
    write.xlsx(sp_overlap, "within-species TROM scores.xlsx", row.names=TRUE, colNames=TRUE)
    
    ########## output overlapping genes ###############
    if (save_overlap_genes == TRUE){
      overlap_genes <- sapply(1:length(sp_associated_idx), FUN=function(k) 
        sapply(1:length(sp_associated_idx), FUN=function(j) {
          temp <- intersect(sp_associated_idx[[k]], sp_associated_idx[[j]])
          if (provide == FALSE) {temp <- sp_gene_expr[,1][temp]}
          if (length(temp)==0){temp <- ""}
          names(temp) <- paste(k,"_",j)
          temp
        }), simplify = FALSE)
      overlap_genes <- unlist(overlap_genes, recursive = FALSE)
      names(overlap_genes) <- sapply(1:length(overlap_genes), FUN=function(x)
        names(overlap_genes[[x]][1]))
      
      listlen<-as.vector(unlist(lapply(overlap_genes,length)))
      mlen<-max(listlen)
      overlap_genes <- lapply(overlap_genes, FUN=function(x,mlen){   
        x<-list(c(as.vector(unlist(x)),rep("",mlen-length(x))))
      },mlen=mlen)
      
      overlap_genes_names <- names( overlap_genes )
      overlap_genes <- unlist(overlap_genes, recursive = FALSE)
      overlap_genes <- data.frame(matrix(unlist(overlap_genes), ncol = length(sp_associated_idx)^2, byrow=FALSE), 
                                  stringsAsFactors=FALSE)
      colnames(overlap_genes) <- overlap_genes_names
      
      write.xlsx(overlap_genes, "within-species overlapping genes between sample pairs", colNames=TRUE)
    }
  }else{ # use two different expr data
    cup <- union(sp_gene_expr[,1], sp_gene_expr2[,1])
    sp1_sp2_orthologs  <- cbind(cup, cup)
    sp_overlap <- ws.trom.two(sp1_gene_expr=sp_gene_expr, sp2_gene_expr=sp_gene_expr2, 
                              sp1_sp2_orthologs=sp1_sp2_orthologs, z_thre=1.5,
                              provide=provide, gene_lists=gene_lists, 
                              save_overlap_genes = save_overlap_genes)
  }
  
  return(sp_overlap)
}
