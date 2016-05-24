ws.trom.orthologs <-
function(sp1_sp2_orthologs, sp_gene_expr = NULL, 
                            single = TRUE, sp_gene_expr2 = NULL,
                            z_thre = 1.5, i, 
                            provide = FALSE, gene_lists = NULL, save_overlap_genes = FALSE){
  if (single == TRUE){
    if (provide == FALSE){
      sp_ortholog_data<-t( sapply(sp1_sp2_orthologs[,i], FUN=function(x) sp_gene_expr[which(sp_gene_expr[,1]==as.character(x)),]) )
      
      sp_ortholog_z<-sp.ortholog.z(sp1_sp2_orthologs,sp_gene_expr,i)
      
      Ind1 <- apply(sp_ortholog_data[,2:ncol(sp_ortholog_data)], 1, function(r){
        sum(unlist(r))>0
      })
      sp_associated_idx <- sapply(2:ncol(sp_ortholog_z), FUN=function(i) which(sp_ortholog_z[,i]>z_thre & Ind1), simplify=FALSE)
      #    sp_associated_idx <- sapply(2:ncol(sp_ortholog_z), FUN=function(i) which(sp_ortholog_z[,i]>z_thre & rowSums(sp_gene_expr[,2:ncol(sp_gene_expr)])>0), simplify=FALSE)
      sp_genes<-as.character(sp_ortholog_data[,1])
      
      overlap<-sapply(1:length(sp_associated_idx), FUN=function(k) {
        x <- unique(sp_genes[sp_associated_idx[[k]]])
        sapply(1:length(sp_associated_idx), FUN=function(j) {
          y <- unique(sp_genes[sp_associated_idx[[j]]])
          num <- length(intersect(x, y))
          m <- length(x)
          n <- length(y)
          sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=length(unique(sp_genes))))
        })
      })
      rownames(overlap) <- colnames(sp_ortholog_data)[2:ncol(sp_ortholog_data)]
      colnames(overlap) <- rownames(overlap)
      overlap_old <- overlap
      overlap[overlap_old==0]=300
      overlap[overlap_old!=0] <- -log(overlap_old[overlap_old!=0]*nrow(overlap_old)*ncol(overlap_old), base=10)
      overlap[overlap<=0] <- 0
      write.xlsx(overlap, "within-species TROM scores (with ortholog genes).xlsx", row.names=TRUE, colNames=TRUE)
      
    }else{
      sp_interested_genes <- read.xlsx(gene_lists, colNames=TRUE)
      names<-colnames(sp_interested_genes)
      sp_interested_genes <- lapply(1:length(sp_interested_genes), FUN=function(x) 
        sp_interested_genes[[x]][!is.na(sp_interested_genes[[x]])])
      sp_interested_genes <- lapply(1:length(sp_interested_genes), FUN=function(x) 
        sp_interested_genes[[x]][!(sp_interested_genes[[x]]=="")])
      
      sp_interested_ortholog <- sapply(1:length(sp_interested_genes), FUN=function(x)
        sp_interested_genes[[x]][which(sp_interested_genes[[x]] %in% sp1_sp2_orthologs[,1])])
      
      overlap<-sapply(sp_interested_ortholog, FUN=function(x) sapply(sp_interested_ortholog, FUN=function(y) {
        num <- length(intersect(x, y))
        m <- length(x)
        n <- length(y)
        N <- length(unique(unlist(sp_interested_genes)))
        sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=N))
      }))
      rownames(overlap) <- names
      colnames(overlap) <- names
      overlap_old <- overlap
      overlap[overlap_old==0]=300
      overlap[overlap_old!=0] <- -log(overlap_old[overlap_old!=0]*nrow(overlap_old)*ncol(overlap_old), base=10)
      overlap[overlap<=0] <- 0
      write.xlsx(overlap, "within-species TROM scores (with ortholog genes).xlsx", row.names=TRUE, colNames=TRUE)
      
      ########## output overlapping genes #######
      if (save_overlap_genes == TRUE){
        if (provide == TRUE){
          sp_associated_idx <- sp_interested_ortholog  
        }
        overlap_genes <- sapply(1:length(sp_associated_idx), FUN=function(k) 
          sapply(1:length(sp_associated_idx), FUN=function(j) {
            temp <- intersect(sp_associated_idx[[k]], sp_associated_idx[[j]])
            if (provide == FALSE){
              temp <- sp_genes[temp]
            }
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
        
        write.xlsx(overlap_genes, "within-species overlapping genes (within ortholog genes) between sample pairs", colNames=TRUE)
      }
      ############################## end output
    }
  }else{
    overlap <- ws.trom.three(sp1_gene_expr=sp_gene_expr, sp2_gene_expr=sp_gene_expr2, 
                             sp1_sp2_orthologs, i, z_thre=z_thre,
                             provide=provide, gene_lists=gene_lists, save_overlap_genes = save_overlap_genes)
  }
  
  
  return(overlap)
}
