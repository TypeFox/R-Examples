ws.trom.three <-
function(sp1_gene_expr=NULL, sp2_gene_expr=NULL, 
                          sp1_sp2_orthologs, i, z_thre=1.5,
                          provide=FALSE, gene_lists=NULL, save_overlap_genes = FALSE){
  if (provide==FALSE){
    
    sp1_associated_idx <- sp.associated.idx(sp1_gene_expr, z_thre=z_thre)
    ind1 <- which(sp1_gene_expr[,1] %in% sp1_sp2_orthologs[,i])
    sp1_associated_idx <- lapply(sp1_associated_idx, function(list){
      intersect(list, ind1)
    })
    
    sp2_associated_idx <- sp.associated.idx(sp2_gene_expr, z_thre=z_thre)
    ind2 <- which(sp2_gene_expr[,1] %in% sp1_sp2_orthologs[,i])
    sp2_associated_idx <- lapply(sp2_associated_idx, function(list){
      intersect(list, ind2)
    })
    #    sp2_associated_idx <- sp.associated.idx(sp2_ortholog_data, z_thre=z_thre)
    
    #    sp1_ortholog_data <- sp1_gene_expr[sp1_gene_expr[,1] %in% sp1_sp2_orthologs[,i], ]
    #    sp2_ortholog_data <- sp2_gene_expr[sp2_gene_expr[,1] %in% sp1_sp2_orthologs[,i], ]
    
    #    sp1_associated_idx <- sp.associated.idx(sp1_ortholog_data, z_thre=z_thre)
    #    sp2_associated_idx <- sp.associated.idx(sp2_ortholog_data, z_thre=z_thre)
    
    sp1_genes <- as.character(sp1_gene_expr[,1])
    sp2_genes <- as.character(sp1_gene_expr[,1])
    
    overlap<-sapply(sp1_associated_idx, FUN=function(x) sapply(sp2_associated_idx, FUN=function(y) {
      num <- length(intersect(sp1_genes[x], sp2_genes[y]))
      m <- length(x)
      n <- length(y)
      sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=length(unique(sp1_sp2_orthologs[,i])) ))
    }))
    overlap_old <- overlap
    overlap[overlap_old==0]=300
    overlap[overlap_old!=0] <- -log(overlap_old[overlap_old!=0]*nrow(overlap_old)*ncol(overlap_old), base=10)
    overlap[overlap<=0] <- 0
    rownames(overlap) <- colnames(sp2_gene_expr)[2:ncol(sp2_gene_expr)]
    colnames(overlap) <- colnames(sp1_gene_expr)[2:ncol(sp1_gene_expr)]
    write.xlsx(overlap, "within-species TROM scores (with ortholog genes).xlsx", row.names=TRUE, colNames=TRUE)
  }
  
  else{
    sp1_interested_genes <- read.xlsx(gene_lists, sheet=1, colNames=TRUE)
    names1<-colnames(sp1_interested_genes)
    sp1_interested_genes <- lapply(1:length(sp1_interested_genes), FUN=function(x) 
      sp1_interested_genes[[x]][!is.na(sp1_interested_genes[[x]])])
    sp1_interested_genes <- lapply(1:length(sp1_interested_genes), FUN=function(x) 
      sp1_interested_genes[[x]][!(sp1_interested_genes[[x]]=="")])
    
    sp2_interested_genes <- read.xlsx(gene_lists, sheet=2, colNames=TRUE)
    names2<-colnames(sp2_interested_genes)
    sp2_interested_genes <- lapply(1:length(sp2_interested_genes), FUN=function(x) 
      sp2_interested_genes[[x]][!is.na(sp2_interested_genes[[x]])])
    sp2_interested_genes <- lapply(1:length(sp2_interested_genes), FUN=function(x) 
      sp2_interested_genes[[x]][!(sp2_interested_genes[[x]]=="")])
    
    sp1_interested_ortholog <- sapply(1:length(sp1_interested_genes), FUN=function(x)
      sp1_interested_genes[[x]][which(sp1_interested_genes[[x]] %in% sp1_sp2_orthologs[,i])])
    sp2_interested_ortholog <- sapply(1:length(sp2_interested_genes), FUN=function(x)
      sp2_interested_genes[[x]][which(sp2_interested_genes[[x]] %in% sp1_sp2_orthologs[,i])])
    
    overlap<-sapply(sp1_interested_ortholog, FUN=function(x) sapply(sp2_interested_ortholog, FUN=function(y) {
      num <- length(intersect(x, y))
      m <- length(x)
      n <- length(y)
      sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=length(unique(sp1_sp2_orthologs[,i])) ))
    }))
    overlap_old <- overlap
    overlap[overlap_old==0]=300
    overlap[overlap_old!=0] <- -log(overlap_old[overlap_old!=0]*nrow(overlap_old)*ncol(overlap_old), base=10)
    overlap[overlap<=0] <- 0
    rownames(overlap) <- names2
    colnames(overlap) <- names1
    write.xlsx(overlap, "within-species TROM scores (with ortholog genes).xlsx", row.names=TRUE, colNames=TRUE)
  }
  
  if (save_overlap_genes == TRUE){
    ################## output overlap genes of species 1
    if (provide == TRUE){
      sp1_associated_idx <- sp1_interested_ortholog
      sp2_associated_idx <- sp2_interested_ortholog
    }
    
    overlap_genes <- sapply(1:length(sp1_associated_idx), FUN=function(k) {
      sapply(1:length(sp2_associated_idx), FUN=function(j){
        if (provide == TRUE){
          temp1 <- intersect(sp1_associated_idx[[k]], sp2_associated_idx[[j]])
        }else{
          temp1 <- intersect(sp1_genes[sp1_associated_idx[[k]] ], sp2_genes[sp2_associated_idx[[j]] ])
        }
        if (length(temp1)==0){temp1 <- ""}
        names(temp1) <- paste(k,"_",j)
        temp1
      })
    }, simplify = FALSE)
    
    overlap_genes <- unlist(overlap_genes, recursive = FALSE)
    names(overlap_genes) <- sapply(1:length(overlap_genes), FUN=function(x)
      names(overlap_genes[[x]][1]))
    
    listlen<-as.vector(unlist(lapply(overlap_genes,length)))
    mlen<-max(listlen)
    overlap_genes <- lapply(overlap_genes, FUN=function(x,mlen){   
      x<-list(c(as.vector(unlist(x)),rep("",mlen-length(x))))
    },mlen=mlen)
    
    overlap_genes_names1 <- names( overlap_genes )
    overlap_genes <- unlist(overlap_genes, recursive = FALSE)
    overlap_genes <- data.frame(matrix(unlist(overlap_genes), ncol = length(sp1_associated_idx)*length(sp2_associated_idx), byrow=FALSE), 
                                stringsAsFactors=FALSE)
    colnames(overlap_genes) <- overlap_genes_names1
    
    write.xlsx(overlap_genes, "within-species overlapping genes (within ortholog genes) between sample pairs.xlsx", colNames=TRUE)
    
    ############################## end output
  }
  return(t(overlap))
}
