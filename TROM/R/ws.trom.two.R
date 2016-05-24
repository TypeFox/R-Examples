ws.trom.two <-
function(sp1_gene_expr=NULL, sp2_gene_expr=NULL, sp1_sp2_orthologs, z_thre=1.5,
                        provide=FALSE, gene_lists=NULL, save_overlap_genes = FALSE){
  if (provide==FALSE){
    sp1_ortholog_data<-t( sapply(sp1_sp2_orthologs[,1], FUN=function(x) sp1_gene_expr[which(sp1_gene_expr[,1]==as.character(x)),]) )
    sp2_ortholog_data<-t( sapply(sp1_sp2_orthologs[,2], FUN=function(x) sp2_gene_expr[which(sp2_gene_expr[,1]==as.character(x)),]) )
    
    sp1_genes <- as.character(sp1_ortholog_data[,1])
    sp2_genes <- as.character(sp2_ortholog_data[,1])
    
    sp1_ortholog_z <- sp.ortholog.z(sp1_sp2_orthologs,sp1_gene_expr,1)
    sp2_ortholog_z <- sp.ortholog.z(sp1_sp2_orthologs,sp2_gene_expr,2)
    
    #    sp1_associated_idx <- sapply(2:ncol(sp1_ortholog_z), FUN=function(i) which(sp1_ortholog_z[,i]>z_thre ), simplify=FALSE)
    #    sp2_associated_idx <- sapply(2:ncol(sp2_ortholog_z), FUN=function(i) which(sp2_ortholog_z[,i]>z_thre ), simplify=FALSE)
    Ind1 <- apply(sp1_ortholog_data[,2:ncol(sp1_ortholog_data)], 1, function(r){
      sum(unlist(r))>0
    })
    sp1_associated_idx <- sapply(2:ncol(sp1_ortholog_z), FUN=function(i) which(sp1_ortholog_z[,i]>z_thre & Ind1), simplify=FALSE)
    Ind2 <- apply(sp2_ortholog_data[,2:ncol(sp2_ortholog_data)], 1, function(r){
      sum(unlist(r))>0
    })
    sp2_associated_idx <- sapply(2:ncol(sp2_ortholog_z), FUN=function(i) which(sp2_ortholog_z[,i]>z_thre & Ind2), simplify=FALSE)
    
    overlap<-sapply(sp1_associated_idx, FUN=function(x) sapply(sp2_associated_idx, FUN=function(y) {
      num <- length(which(sp1_sp2_orthologs[,1] %in% sp1_genes[x] & sp1_sp2_orthologs[,2] %in% sp2_genes[y]))
      m <- length(x)
      n <- length(y)
      sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=nrow(sp1_sp2_orthologs)))
    }))
    overlap_old <- overlap
    overlap[overlap_old==0]=300
    overlap[overlap_old!=0] <- -log(overlap_old[overlap_old!=0]*nrow(overlap_old)*ncol(overlap_old), base=10)
    overlap[overlap<=0] <- 0
    rownames(overlap) <- colnames(sp2_ortholog_data)[2:ncol(sp2_ortholog_data)]
    colnames(overlap) <- colnames(sp1_ortholog_data)[2:ncol(sp1_ortholog_data)]
    write.xlsx(overlap, "within-species TROM scores.xlsx", row.names=TRUE, colNames=TRUE)
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
      sp1_interested_genes[[x]][which(sp1_interested_genes[[x]] %in% sp1_sp2_orthologs[,1])])
    sp2_interested_ortholog <- sapply(1:length(sp2_interested_genes), FUN=function(x)
      sp2_interested_genes[[x]][which(sp2_interested_genes[[x]] %in% sp1_sp2_orthologs[,2])])
    
    overlap<-sapply(sp1_interested_ortholog, FUN=function(x) sapply(sp2_interested_ortholog, FUN=function(y) {
      num <- length(which(sp1_sp2_orthologs[,1] %in% x & sp1_sp2_orthologs[,2] %in% y))
      m <- length(x)
      n <- length(y)
      sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=nrow(sp1_sp2_orthologs)))
    }))
    overlap_old <- overlap
    overlap[overlap_old==0]=300
    overlap[overlap_old!=0] <- -log(overlap_old[overlap_old!=0]*nrow(overlap_old)*ncol(overlap_old), base=10)
    overlap[overlap<=0] <- 0
    rownames(overlap) <- names2
    colnames(overlap) <- names1
    write.xlsx(overlap, "within-species TROM scores.xlsx", row.names=TRUE, colNames=TRUE)
  }
  
  if (save_overlap_genes == TRUE){
    ################## output overlap genes of species 1
    if (provide == TRUE){
      sp1_associated_idx <- sp1_interested_ortholog
      #      sp2_associated_idx <- sp2_interested_ortholog
    }
    overlap_genes_sp1 <- sapply(1:length(sp1_associated_idx), FUN=function(k) 
      sapply(1:length(sp2_associated_idx), FUN=function(j) {
        if (provide == TRUE){
          index <- which(sp1_sp2_orthologs[,1] %in% sp1_interested_ortholog[[k]] & 
                           sp1_sp2_orthologs[,2] %in% sp2_interested_ortholog[[j]])
          temp1 <- sp1_sp2_orthologs[,1][index]
        }else {
          index <- which(sp1_sp2_orthologs[,1] %in% sp1_genes[sp1_associated_idx[[k]]] & 
                           sp1_sp2_orthologs[,2] %in% sp2_genes[sp2_associated_idx[[j]]])
          temp1 <- sp1_genes[index]
        }
        if (length(temp1)==0){temp1 <- ""}
        names(temp1) <- paste(k,"_",j)
        temp1
      }), simplify = FALSE)
    
    overlap_genes_sp1 <- unlist(overlap_genes_sp1, recursive = FALSE)
    names(overlap_genes_sp1) <- sapply(1:length(overlap_genes_sp1), FUN=function(x)
      names(overlap_genes_sp1[[x]][1]))
    
    listlen<-as.vector(unlist(lapply(overlap_genes_sp1,length)))
    mlen<-max(listlen)
    overlap_genes_sp1 <- lapply(overlap_genes_sp1, FUN=function(x,mlen){   
      x<-list(c(as.vector(unlist(x)),rep("",mlen-length(x))))
    },mlen=mlen)
    
    overlap_genes_names1 <- names( overlap_genes_sp1 )
    overlap_genes_sp1 <- unlist(overlap_genes_sp1, recursive = FALSE)
    overlap_genes_sp1 <- data.frame(matrix(unlist(overlap_genes_sp1), ncol = length(sp1_associated_idx)*length(sp2_associated_idx), byrow=FALSE), 
                                    stringsAsFactors=FALSE)
    colnames(overlap_genes_sp1) <- overlap_genes_names1
    
    write.xlsx(overlap_genes_sp1, "within-species overlapping genes between sample pairs", colNames=TRUE)
    
    ############################## end output
  }
  return(t(overlap))
}
