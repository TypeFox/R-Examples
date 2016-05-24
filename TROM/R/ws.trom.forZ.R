ws.trom.forZ <-
function(gene_lists, N=N){
  
  sp_associated_idx <- gene_lists
  #names <- colnames(sp_gene_expr)[2:ncol(sp_gene_expr)]
  
  #### calculate TROM scores
  sp_overlap<-sapply(1:length(sp_associated_idx), FUN=function(k) 
    sapply(1:length(sp_associated_idx), FUN=function(j) {
      if(j<k){ return(1)}
      else{
        num <- length(intersect(sp_associated_idx[[k]], sp_associated_idx[[j]]))
        m <- length(sp_associated_idx[[k]])
        n <- length(sp_associated_idx[[j]])
        sum(sapply(num:min(m,n), P_mass, m=m, n=n, N=N))
      }
    }))
  sp_overlap_old <- sp_overlap
  sp_overlap[sp_overlap_old==0]=300
  sp_overlap[sp_overlap_old!=0] <- -log(sp_overlap_old[sp_overlap_old!=0]*nrow(sp_overlap_old)*ncol(sp_overlap_old), base=10)
  sp_overlap[sp_overlap<=0] <- 0
  
  return(sp_overlap)
}
