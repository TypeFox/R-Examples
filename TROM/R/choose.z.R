choose.z <-
function(sp_gene_expr){
  sp_gene_expr <- sp_gene_expr[apply(sp_gene_expr[2:ncol(sp_gene_expr)],1,sum)!=0, ]
  N = nrow(sp_gene_expr)
  
  TROM.list <- list()
  zseq <- seq(-2,3,by=0.1)
  
  for (i in 1:length(zseq)){
    z <- zseq[i]
    print(paste("z=",z,sep=""))
    associated_idx <- sp.associated.idx(sp_gene_expr,z_thre=z)
    
    #trueTROM <- ws.trom.bi(associated_idx,N)
    trueTROM <- ws.trom.forZ(associated_idx,N)
    TROM.list[[i]] <- trueTROM
    
  }
  
  
  pat <- unlist(lapply(TROM.list, function(A){
    mean(A[row(A)!=col(A)])
  }))
  logpat <- log(pat+1,base=10)
  dens <- density(logpat,from=min(logpat), to=max(logpat))
  upy <- max(dens$y)
  lowy <- upy*90/100
  lowx <- range(dens$x[dens$y>=lowy])
  z_range_temp <-zseq[logpat >= lowx[1] & logpat <= lowx[2]]
  z_range <- z_range_temp[-(1:which(diff(z_range_temp)>0.2))]
  
  return(max(z_range))
  
}
