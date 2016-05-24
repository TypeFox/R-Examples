toTable <-
function(mGSZobj, sample = FALSE, m = c("mGSZ","mGSA","mAllez","WRS","SS","SUM","KS","wKS"), n = 5){

  if(!mGSZobj$gene.perm.log & !sample & !mGSZobj$other.methods){
    table <- mGSZobj$mGSZ[1:n,]
  }
  
  else if(mGSZobj$gene.perm.log & !sample & !mGSZobj$other.methods){
    table <- mGSZobj$mGSZ.gene.perm[1:n,]
  }
  
  else if(mGSZobj$gene.perm.log & sample & !mGSZobj$other.methods){
    table <- mGSZobj$mGSZ.sample.perm[1:n,]
  }
  
  else if(!mGSZobj$gene.perm.log & !sample & mGSZobj$other.methods){
    m <- match.arg(m)
    table <- mGSZobj[[m]][1:n,]
  }
  
  else if(mGSZobj$gene.perm.log & !sample & mGSZobj$other.methods){
    m <- match.arg(m)
    table <- mGSZobj$gene.perm[[m]][1:n,]
  }
  
  else if(mGSZobj$gene.perm.log & sample & mGSZobj$other.methods){
    m <- match.arg(m)
    table <- mGSZobj$sample.perm[[m]][1:n,]
  }
  return(table)
}
