#################################################################################################
qdg.perm.test <- function(cross,nperm,node1,node2,common.cov=NULL,DG,QTLs,addcov=NULL,intcov=NULL)
{
  ##################################################
  permuta.block.pheno <- function(cross,node1,node2,common.cov,addcov,intcov)
    {
      pheno <- cross$pheno
      le <- length(pheno[,1])
      aux <- sample(c(1:le),le,replace=FALSE)
      perm.pheno <- pheno
      block <- c(node1,node2,common.cov,addcov,intcov)
      perm.pheno[,block] <- pheno[aux,block]
      perm.cross <- cross
      perm.cross$pheno <- perm.pheno
      perm.cross
    }
  #######################################################
  intersect <- function(x, y) y[match(x, y, nomatch = 0)]
  ############################################################
  pair <- intersect(x = which(DG[,1] == node1), y = which(DG[,3] == node2))
  if(!length(pair))
    stop(paste(node1, "and", node2, "do not have a directed edge"))
  
  cov <- get.covariates(pair=pair, DG=DG)
  obs.lod <- DG[pair,4]
  ps <- rep(0,nperm)
  for(i in 1:nperm){
    perm.cross <- permuta.block.pheno(cross,node1,node2,common.cov,addcov,intcov)
    ps[i] <- lod.score(cross=perm.cross, node1=node1, node2=node2, 
                       qtl.node1=QTLs[[node1]], qtl.node2=QTLs[[node2]], 
                       cov.node1=c(cov[[1]],addcov), 
                       cov.node2=c(cov[[2]],addcov),
                       intcov=intcov,artfact.qtl=QTLs[[1]])
  }
  if(obs.lod > 0) pvalue <- length(which(ps >= obs.lod))/nperm
  else pvalue <- length(which(ps <= obs.lod))/nperm	
  mylist <- list(pvalue,obs.lod,ps,node1,node2)
  names(mylist) <- c("pvalue","obs.lod","permSample","node1","node2")
  class(mylist) <- c("qdg.perm.test", "list")

  mylist
}

summary.qdg.perm.test <- function(object, ...)
{
  cat("\nNodes:\n")
  print(c(object$node1,object$node2))
  cat("\nPermutation p-value for direction:\n")
  print(object$pvalue)
  cat("\nObserved direction LOD score:\n")
  print(object$obs.lod)
  invisible()
}

print.qdg.perm.test <- function(x, ...) summary(x, ...)
