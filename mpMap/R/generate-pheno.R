generate_pheno <- function(n.founders, qtlgeno, qtleffects, vare, n.ind)
{
  if (is.null(qtleffects)) n.qtl <- 0 else n.qtl <- nrow(qtleffects)
 
  pheno <- rnorm(n.ind, 0, vare)
  
  if (n.qtl>0)
  for (i in 1:n.qtl)
  for (j in 1:n.founders)
  {
    pheno[qtlgeno[,i]==j] <- pheno[qtlgeno[,i]==j]+qtleffects[i,j+2]
    pheno[qtlgeno[,i+n.qtl]==j] <- pheno[qtlgeno[,i+n.qtl]==j]+qtleffects[i,j+2]  }

  return(pheno)
}

