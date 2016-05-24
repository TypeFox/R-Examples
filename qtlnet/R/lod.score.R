lod.score <- function(cross, node1, node2, qtl.node1, qtl.node2, cov.node1=NULL, cov.node2=NULL,
                      intcov = NULL, artfact.qtl) {
  nQ1 <- length(qtl.node1$chr)
  nQ2 <- length(qtl.node2$chr)
  cov.node1 <- unique(c(intcov,cov.node1))
  cov.node2 <- unique(c(intcov,cov.node2))
  if(nQ1 == 0 & nQ2 > 0) qtl.node1 <- qtl.node2
  if(nQ1 > 0 & nQ2 == 0) qtl.node2 <- qtl.node1
  if(nQ1 == 0 & nQ2 == 0) qtl.node1 <- qtl.node2 <- artfact.qtl
  
  old.fitqtl <- compareVersion(qtlversion(), "1.08-43") < 0
  if(old.fitqtl)
    myfitqtl <- function(cross, pheno.col, ...)
      fitqtl(cross$pheno[[pheno.col]], ...)
  else
    myfitqtl <- function(cross, pheno.col, ...)
      fitqtl(cross, pheno.col, ...)
  
  node1.col <- find.pheno(cross, pheno=node1)
  node1 <- data.frame(node1 = cross$pheno[,node1.col])

  node2.col <- find.pheno(cross, pheno=node2)
  node2 <- data.frame(node2 = cross$pheno[,node2.col])
  
  if(is.null(cov.node1)) {
    dat1 <- add1 <- NULL
    dat2 <- node2
  }
  else {
    add1 <- cov.node1
    dat1 <- data.frame(cross$pheno[,add1])
    names(dat1) <- add1
    dat2 <- data.frame(node2, dat1)
  }
  add2 <- c("node2", add1)

  ## Would be better to replace hard-coded 10 with name of this entry.
  aux <- s.formula(addcov = add1, intcov = intcov, nQ = nQ1, dat = dat1) 
  out <- myfitqtl(cross, node1.col, qtl.node1, formula = aux[[1]], cov = aux[[2]],
                  dropone = FALSE)$result.full[10]
  aux <- s.formula(addcov = add2, intcov = intcov, nQ = nQ1, dat = dat2)
  out <- out - myfitqtl(cross, node1.col, qtl.node1, formula = aux[[1]], cov = aux[[2]], 
                        dropone = FALSE)$result.full[10]
  
  if(is.null(cov.node2)) {
    dat1 <- add1 <- NULL
    dat2 <- node1
  }
  else {
    add1 <- cov.node2
    dat1 <- data.frame(cross$pheno[,add1])
    names(dat1) <- add1
    dat2 <- data.frame(node1, dat1)
  }
  add2 <- c("node1", add1)

  aux <- s.formula(addcov = add1, intcov = intcov, nQ = nQ2, dat = dat1)
  out <- out - myfitqtl(cross, node2.col, qtl.node2, formula = aux[[1]], cov = aux[[2]], 
                        dropone = FALSE)$result.full[10]

  aux <- s.formula(addcov = add2, intcov = intcov, nQ = nQ2, dat = dat2)
  out <- out + myfitqtl(cross, node2.col, qtl.node2, formula = aux[[1]], cov = aux[[2]], 
                        dropone = FALSE)$result.full[10]
  out
}
