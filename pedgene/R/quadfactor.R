
## Functions: quadfactor for getting quadratic factor over all pedigrees
## for both the burden and kernel tests on autosomes and X-chrom
## Author: Dan Schaid 5/2013, revised: Jason Sinnwell 7/2013

quadfactor <- function(kinmat, chrom, resid, sex, male.dose) {

  if(chrom=="X") {
    sex.factor <- 2*(sex==2) %o% (sex==2) +
        male.dose^2*(sex==1) %o% (sex==1) + 
        sqrt(2)*male.dose*((sex==2) %o% (sex==1)+(sex==1) %o% (sex==2)) 

  } else {
    sex.factor <- 2
  }

  genet.var <- 2*kinmat
  ## if no inbreeding, d = 1 so genet.var = genet.cor, but
  ## below allows for inbred subjects

  d <- sqrt(diag(genet.var))
  genet.cor <- sex.factor * (genet.var / (d%o%d))
  c.factor <-  t(resid) %*% genet.cor %*% resid 
  
  return(c.factor)
}
