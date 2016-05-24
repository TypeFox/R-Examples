#8-9-2014 MRC-Epid JHZ

ReadGRMPLINK <- function(prefix, diag=1)
{
  fn <- paste(prefix,".genome",sep="")
  genome <- read.table(fn,header=TRUE)
  require(bdsmatrix)
  genome <- with(genome, {
     ID <- unique(c(FID1,FID2))
     t1 <- cbind(FID1,FID2,PI_HAT)
     t2 <- cbind(FID1=ID,FID2=ID,PI_HAT=diag)
     trio <- rbind(t1,t2)
     bdsmatrix.ibd(trio)
  })
}
# adapted from cross.R of the NSHD analysis
# IDs have to be integer, so the version in R/gap NAMESPACE is more general
