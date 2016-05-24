# Read a pedfile (which includes the variable names in the top line) 
# and build haplotypes using the markers which appear third, second, and 
# first in the pedfile. 

  filespec <- file.path(.path.package("tdthap"),"tests","test.ped")
  ped <- read.table(filespec)

# call hap.transmit, tdt.select, and tdt.rr

  haps <- hap.transmit(ped, markers=c(3,2,1))
  hap.use <- tdt.select(haps, markers=1:2)
  table(hap.use$trans)
  table(hap.use$untrans)
  rr <- tdt.rr(hap.use)
  rr

# To do a Geary_Moran test on a 10 marker haplotype
  gaps <- c(0, 50, 60, 80, 20, 30, 50, 40, 50, 100, 0)
  set.similarity(nloci=10, spacing=gaps, power=0.5)
  get.similarity(nloci=10)

  test <- tdt.quad(hap.use, nsim=10000, keep=T)
  names(test) 
  summary(test$sim)
