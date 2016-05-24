
## package: haplo.stats
## test script: haplo.em

## settings
verbose=TRUE

require(haplo.stats)
options(stringsAsFactors=FALSE)
tmp <- Sys.setlocale("LC_COLLATE", "C")
tmp <- Sys.getlocale("LC_COLLATE")

if(verbose) cat("prepare two datasets, one with char alleles, the other 3 loci from hla data\n")
  
  # make ficticious data set with an intention of some trend in
  # haplotypes having H-allele at locus-H with F-allele at locus-F
  geno.char <- matrix(c('F','f','g','G','h1','h1',
                     'F','F','g','G','H','h1',
                     'F','f','g','G','h2','h2',
                     'f','f','g','G','h2','h1',
                     'F','F','g','G','H','h2',
                     'f','f','G','G','h1','h2',
                     'F','f','G','g','h2','h2',
                     'F','F','g','G','h1','z',
                     'F','f','z','z','h1','h1',
                     'F','f','G','g','h1','h2',
                     'F','f','G','G','h1','h2',
                     'F','F','g','G','h1','z',
                     'F','f','z','z','h1','h1',
                     'f','f','G','g','h1','h2'), nrow=14,byrow=T)
char.label <- c("F","G","H")

data(hla.demo)
  
hla.sub <- hla.demo[,c(1,2,3,4,17,18,21:24)]
geno.hla <- hla.sub[,-c(1:4)]
hla.label=c("DQB","DRB","HLA.B")

seed <- c(33, 10, 39,  6, 16,  0, 40, 24, 12, 60,  7,  1)

if(verbose) cat("character alleles\n")

geno.char.recode <- setupGeno(geno.char, miss.val="z")

geno.char.recode

set.seed(seed)
em.char <- haplo.em(geno.char, miss.val='z',locus.label=char.label,
                    control = haplo.em.control())

print.haplo.em(em.char, digits=3)
summary(em.char, digits=2, show.haplo=TRUE)

summary.haplo.em(em.char, digits=2)

  
if(verbose) cat("hla data, 3 loci\n")
set.seed(seed)
em.hla3 <- haplo.em(geno.hla, locus.label=hla.label, miss.val=0,
                    control = haplo.em.control())

print.haplo.em(em.hla3, digits=3)
 
if(verbose) cat("snap SNP data, unphased\n")
snapDF <- read.table("snapData.csv",header=TRUE, sep=",", stringsAsFactors=FALSE)

geno.snap <- setupGeno(geno=snapDF[,-1])
set.seed(seed)
em.snap <- haplo.em(geno=geno.snap)

print(em.snap, digits=3)


if(verbose) cat("Check Phase against SNaP data, phased\n")
snapfile <- "snap.sim.phased.dat"
source("snapFUN.s")
set.seed(seed)
block1.phase <- checkPhase(snapfile, blocknum=1)
set.seed(seed)
block2.phase <- checkPhase(snapfile, blocknum=2)


print(block1.phase[[2]], digits=4)
print(block1.phase[[1]][1:50,],digits=4)  # long output 1120 lines
print(block2.phase[[2]],digits=4)




