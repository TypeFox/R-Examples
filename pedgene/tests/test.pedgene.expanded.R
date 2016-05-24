
## Expanded test script for pedgene package
devel=FALSE
if(devel) {
  rfiles <- list.files(path="../R/", pattern="*.R$")
  for(f in rfiles) source(paste("../R/",f,sep=""))
  dfiles <- list.files(path="../data/", pattern="*.RData$")
  for(d in dfiles) load(paste("../data/",d,sep=""))
  library(kinship2)
  library(survey)
} else {
  
  require(pedgene)
  data(example.ped)
  data(example.geno)
  data(example.map)
  
}

#require(survery)


######################################
## From Dan Weeks, issues to check
######################################
## 1) missid ="0" when the rest of the ids are character
## 2) skip pedigree checking, checkpeds=TRUE/FALSE
## 3) character alleles
## 4) disconnected pedigrees
## 5) "flipped" 0/2 geno counts

#########################################################
## Original results for this test case, on non-X chrom.
#########################################################
## X-chrom Burden p-val is .00515
## we also do the Davies p-value for the kernel test, so will be slightly different
#$stat.kernel pedl
#[1] 80.10206
#$pval.kernel
#[1] 0.4039026
#$stat.burden
#[1,] 4.895617
#$pval.burden
#[1,] 0.02692495

## base case, m-b weights
pg.out.m2 <- pedgene(ped=example.ped, geno=example.geno, map=example.map, male.dose=2,
                     weights.mb=TRUE,checkpeds=TRUE)

# summary/print and plot methods for this object
print.pedgene(pg.out.m2,digits=3)
## standard result to compare against, note pval.kernel.davies different
##    gene chrom stat.kernel pval.kernel stat.burden pval.burden
## AA   AA     1        80.1       0.404        4.90     0.02692
## AX   AX     X       198.2       0.186        7.82     0.00515

## base case, beta weights, no pedcheck
pg.beta.m2 <- pedgene(ped=example.ped, geno=example.geno, map=example.map, male.dose=2, verbose.return=TRUE)
names(pg.beta.m2)
lapply(pg.beta.m2$save, dim)

print(pg.beta.m2, digits=4)

## base case, mb weights, method=kounen, no pedcheck
pg.kounen.m2 <- pedgene(ped=example.ped, geno=example.geno, map=example.map, male.dose=2, weights.mb=TRUE,method="kounen")
print(pg.kounen.m2,digits=4)

## try making ped1 disconeeded by taking 2nd-generation parents away
## results will differ a little
pg.out.m2.rm34 <- pedgene(example.ped[-(3:4),], example.geno, example.map, male.dose=2, checkpeds=FALSE, weights.mb=TRUE)
pg.out.m2.rm34

## Test character ids, which is robust now because we're now making super-ids by
## pasting ped-person together within the function
options(stringsAsFactors=FALSE)
char.ped <- with(example.ped, data.frame(ped=as.character(ped), person=as.character(person), father=as.character(father), mother=as.character(mother), sex=sex, trait=trait))
options(stringsAsFactors=TRUE)

## as long as subject and ped ids are character, not factor, this will work
## pedgene makes sure to not treat character as factor 
pg.out.m2.char <- pedgene(char.ped, example.geno, example.map, male.dose=2, checkpeds=FALSE)
pg.out.m2.char

## show that it accepts 23 as X, but recodes 23 to X within the function
map23 <- example.map
map23$chrom[map23$chrom=="X"] <- 23
pg.X23.m2 <- pedgene(ped=example.ped, geno=example.geno, map=map23, male.dose=2,
                     weights=NULL, checkpeds=TRUE)

print(pg.X23.m2, digits=3)


## geno row with all NA
geno.narow <- example.geno
geno.narow[4,3:ncol(example.geno)] <- NA
# to check if male dose>1 for males on X chrom -- works 
#geno.narow[3,3:ncol(example.geno)] <- ifelse(geno.narow[3,2:ncol(example.geno)]==0,0,2)
pg.narow.m2 <- pedgene(ped=example.ped, geno=geno.narow, map=example.map, male.dose=2,
                     weights=NULL, weights.mb=TRUE,checkpeds=TRUE)
print(pg.narow.m2,digits=3)

## choose marker 4 as the 1-snp to represent the gene
## single-snp genes reduce kernel test to burden, check that p-vals agree, stat.kern=stat.burd^2
## This also caused indexing problems in the past.
pg.g1.m2 <- pedgene(ped=example.ped, geno=example.geno[,c(1:2,6,13:22)],
                    map=example.map[c(1,11:20),], male.dose=2,weights.mb=TRUE)
pg.g1.m2

# male dose=1
pg.out.m1 <- pedgene(example.ped, example.geno, example.map, male.dose=1, )

print(pg.out.m1, digits=3)
##    gene chrom stat.kernel pval.kernel stat.burden pval.burden
## AA   AA     1      80.102      0.4039      4.8956    0.026925
## AX   AX     X      49.140      0.3568      5.2908    0.021438


## test with no map arg (all variants in one gene columns 3:12)
pg.out.nomap <- pedgene(example.ped, example.geno[,1:12])
pg.out.nomap
#      gene   chrom stat.kernel pval.kernel.davies stat.burden pval.burden
# 1 unknown unknown      80.102             0.3926      4.8956    0.026925

## test with extra subject in geno, make sure it is removed
example2.geno <- rbind(example.geno[1,],example.geno)
pg.out <- pedgene(ped=example.ped, geno=example2.geno, map=example.map, male.dose=2,
                     weights.mb=TRUE,checkpeds=TRUE, method=NA)
warnings()
example2.geno[1,1:2] <- c(10,5)
pg.out <- pedgene(ped=example.ped, geno=example2.geno, map=example.map, male.dose=2,
                     weights.mb=TRUE,checkpeds=TRUE)
warnings()
pg.out

## Testing first gene with dose=2-dose
geno.recode <- cbind(example.geno[,1:2], 2-example.geno[,grep("AA", names(example.geno))])
pg.recode.mb <- pedgene(example.ped, geno.recode, male.dose=2, weights.mb=TRUE)
## note when map not given, assumes all 1 gene, and assigns "unknown" gene/chrom
pg.recode.mb

pg.recode.beta <- pedgene(example.ped, geno.recode, male.dose=2)
## note when map not given, assumes all 1 gene, and assigns "unknown" gene/chrom
pg.recode.beta


## weights, Madsen-Browning
maf <- colMeans(example.geno[,-(1:2)]/2)
## maf not correct for X matrix, b/c n-alleles for males is not 2
## so these results will be a little different for X-chrom

pg.out.wts.m2 <- pedgene(example.ped, example.geno, map=example.map,
         weights=1/sqrt((maf*(1-maf))))
# note stat, pval for AX gene don't match pg.out.m2
print(pg.out.wts.m2)


## one column genotype
pg.out.1snp <- pedgene(example.ped, example.geno[,c(1,2,4),drop=FALSE], map=example.map[2,,drop=FALSE])
pg.out.1snp

## plot, consider using the unrelated kernel-clustering plot method to show
##       regions of clustering more than expected,
##       plot gene regions separately


## Testing many genes at once:

genobig <- example.geno
mapbig <- example.map
for(k in 2:10) {
  genobig <- cbind(genobig, example.geno[,-(1:2)])
  mapbig <- rbind(mapbig, example.map)
  mapbig$gene[((k-1)*20+1):(20*k)] <- paste(example.map$gene[1:20],k,sep="")
}

## Add two genes: one with 1 variant. Another with no markers with variance
genobig <- cbind(genobig, example.geno[,6], rep(1, nrow(example.geno)), rep(2, nrow(example.geno)))
mapbig <- rbind(mapbig, c(10, "onevar"), c(11,"novar"), c(11, "novar"))
                

pgbig.m2 <- pedgene(example.ped, genobig, mapbig, male.dose=2, acc.davies=1e-4)
pgbig.m1 <- pedgene(example.ped, genobig, mapbig, male.dose=1, acc.davies=1e-4)

print(pgbig.m2, digits=3)
print(pgbig.m1, digits=3)

