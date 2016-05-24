### R code from vignette source 'vignette.Rnw'

###################################################
### code chunk number 1: vignette.Rnw:11-12
###################################################
library(ARTP)


###################################################
### code chunk number 2: vignette.Rnw:20-24
###################################################
pheno_file <- system.file("sampleData", "pheno_data.txt", package="ARTP")
geno_file  <- system.file("sampleData", "geno_data.txt", package="ARTP")
print(pheno_file)
print(geno_file)


###################################################
### code chunk number 3: vignette.Rnw:30-32
###################################################
pheno.list <- list(file=pheno_file, delimiter="\t", header=1, id.var="ID", 
                   response.var="Y", main.vars=c("X1", "X2"))


###################################################
### code chunk number 4: vignette.Rnw:38-39
###################################################
geno.list <- list(file=geno_file, delimiter="\t", file.type=2)


###################################################
### code chunk number 5: vignette.Rnw:44-46
###################################################
out.dir <- getwd()
print(out.dir)


###################################################
### code chunk number 6: vignette.Rnw:51-53
###################################################
gs_file  <- system.file("sampleData", "gene_SNP_data.txt", package="ARTP")
print(gs_file)


###################################################
### code chunk number 7: vignette.Rnw:57-58
###################################################
gs.list <- list(file=gs_file, snp.var="SNP", gene.var="Gene", delimiter="\t", header=1)


###################################################
### code chunk number 8: vignette.Rnw:65-67
###################################################
obs.outfile  <- paste(out.dir, "/", "obs.txt", sep="")
perm.outfile <- paste(out.dir, "/", "perm.txt", sep="") 


###################################################
### code chunk number 9: vignette.Rnw:72-74
###################################################
nperm   <- 50
op.list <- list(nperm=nperm, obs.outfile=obs.outfile, perm.outfile=perm.outfile, perm.method=2)


###################################################
### code chunk number 10: vignette.Rnw:78-79
###################################################
runPermutations(geno.list, pheno.list, 1, op=op.list)


###################################################
### code chunk number 11: vignette.Rnw:85-88
###################################################
set.seed(76523)
ret <- ARTP_pathway(obs.outfile, perm.outfile, nperm, out.dir, gene.list=gs.list)
print(ret)


###################################################
### code chunk number 12: vignette.Rnw:93-96
###################################################
set.seed(76523)
ret <- ARTP_pathway(obs.outfile, perm.outfile, nperm, out.dir)
print(ret)


