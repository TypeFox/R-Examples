## Raw data may be retrieved from https://tcga-data.nci.nih.gov/docs/publications/brca_2012/:
##
## http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.methylation.27k.450k.466.zip
## http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.GISTIC2.tar.gz
## http://tcga-data.nci.nih.gov/docs/publications/brca_2012/BRCA.exp.466.med.txt
##
## http://supportres.illumina.com/documents/myillumina/b78d361a-def5-4adb-ab38-e8990625f053/humanmethylation450_15017482_v1-2.csv

dataSet <- "tcga2012brca"

library("R.utils")
path <- "data"
path <- Arguments$getReadablePath(path)

## - - - - - - - - - - - - - - - - - - - - - - -
## 1. Methylation
## - - - - - - - - - - - - - - - - - - - - - - -

## 1a. data
filename <- "BRCA.methylation.27k.450k.466.txt"
pathname <- file.path(path, filename)

mdat <- read.table(pathname, nr=2, sep="\t", check.names=FALSE)
dim(mdat)
## [1] 21986   466
## one line per probe present in both 27k and 450k datasets (?)

mpnames <- names(mdat)
mpnames <- substr(mpnames, 1, 15)
str(mpnames)

rownames(mdat)

mdat <- read.table(pathname, sep="\t", check.names=FALSE)
names(mdat) <- mpnames
methMat <- as.matrix(mdat)

mpnames <- names(mdat)
mpnames <- substr(mpnames, 1, 15)
str(mpnames)

## 1b. annotation
afilename <- "HumanMethylation450_15017482_v.1.2.csv"
apathname <- file.path(path, afilename)

adat <- read.table(apathname, sep=",", skip=7, quote="\"", nr=485577, header=TRUE, as.is=TRUE)
names(adat)

## focus on probes that are in the microarray
cgNames <- rownames(mdat)  ## names of probes on the microarray
mm <- match(cgNames, adat[["Name"]])
(length(mm)==length(cgNames))  ## TRUE: all probes in 'mdat' are also in the annotation data file

adat2 <- adat[mm, ]
geneNames <- adat2[, "UCSC_RefGene_Name"]
str(geneNames)
## chr [1:21986] "ATP2A1;ATP2A1" "SLMAP" "MEOX2" "HOXD3" "PANX1" ...
mgnames <- unique(unlist(strsplit(geneNames, ";")))
str(mgnames)
## chr [1:13956] "ATP2A1" "SLMAP" "MEOX2" "HOXD3" "PANX1" "COX8C" ...

## - - - - - - - - - - - - - - - - - - - - - - -
## 2. Copy number
## - - - - - - - - - - - - - - - - - - - - - - -

## data and annotation
cpathname <- "data/GISTIC2/all_data_by_genes.txt"
cdat <- read.table(cpathname, sep="\t", as.is=TRUE, header=TRUE, check.names=FALSE)

cnMat <- as.matrix(cdat[, -c(1:3)])

cpnames <- colnames(cnMat)
## cpnames <- substr(cpnames, 1, 15)
str(cpnames)

cgnames <- cdat[["Gene Symbol"]]
str(cgnames)

## - - - - - - - - - - - - - - - - - - - - - - -
## 3. Expression
## - - - - - - - - - - - - - - - - - - - - - - -
## expression gene and patient names
epathname <- "data/BRCA.exp.466.med.txt"
edat <- read.table(epathname, sep="\t", as.is=TRUE, header=TRUE, quote="\"", check.names=FALSE)
egnames <- edat[["NAME"]]
str(egnames)

## data
exprMat <- as.matrix(edat[, -1])
colnames(exprMat) <- substr(colnames(exprMat), 1, 15)
epnames <- colnames(exprMat)
str(epnames)

## gene names
length(intersect(cgnames, mgnames))  ## 13261!
length(intersect(cgnames, egnames))  ## 15466!
length(intersect(intersect(cgnames, mgnames), egnames))  ## 11943

gids <- intersect(intersect(cgnames, mgnames), egnames)  ## 15433!!

eg <- match(gids, egnames)
cg <- match(gids, cgnames)

## patient names
length(intersect(cpnames, mpnames)) 
length(intersect(cpnames, epnames)) 
length(intersect(intersect(cpnames, mpnames), epnames))

pids <- intersect(intersect(cpnames, mpnames), epnames)

ep <- match(pids, epnames)
cp <- match(pids, cpnames)
mp <- match(pids, mpnames)

## methylation: mapping 'cg' ids to gene names
x <- strsplit(geneNames, ";")
x <- lapply(x, unique)
xx <- sapply(x, paste, collapse=";")
length(xx)

## a bit ad hoc
tf1 <- tempfile()
tf2 <- tempfile()
tf2 <- tempfile(tmpdir=".")
write.table(xx, file=tf1, row.names=FALSE, quote=FALSE, col.names=FALSE)
pl <- system.file("testScripts/tcga2012brca/findGenesIndices.pl", package="tmle.npvi")
system(paste(pl, tf1, ">", tf2))
y <- read.table(tf2, sep="\t", header=FALSE, as.is=TRUE)
names(y) <- c("name", "char")
head(y)
## /ad hoc

dim(y)
length(intersect(y$name, gids))==length(gids)

mg <- match(gids, y$name)

if (FALSE) { ## sanity check
  ids <- as.numeric(unlist(strsplit(y$char[mg[1]], " ")))
  unique(unlist(x[ids]))==gids[1]

  lens <- sapply(y$char, FUN=function(x) length(unlist(strsplit(x, " "))))
  table(lens)
}

## - - - - - - - - - - - - - - - - - - - - - - -
## 4. Data export: one gene
## - - - - - - - - - - - - - - - - - - - - - - -
gid <- gids[1]

## gene expression
idxE <- eg[match(gid, gids)]
geneExpr <- exprMat[idxE, ep]

## DNA copy number
idxC <- cg[match(gid, gids)]
copyNumber <- cnMat[idxC, cp]

## methylation
idxM <- mg[match(gid, gids)]
idxsM <- as.numeric(unlist(strsplit(y[idxM, 2], " ")))
methyl <- methMat[idxsM, mp, drop=FALSE]
dim(methyl)
stopifnot(identical(names(copyNumber), colnames(methyl)))

rownames(methyl) <- paste("W", 1:nrow(methyl), sep="")

obs <- cbind(Y=geneExpr, X=copyNumber, W=t(methyl))
str(obs)

## pairs(obs)

saveObject(obs, file="obs,ATAD3A,2.xdr")

## - - - - - - - - - - - - - - - - - - - - - - -
## 5. Data export: many genes
## - - - - - - - - - - - - - - - - - - - - - - -

chr <- 1:22

## Retrieve gene names and positions from biomaRt
if (!require("biomaRt")) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
}
library("biomaRt")
ensembl <- useMart("ensembl",dataset="hsapiens_gene_ensembl")
bdat <- getBM(attributes=c("hgnc_symbol", "chromosome_name", "start_position", "end_position"),
      filters=c("chromosome_name"),
      values=list(1:22), mart=ensembl)
names(bdat) <- c("name", "chr", "start", "end") 
head(bdat, 10)

geneNames <- bdat[["name"]]
inter <- intersect(geneNames, gids)

length(inter)
## [1] 11314
length(gids)
## [1] 11943

idxs <- match(inter, geneNames)
stopifnot(all(!is.na(idxs)))
bdat <- bdat[idxs, ]

outPath <- file.path("geneData", dataSet)
outPath <- Arguments$getWritablePath(outPath)
  
for (ch in chr) {
  print(ch)
  ## for one chromosome
  bdatCC <- subset(bdat, (chr==ch) & (name!=""))
  dim(bdatCC)

  for (ii in 1:nrow(bdatCC)) {
    gid <- bdatCC[ii, "name"]
    pos <- bdatCC[ii, "start"]
    me <- match(gid, gids)
    if (is.na(me)) {
      warning("Gene name not found: ", gid)
      next;
    }  

    ## gene expression
    idxE <- eg[me]
    geneExpr <- exprMat[idxE, ep]

    ## DNA copy number
    idxC <- cg[match(gid, gids)]
    copyNumber <- cnMat[idxC, cp]

    ## methylation
    idxM <- mg[match(gid, gids)]
    idxsM <- as.numeric(unlist(strsplit(y[idxM, 2], " ")))
    methyl <- methMat[idxsM, mp, drop=FALSE]
    dim(methyl)
    stopifnot(identical(names(copyNumber), colnames(methyl)))
    if (nrow(methyl)==1) {
      rownames(methyl) <- "W"
    } else {
      rownames(methyl) <- paste("W", 1:nrow(methyl), sep="")
    }
    obs <- cbind(Y=geneExpr, X=copyNumber, W=t(methyl))
    ## str(obs)
    ## pairs(obs)
    filename <- sprintf("chr%s,%06d,%s.xdr", ch, round(pos/1e3), gid)
    pathname <- file.path(outPath, filename)
    
    saveObject(obs, file=pathname)
  }
}

## discretization of copy numbers
thr <- 2e-2
nz <- apply(cnMat[, cp], 1, function(x) sum(abs(x)<=thr))
summary(nz)
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
##  49.0   111.0   180.0   162.2   200.0   257.0
