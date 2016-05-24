## ----htmlvignette, echo=FALSE--------------------------------------------
html_vignette <- function(fig_width = 6,
                          fig_height = 6,
                          dev = 'png',
                          css = NULL,
                          ...) {
  
  if (is.null(css)) {
    css <- system.file("rmarkdown", "templates", "html_vignette" ,"resources", 
      "vignette.css", package = "rmarkdown")
  }
  
  html_document(fig_width = fig_width, 
                fig_height = fig_height, 
                dev = dev,
                fig_retina = FALSE,
                css = css, 
                theme = NULL,
                highlight = "pygments",
                ...)
}


## ----loadlibs, echo=FALSE, error=FALSE, warning=FALSE, message=FALSE, include=FALSE----
#options(RCurlOptions=list(proxy="www-int2.inet.dkfz-heidelberg.de:80", http.version=1))
options(repos = c(CRAN="http://cran.r-project.org"))
if("caRpools" %in% rownames(installed.packages()) == FALSE) {install.packages("caRpools")}
  library(caRpools,warn.conflicts = FALSE, quietly = TRUE,verbose =FALSE)
load.packages()
data(caRpools)

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  ## for data extraction of a third replicate or more
#  
#  # First edit the MIACCS file loading
#  fileCONTROL3 = as.character(miaccs["carpools.untreated3",3])
#  d.CONTROL3 = as.character(miaccs["carpools.untreated3.desc",3])
#  
#  fileTREAT3 = as.character(miaccs["carpools.treated3",3])
#  d.TREAT3 = as.character(miaccs["carpools.treated3.desc",3])
#  
#  # Additional CONTROL replicate
#  fileCONTROL3 = data.extract(scriptpath, datapath, fastqfile=fileCONTROL3, extract, seq.pattern, maschine.pattern, createindex, bowtie2file, referencefile, mapping, reversecomplement, threads, bowtieparams, sensitivity, match)
#  # Additional Treatment replicate
#  fileTREAT3 = data.extract(scriptpath, datapath, fastqfile=fileTREAT3, extract, seq.pattern, maschine.pattern, createindex, bowtie2file, referencefile, mapping, reversecomplement, threads, bowtieparams, sensitivity, match)
#  
#  ## If you do not want to extract FASTQ files, but just load a read-count file this is the way to go
#  
#  CONTROL3 = load.file(paste(datapath, fileCONTROL3, sep="/"))
#  TREAT3 = load.file(paste(datapath, fileTREAT3, sep="/"))
#  

## ---- eval=FALSE, echo=TRUE----------------------------------------------
#  
#  ## WILCOX
#  # Provide more replicates as list
#  data.wilcox = stat.wilcox(untreated.list = list(CONTROL1, CONTROL2, CONTROL3, CONTROL4), treated.list = list(TREAT1,TREAT2, TREAT3, TREAT4), namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, normalize=normalize, norm.fun=norm.function, sorting=FALSE, controls=controls.nontarget, control.picks=control.picks)
#  
#  ## DESEQ2
#  data.deseq = stat.DESeq(untreated.list =  list(CONTROL1, CONTROL2, CONTROL3, CONTROL4), treated.list = list(TREAT1,TREAT2, TREAT3, TREAT4), namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, extractpattern=g.extractpattern, sorting=FALSE, filename.deseq = paste(analysis.name, "-ANALYSIS-DESeq2-sgRNA.tab", sep=""), sgRNA.pval = sig.pval.deseq, fitType="mean")
#  
#  ## MAGECK
#  data.mageck = stat.mageck(untreated.list = list(CONTROL1, CONTROL2, CONTROL3, CONTROL4), treated.list = list(TREAT1,TREAT2, TREAT3, TREAT4), namecolumn=namecolumn, fullmatchcolumn=fullmatchcolumn, norm.fun="median", extractpattern=g.extractpattern, mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr", filename = paste(analysis.name, "-ANALYSIS-MAGeCK-RAW", sep=""), fdr.pval=sig.pval.mageck)
#  

## ----load-file, echo=TRUE, eval=FALSE------------------------------------
#  CONTROL1 = load.file("./data/untreated1.txt", header= TRUE, sep="\t")

## ----load-files-read-count, echo=TRUE, eval=FALSE------------------------
#  CONTROL1 = load.file("./data/untreated1.txt", header= TRUE, sep="\t")
#  CONTROL = load.file("./data/untreated2.txt", header= TRUE, sep="\t")
#  TREAT1 = load.file("./data/treated1.txt", header= TRUE, sep="\t")
#  TREAT2 = load.file("./data/treated2.txt", header= TRUE, sep="\t")
#  
#  # Don't forget the library reference
#  libFILE = load.file(
#    paste(datapath, paste(referencefile,".fasta",sep=""), sep="/"),header = FALSE,
#    type="fastalib")

## ----load-files-fastq, echo=TRUE, eval=FALSE-----------------------------
#  fileCONTROL1 = data.extract(
#    scriptpath="path.to.scripts", datapath="path.to.FASTQ",fastqfile="filename1",
#    extract=TRUE, seq.pattern, maschine.pattern, createindex=TRUE,
#    bowtie2file=filename.lib.reference, referencefile="filename.lib.reference",
#    mapping=TRUE, reversecomplement=FALSE, threads, bowtieparams,
#    sensitivity="very-sensitive-local",match="perfect")
#  # Now we can load the generated Read-Count file directly!
#  CONTROL1 = load.file(paste(datapath, fileCONTROL1, sep="/")) # Untreated sample 1 loaded
#  
#  # Don't forget the library reference
#  libFILE = load.file(
#    paste(datapath, paste(referencefile,".fasta",sep=""), sep="/"),header = FALSE,
#    type="fastalib")

## ----aggregate1, echo=TRUE, eval=FALSE-----------------------------------
#  CONTROL1.g=aggregatetogenes(data.frame = CONTROL1, agg.function=sum,
#          extractpattern = expression("^(.+?)(_.+)"),type = "aggregate")

## ----aggregate-example1, echo=TRUE, eval=FALSE---------------------------
#  CONTROL1.g=aggregatetogenes(data.frame = CONTROL1, agg.function=sum,
#              extractpattern = expression("^(.+?)(_.+)"), type="aggregate")
#  knitr::kable(CONTROL1.g[1:10,])

## ----aggregate-example2, echo=TRUE, eval=TRUE----------------------------
CONTROL1.g.annotate=aggregatetogenes(data.frame = CONTROL1, agg.function=sum,
            extractpattern = expression("^(.+?)(_.+)"), type="annotate")
knitr::kable(CONTROL1.g.annotate[1:10,])

## ----get.gene.info, echo=TRUE,eval=FALSE---------------------------------
#  # Convert HGNC gene symbol to EnsemblID
#  CONTROL1 = get.gene.info(
#    CONTROL1, namecolumn=1, extractpattern=expression("^(.+?)(_.+)"),
#    host="www.ensembl.org",
#    database="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", filters="hgnc_symbol",
#    attributes = "ensembl_gene_id", return.val = "dataset")
#  # Also convert the library reference
#  libFILE = get.gene.info(
#    libFILE, namecolumn=1, extractpattern=expression("^(.+?)(_.+)"),
#    host="www.ensembl.org",
#    database="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", filters="hgnc_symbol",
#    attributes = "ensembl_gene_id", return.val = "dataset")

## ----geneinfo-replace, echo=TRUE, eval=TRUE------------------------------
CONTROL1.replaced = get.gene.info(CONTROL1, namecolumn=1, host="www.ensembl.org",
      extractpattern=expression("^(.+?)(_.+)"), database="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
      filters="hgnc_symbol", attributes =c("ensembl_gene_id"), return.val = "dataset")

knitr::kable(CONTROL1.replaced[1:10,])

## ----geneinfo-replace2, echo=TRUE, eval=TRUE-----------------------------
CONTROL1.replaced.info = get.gene.info(CONTROL1, namecolumn=1, host="www.ensembl.org",
      extractpattern=expression("^(.+?)(_.+)"), database="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl",
      filters="hgnc_symbol", attributes = c("ensembl_gene_id","description"), return.val = "info")

knitr::kable(CONTROL1.replaced.info[1:10,])

## ----stats-general-pre, echo=FALSE, eval=TRUE----------------------------
# General
U1.stats = stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
                      extractpattern=expression("^(.+?)_.+"), type="stats")
U2.stats = stats.data(dataset=CONTROL2, namecolumn = 1, fullmatchcolumn = 2,
                      extractpattern=expression("^(.+?)_.+"), type="stats")
T1.stats =stats.data(dataset=TREAT1, namecolumn = 1, fullmatchcolumn = 2,
                     extractpattern=expression("^(.+?)_.+"), type="stats")
T2.stats =stats.data(dataset=TREAT2, namecolumn = 1, fullmatchcolumn = 2,
                     extractpattern=expression("^(.+?)_.+"), type="stats")
# Combine Stats
combined.stats = cbind.data.frame(U1.stats[,1:2], U2.stats[,2], T1.stats[,2], T2.stats[,2])
colnames(combined.stats) = c("Readcount", d.CONTROL1, d.CONTROL2, d.TREAT1, d.TREAT2)

## ----stats-general, echo=TRUE, eval=FALSE--------------------------------
#  # General
#  U1.stats = stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
#                        extractpattern=expression("^(.+?)_.+"), type="stats")
#  U2.stats = stats.data(dataset=CONTROL2, namecolumn = 1, fullmatchcolumn = 2,
#                        extractpattern=expression("^(.+?)_.+"), type="stats")
#  T1.stats =stats.data(dataset=TREAT1, namecolumn = 1, fullmatchcolumn = 2,
#                       extractpattern=expression("^(.+?)_.+"), type="stats")
#  T2.stats =stats.data(dataset=TREAT2, namecolumn = 1, fullmatchcolumn = 2,
#                       extractpattern=expression("^(.+?)_.+"), type="stats")
#  # Combine Stats
#  combined.stats = cbind.data.frame(U1.stats[,1:2], U2.stats[,2], T1.stats[,2], T2.stats[,2])
#  colnames(combined.stats) = c("Readcount", d.CONTROL1, d.CONTROL2, d.TREAT1, d.TREAT2)
#  
#  # output all to a file
#  xlsx::write.xlsx(combined.stats, file="General-stats.xls", sheetName="Combined Stats", row.names=FALSE)
#  
#  # Create stats for whole data set
#  U1.allstats = stats.data(dataset=CONTROL1, namecolumn = namecolumn,
#      fullmatchcolumn = fullmatchcolumn, extractpattern=g.extractpattern, type="dataset")

## ----unmapped-readcount, echo=TRUE, eval=TRUE----------------------------
knitr::kable(stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
  extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = 1786217, type="mapping"))

## ----dataset-readcount, echo=TRUE, eval=TRUE-----------------------------
knitr::kable(stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
  extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = 1786217,
  type="stats"))

## ----list-readcount, echo=TRUE, eval=TRUE--------------------------------
knitr::kable(stats.data(dataset=CONTROL1, namecolumn = 1, fullmatchcolumn = 2,
  extractpattern=expression("^(.+?)_.+"), readcount.unmapped.total = 1786217,
  type="dataset")[1:10,1:5])

## ----stats-dropouts, echo=TRUE-------------------------------------------
# Calculate the number
U1.unmapped = unmapped.genes(data=CONTROL1, namecolumn=1, fullmatchcolumn=2, genes=NULL, extractpattern=expression("^(.+?)_.+"))

# output to a file
xlsx::write.xlsx(U1.unmapped, file="DROPOUT.xls", sheetName=d.CONTROL1, row.names=FALSE)

## ----stats-dropouts1, echo=TRUE------------------------------------------
# Calculate the number
U1.unmapped = unmapped.genes(data=CONTROL1, namecolumn=1, fullmatchcolumn=2, genes="random", extractpattern=expression("^(.+?)_.+"))

## ----carpools.read.distribution, echo=TRUE, eval=FALSE, sanitize=TRUE----
#  carpools.read.distribution(CONTROL1, fullmatchcolumn=2, breaks=200,
#    title=d.CONTROL1, xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE)
#  carpools.read.distribution(CONTROL1, fullmatchcolumn=2, breaks=200,
#    title=d.CONTROL1, xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE, type="whisker")

## ----read-distribution, echo=TRUE,fig.height=4, fig.width=8, eval=TRUE----
# Histogram
carpools.read.distribution(CONTROL1, fullmatchcolumn=2,breaks=200,
  title=d.CONTROL1, xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE) 

## ----read-distribution-whisker, echo=TRUE,fig.height=4, fig.width=8, eval=TRUE----
# Whisker
carpools.read.distribution(CONTROL1, fullmatchcolumn=2,breaks=200,
  title=d.CONTROL1, xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE,
  type="whisker") 

## ----read-distribution-onlygene, echo=TRUE, fig.height=4, fig.width=8, eval=TRUE----
# Histogram
carpools.read.distribution(CONTROL1, fullmatchcolumn=2,breaks=200,
  title=d.CONTROL1, xlab="log2 Readcount", ylab="# sgRNAs",statistics=TRUE, plotgene="CASP8") 

## ----QC-readdepth, echo=TRUE, eval=FALSE, sanitize=TRUE------------------
#  # Plot for a singe sample
#  carpools.read.depth(datasets = list(CONTROL1), namecolumn=1 ,fullmatchcolumn=2,
#    dataset.names=list(d.CONTROL1), extractpattern=expression("^(.+?)_.+"),
#    xlab="Genes", ylab="Read Count per sgRNA",statistics=TRUE, labelgenes = NULL,
#    controls.target = NULL, controls.nontarget="random")
#  
#  # Plot for 4 samples at once
#  carpools.read.depth(datasets = list(CONTROL1,CONTROL2,TREAT1,TREAT2), namecolumn=1,
#    fullmatchcolumn=2, dataset.names=list(d.CONTROL1,d.CONTROL2,d.TREAT1,d.TREAT2),
#    extractpattern=expression("^(.+?)_.+"), xlab="Genes", ylab="Read Count per sgRNA",
#    statistics=TRUE, labelgenes = NULL, controls.target = NULL, controls.nontarget="random")

## ----read-depth-example1, echo=TRUE, eval=TRUE,fig.height=4, fig.width=8----
# Plot for a singe sample
carpools.read.depth(datasets = list(CONTROL1), namecolumn=1 ,fullmatchcolumn=2,
  dataset.names=list(d.CONTROL1), extractpattern=expression("^(.+?)_.+"),
  xlab="Genes", ylab="Read Count per sgRNA",statistics=TRUE, labelgenes = NULL,
  controls.target = "CASP8", controls.nontarget="random", waterfall=FALSE)

## ----read-depth-example2, echo=TRUE, eval=TRUE,fig.height=4, fig.width=8----
# Plot for a singe sample sorted to read depth
carpools.read.depth(datasets = list(CONTROL1), namecolumn=1 ,fullmatchcolumn=2,
  dataset.names=list(d.CONTROL1), extractpattern=expression("^(.+?)_.+"),
  xlab="Genes", ylab="Read Count per sgRNA",statistics=TRUE, labelgenes = NULL,
  controls.target = "CASP8", controls.nontarget="random", waterfall=TRUE)

## ----QC-despergene, echo=TRUE, eval=FALSE, sanitize=TRUE-----------------
#  # Plot sgRNA presence per gene for control 1
#  control1.readspergene = carpools.reads.genedesigns(CONTROL1, namecolumn=1, fullmatchcolumn=2, title=paste("% sgRNAs:", d.CONTROL1, sep=" "), xlab="% of sgRNAs present", ylab="# of Genes")

## ----reads-genedesigns1, echo=TRUE, eval=TRUE----------------------------
control1.readspergene = carpools.reads.genedesigns(CONTROL1, namecolumn=1, fullmatchcolumn=2, title=paste("% sgRNAs:", d.CONTROL1, sep=" "), xlab="% of sgRNAs present", ylab="# of Genes")

## ----QC-targeting-scatter, echo=TRUE, eval=FALSE, sanitize=TRUE, fig.width=8, fig.height=8----
#  # Highlight the non-targeting control in tweo dataset
#  carpools.read.count.vs(dataset=list(TREAT1,CONTROL1), dataset.names = c(d.TREAT1, d.CONTROL1),
#    pairs=FALSE, namecolumn=1, fullmatchcolumn=2, title="", pch=16,
#    normalize=TRUE, norm.function="median", labelgenes="random", labelcolor="blue",
#    center=FALSE, aggregated=FALSE)
#  # Highlight non-targeting control in all datasets
#  carpools.read.count.vs(dataset=list(TREAT1, TREAT2, CONTROL1, CONTROL2),
#    dataset.names = c(d.TREAT1, d.TREAT2, d.CONTROL1, d.CONTROL2), pairs=TRUE, namecolumn=1,
#    fullmatchcolumn=2, title="", pch=16, normalize=TRUE, norm.function="median",
#    labelgenes="random", labelcolor="blue", center=FALSE, aggregated=FALSE)

## ----readcount-vs1, echo=TRUE, eval=TRUE, fig.width=7, fig.height=7------
# Highlight the non-targeting control in two datasets
carpools.read.count.vs(dataset=list(TREAT1,CONTROL1), dataset.names = c(d.TREAT1, d.CONTROL1),
  pairs=FALSE, namecolumn=1, fullmatchcolumn=2, title="", pch=16,
  normalize=TRUE, norm.function=median, labelgenes="random", labelcolor="blue",
  center=FALSE, aggregated=FALSE)

## ----readcount-vs12, echo=TRUE, eval=TRUE, fig.width=7, fig.height=7-----
# Highlight the non-targeting control in two datasets with CENTERING
carpools.read.count.vs(dataset=list(TREAT1,CONTROL1), dataset.names = c(d.TREAT1, d.CONTROL1),
  pairs=FALSE, namecolumn=1, fullmatchcolumn=2, title="", pch=16,
  normalize=TRUE, norm.function=median, labelgenes="random", labelcolor="blue",
  center=TRUE, aggregated=FALSE)

## ----readcount-vs13, echo=TRUE, eval=TRUE, fig.width=7, fig.height=7-----
# Highlight the targeting control in two datasets on GENE level
carpools.read.count.vs(dataset=list(TREAT1.g,CONTROL1.g), dataset.names = c(d.TREAT1, d.CONTROL1),
  pairs=FALSE, namecolumn=1, fullmatchcolumn=2, title="", pch=16,
  normalize=TRUE, norm.function=median, labelgenes="CASP8", labelcolor="red",
  center=FALSE, aggregated=TRUE)

## ----readcount-vs2, echo=TRUE, eval=TRUE, fig.width=8, fig.height=8------
# Highlight non-targeting control in all datasets
carpools.read.count.vs(dataset=list(TREAT1, TREAT2, CONTROL1, CONTROL2),
  dataset.names = c(d.TREAT1, d.TREAT2, d.CONTROL1, d.CONTROL2), pairs=TRUE, namecolumn=1,
  fullmatchcolumn=2, title="", pch=16, normalize=TRUE, norm.function=median,
  labelgenes="random", labelcolor="blue", center=FALSE, aggregated=FALSE)

## ----readcount-vs22, echo=TRUE, eval=TRUE, fig.width=8, fig.height=8-----
# Highlight targeting control in all datasets on Gene Level
carpools.read.count.vs(dataset=list(TREAT1.g, TREAT2.g, CONTROL1.g, CONTROL2.g),
  dataset.names = c(d.TREAT1, d.TREAT2, d.CONTROL1, d.CONTROL2), pairs=TRUE, namecolumn=1,
  fullmatchcolumn=2, title="", pch=16, normalize=TRUE, norm.function=median,
  labelgenes="CASP8", labelcolor="blue", center=FALSE, aggregated=TRUE)

## ----plot.raw.genes, echo=TRUE, eval=FALSE-------------------------------
#  carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
#    treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
#    fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"),
#    do.plot=TRUE, log=FALSE, put.names=TRUE, type="foldchange" )

## ----plot.raw.genes-example1, echo=TRUE, eval=TRUE , fig.width=8, fig.height=4----
# Foldchange
p1 = carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
  fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), 
  do.plot=TRUE, log=FALSE, put.names=TRUE, type="foldchange" )

# Z-Ratio
p2 = carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
  fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), 
  do.plot=TRUE, log=FALSE, put.names=TRUE, type="z-ratio" )

# Read Count
p3 = carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
  fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), 
  do.plot=TRUE, log=FALSE, put.names=TRUE, type="readcount" )

# Violine plot
p4 = carpools.raw.genes(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1, TREAT2), genes="CASP8", namecolumn=1,
  fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)_.+"), 
  do.plot=TRUE, log=FALSE, put.names=TRUE, type="vioplot" )

## ----HT-method-wilcox, echo=TRUE, eval=TRUE------------------------------
data.wilcox = stat.wilcox(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1,TREAT2), namecolumn=1, fullmatchcolumn=2,
  normalize=TRUE, norm.fun=median, sorting=FALSE, controls="random",
  control.picks=NULL)

## ----HT-method-DESeq2, echo=TRUE, eval=TRUE, message=FALSE, warning=FALSE----
# DESeq2
data.deseq = stat.DESeq(untreated.list = list(CONTROL1, CONTROL2),
  treated.list = list(TREAT1,TREAT2), namecolumn=1,
  fullmatchcolumn=2, extractpattern=expression("^(.+?)(_.+)"),
  sorting=FALSE, filename.deseq = "ANALYSIS-DESeq2-sgRNA.tab",
  fitType="parametric")

## ----HT-method-MAGeCK, echo=TRUE-----------------------------------------
# MAGeCK
data.mageck = stat.mageck(untreated.list = list(CONTROL1, CONTROL2), treated.list = list(TREAT1,TREAT2), namecolumn=1, fullmatchcolumn=2, norm.fun="median", extractpattern=expression("^(.+?)(_.+)"), mageckfolder=NULL, sort.criteria="neg", adjust.method="fdr", filename = "TEST" , fdr.pval = 0.05)


## ----p-val-visualization, echo=TRUE, eval=TRUE, fig.width=7, fig.height=7----
carpools.waterfall.pval(type="mageck", dataset=data.mageck, pval=0.05, log=TRUE)

## ----HT-visualization-mageck, echo=TRUE, eval=FALSE, sanitize=TRUE-------
#  mageck.result = carpools.hitident(data.mageck, type="mageck",
#        title="MAgECK P-values", xlab="Genes", ylab="log10(p-value)",
#        print.names=TRUE, offsetplot=1.2, plot.p=0.05, separate=TRUE)

## ----HT-visualization-mageck-example1, echo=TRUE, eval=TRUE, sanitize=TRUE, fig.width=10, fig.height=6----
mageck.result = carpools.hitident(data.mageck, type="mageck", title="MAGeCK", inches=0.1, print.names=TRUE, plot.p=0.05, offsetplot=1.2, sgRNA.top=1)

## ----HT-visualization-wilcox-example1, echo=TRUE, eval=TRUE, sanitize=TRUE, fig.width=10, fig.height=6----
wilcox.result = carpools.hitident(data.wilcox, type="wilcox", title="Wilcox", inches=0.1, print.names=TRUE, plot.p=0.05, offsetplot=1.2, sgRNA.top=1)

## ----HT-overview-hits, echo=TRUE, eval=TRUE, sanitize=TRUE, fig.width=7, fig.height=7----
carpools.hit.overview(wilcox=data.wilcox, deseq=data.deseq, mageck=data.mageck,
    cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05,
    cutoff.override=FALSE, cutoff.hits=NULL, plot.genes="overlapping", type="enriched")

## ----HT-compare-enriched-overlap, echo=TRUE, sanitize=TRUE, fig.width=7, fig.height=7----
venn.enriched = compare.analysis(wilcox=data.wilcox, deseq=data.deseq, mageck=data.mageck, type="enriched", cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL, output="venn")
require(VennDiagram)
plot.new()
grid::grid.draw(VennDiagram::venn.diagram(venn.enriched, file=NULL, fill=c("lightgreen","lightblue2","lightgray"), na="remove", cex=2,lty=2, cat.cex=2))

## ----HT-comparetable-enriched, echo=TRUE, sanitize=TRUE------------------
# Perform the comparison
data.analysis.enriched = compare.analysis(wilcox=data.wilcox,
    deseq=data.deseq, mageck=data.mageck, type="enriched",
    cutoff.override = FALSE, cutoff.hits=NULL, output="list",
    sort.by=c("mageck","fdr","rank"))
## Write to a file
xlsx::write.xlsx(data.analysis.enriched,
    file="COMPARE-HITS.xls",
    sheetName="Enriched")
# Print to console
knitr::kable(data.analysis.enriched[1:10,c(2:7)])

## ----HT-comparetable-enriched-ranked, echo=TRUE, sanitize=TRUE, fig.width=7, fig.height=7----
# Perform the comparison
data.analysis.enriched = compare.analysis(wilcox=data.wilcox,
    deseq=data.deseq, mageck=data.mageck, type="enriched",
    cutoff.override = FALSE, cutoff.hits=NULL, output="rank",
    sort.by=c("mageck","fdr","rank"))
## Write to a file
xlsx::write.xlsx(data.analysis.enriched,
    file="COMPARE-RANKs.xls",
    sheetName="Enriched")
# Print to console
knitr::kable(data.analysis.enriched[1:10,])

## ----HT-comparetable-enriched-3d, echo=TRUE, sanitize=TRUE, fig.width=7, fig.height=7----
# Perform the comparison and create 3D plot only for enriched
compared = compare.analysis( wilcox=data.wilcox, deseq=data.deseq, mageck=data.mageck, type="enriched", cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL, output="3dplot", sort.by=c("mageck","rank","rank"), plot.method=c("wilcox","deseq","mageck"), plot.feature=c("pval","pval","fdr"))

## ----final-table, echo=TRUE, sanitize=TRUE-------------------------------
# create final table with all genes + hit analysis results
final.tab = final.table(wilcox=data.wilcox, deseq=data.deseq, mageck=data.mageck, dataset=CONTROL1.g, namecolumn=1, type="genes")

# Write this to an EXCEL sheet
xlsx::write.xlsx(final.tab, "FINAL.xls", sheetName="Genelist", 
  col.names=TRUE, row.names=FALSE, append=FALSE, showNA=TRUE)

## ---- calculate-overlap, echo=TRUE, eval=TRUE----------------------------
overlap.enriched = generate.hits(wilcox=data.wilcox, deseq=data.deseq, mageck=data.mageck, type="enriched", cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL , plot.genes="overlapping")

## ----hit-scatter, echo=TRUE, eval=TRUE, fig.height=8, fig.width=8--------
plothitsscatter.enriched = carpools.hit.scatter(wilcox=data.wilcox, deseq=data.deseq, mageck=data.mageck, dataset=list(TREAT1, TREAT2, CONTROL1, CONTROL2), dataset.names = c(d.TREAT1, d.TREAT2, d.CONTROL1, d.CONTROL2), namecolumn=1, fullmatchcolumn=2, title="Title", labelgenes="CASP8", labelcolor="orange", extractpattern=expression("^(.+?)(_.+)"), normalize=TRUE, norm.function=median, offsetplot=1.2, center=FALSE, aggregated=FALSE, type="enriched", cutoff.deseq = 0.001, cutoff.wilcox = 0.05, cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL,  pch=16)

## ----carpools.hit.sgrna, echo=TRUE, eval=TRUE, fig.width=8, fig.height=6----
sgrnas.en = carpools.hit.sgrna(wilcox=data.wilcox, deseq=data.deseq,
    mageck=data.mageck, dataset=list(CONTROL1, CONTROL2, TREAT1, TREAT2),
    dataset.names = c(d.CONTROL1, d.CONTROL2, d.TREAT1, d.TREAT2), namecolumn=1,
    fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)(_.+)"),
    put.names=TRUE, type="enriched", labelgenes="CASP8", plot.type=NULL, 
    cutoff.deseq = 0.001, cutoff.wilcox=0.05, cutoff.mageck = 0.05,
    cutoff.override=FALSE, cutoff.hits=NULL, controls.target="CASP8", controls.nontarget="random")

## ----sgrna.table, echo=TRUE, eval=TRUE-----------------------------------
sgrnas.en.table = carpools.sgrna.table(wilcox=data.wilcox, deseq=data.deseq,
    mageck=data.mageck, dataset=list(CONTROL1, CONTROL2, TREAT1, TREAT2),
    dataset.names = c(d.CONTROL1, d.CONTROL2, d.TREAT1, d.TREAT2), namecolumn=1,
    fullmatchcolumn=2, norm.function=median, extractpattern=expression("^(.+?)(_.+)"),
    type="enriched", labelgenes="CASP8", cutoff.deseq = 0.001, cutoff.wilcox=0.05,
    cutoff.mageck = 0.05, cutoff.override=FALSE, cutoff.hits=NULL, sgrna.file = libFILE,
    write=FALSE)

