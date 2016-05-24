
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## Tests read.gtf features
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##
## hs.ensembl.76.small.gtf
## contains the first 100 lines of file
## ftp://ftp.ensembl.org/pub/release-76/gtf/
##                                  homo_sapiens/Homo_sapiens.GRCh38.76.gtf.gz
## + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + ##

ens76 <- system.file("extdata",
            "hs.ensembl.76.small.gtf", package = "refGenome", mustWork=TRUE)

en76 <- ensemblGenome()
basedir(en76) <- dirname(ens76)
read.gtf(en76, basename(ens76))

if(!exists("genes", where=en76@ev, inherits=FALSE))
    stop("[test_read_gtf] 'genes' table does not exist")

if(!all(getGeneTable(en76)$gene_name[1:2]==c("DDX11L1","WASH7P")))
    stop("[test_read_gtf] Wrong gene names in 'genes' table.")
