## package: megaptera
## called by: markersGenbank
## author: Christoph Heibl (at gmx.net)
## last update: 2014-07-15

EFetchLocus <- function (gi){
  
  ## retrieve sequence information
  ## -----------------------------
  x <- paste("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/", 
             "efetch.fcgi?db=nucleotide&id=", gi, 
             "&rettype=gb&retmode=xml", sep = "")
#   xml <- scan(x, what = "c", quiet = TRUE, sep = "\n")
#   write(xml, "test2.xml"); system("open -t test2.xml")
#       x  <- "test.xml"
  x <- xmlTreeParse(x, getDTD = FALSE, 
                    useInternalNodes = TRUE)
  
  ## get locus
  ## ---------
  gene <- xpathSApply(x, "//GBQualifier[GBQualifier_name='gene']/GBQualifier_value", xmlValue)
  if ( is.null(gene) ) 
    gene <- xpathSApply(x, "//GBQualifier[GBQualifier_name='product']/GBQualifier_value", xmlValue)
  if ( is.null(gene) ) 
    gene <- xpathSApply(x, "//GBQualifier[GBQualifier_name='note']/GBQualifier_value", xmlValue)
  unique(gene)
}