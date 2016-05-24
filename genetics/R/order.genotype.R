### order.genotype.R
###------------------------------------------------------------------------
### What: Get order of genotype object according to order of allele names or
###       genotype names
### Time-stamp: <2007-07-20 10:47:44 ggorjan>
###------------------------------------------------------------------------

## Define order as generic function
order <- function(..., na.last=TRUE, decreasing=FALSE)
  UseMethod("order")

## Get default method for order from base package
order.default <- function(..., na.last=TRUE, decreasing=FALSE)
  base::order(..., na.last=na.last, decreasing=decreasing)

## Genotype/Haplotype methods
order.genotype <- function(..., na.last=TRUE, decreasing=FALSE,
                           alleleOrder=allele.names(x), genotypeOrder=NULL)
{
  x <- list(...)[[1]]
  isHap <- is.haplotype(x)
  reorder <- ifelse(isHap, "no", "yes")
  if (is.null(genotypeOrder)) {

    ## --- Sort by alleleOrder ---

    alleleOrder <- unique(alleleOrder)

    ## Prepair working genotype data for given alleles only
    if(!identical(alleleOrder, allele.names(x)))
      x <- genotype(x, alleles=alleleOrder, reorder=reorder)

    ## Prepair working data for sort i.e. 01_allele1/15_allele15
    tmp <- as.character(x)
    id <- seq(along=alleleOrder)
    id <- formatC(id, width=ceiling(log10(max(id))), flag="0")
    for(i in seq(along=alleleOrder)) {
      tmp <- sub(pattern=paste(alleleOrder[i], "/", sep=""),
                 replacement=paste(i, "_", alleleOrder[i], "/", sep=""),
                 x=tmp)
      tmp <- sub(pattern=paste("/", alleleOrder[i], sep=""),
                 replacement=paste("/", i, "_", alleleOrder[i], sep=""),
                 x=tmp)
    }
  } else {

    ## --- Sort by genotypeOrder ---

    genotypeOrder <- unique(genotypeOrder)

    if(!isHap) {
      ## Match both A/B and B/A
      genotypeOrder <- .genotype2Haplotype(x=genotypeOrder)
    }
    tmp <- match(x, genotypeOrder)
  }
  ## print(tmp)
  return(order(tmp, na.last=TRUE, decreasing=FALSE))
}

sort.genotype <- function(x, decreasing=FALSE, na.last=NA, ...,
                          alleleOrder=allele.names(x), genotypeOrder=NULL)
{
  x[order(x, decreasing=decreasing, na.last=na.last,
          alleleOrder=alleleOrder, genotypeOrder=genotypeOrder)]
}

## No need for haplotype methods as they are exactly the same and haplotype
## is extended class of genotype sort.haplotype

###------------------------------------------------------------------------
### order.genotype.R ends here
