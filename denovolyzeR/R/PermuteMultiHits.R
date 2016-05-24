#' Permutes x variants across a genelist, and counts genes with multiple hits
#'
#' An internal function called by denovolyzeMultiHits
#'
#' @inheritParams denovolyze
#' @param x Total number of de novo variants observed in dataset
#' @param y Number of genes with >1 de novo variant (of class "class") in the population
#' @param nperms Number permutations
#' @param class In c("lof","mis","syn","prot")
#'
#' @seealso \code{\link{denovolyzeMultiHits}}
#'
#' @return Returns a named vector of 5 values
#'

# --------------------

PermuteMultiHits <- function(x,y,nperms=100,
                             class="lof",
                             geneId="geneName",
                             includeGenes="all",
                             probTable=pDNM) {

  #x = total number of DNM observed
  #y = no of genes with >1 DNM in class of interest
  #nperms = number permutations, defaults to 100
  #class = type of DNM assessed, defaults to "lof"

  # --------------------
  # Use specified gene ID
  names(probTable)[names(probTable)==geneId] <- "gene"
  probTable$gene <- toupper(as.character(probTable$gene))
  includeGenes <- toupper(as.character(includeGenes))

  # --------------------
  # If a list of genes for inclusion is specified, restrict analysis to these genes
  if(includeGenes[1]!="ALL"){
    probTable <- probTable[probTable$gene %in% includeGenes,]
  }
  probtable <- probTable[probTable$class==class,c("gene","value")]
  mycounts<- NA
  for (i in 1:nperms) {
    DNMsim <- sample(probtable$gene,x,replace=T,prob=probtable$value)
    mycounts[i] <- length(unique(DNMsim[duplicated(DNMsim)]))
  }
  empirical.p<-length(which(mycounts >= y))/nperms

  # --------------------
  output <- vector(length=5)
  output <-c(y,mean(mycounts),max(mycounts),empirical.p,x)
  names(output) <- c("obs","expMean","expMax","pValue","nVars")
  return(output)

}
