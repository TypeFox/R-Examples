#' An adjacency matrix of a random graph
#' 
#' An adjacency matrix of a random graph with some random Refseq
#' Protein IDs for use in example files and the vignette
#'
#' @name adjacency.matrix
#' @docType data
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @keywords data
NULL

#' Example gene expression data
#' 
#' A data matrix \code{x} containing gene expression data of 10 patients
#'
#' @name x
#' @docType data
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @keywords data
NULL

#' Example class labels for the gene expression data
#' 
#' Class labels for the 10 patients contained the data matrix \code{x}
#'
#' @name y
#' @docType data
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @keywords data
NULL

#' A mapping of Refseq Protein IDs to probe set IDs for the gene expression data
#' 
#' A mapping of the hgu95av2 probe set IDs in x to the Refseq protein IDs contained
#' in adjacency.matrix
#'
#' @name mapping
#' @docType data
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @keywords data
NULL


#' Summarize probe sets
#'
#' Summarize multiple probe sets targeting one gene into one value for that gene.
#' On most microarays there will be more than one probe set for a gene. However, in
#' the underlying network the gene will only be present one time. Therefore, in
#' order to calculate a Gene(Page)Rank weight for this gene, all expression measurements
#' have to be summarized.
#'
#' summarizes all probes of a gene to one
#' value for that gene
#' if the summarization method is 'none' then the only
#' thing which is done is that all probesets for which no
#' pathway is available are discarded.

#' @param exprs \eqn{n \times p}{n x p} matrix with \eqn{n} probe sets and \eqn{p} samples.
#' @param mapping a matrix or data.frame with 2 columns. The colnames of mapping have to contain at least 'graphID' and 'probesetID'.
#' These two columns define the mapping between the probe sets on the microarray and the nodes of the graph.
#' @param method defines how several probe sets should be combined. One of \code{median}, \code{mean}, \code{foldChange} or \code{none}.
#' @param groups defines the grouping of samples. Only needed if \code{method} is \code{foldChange}.
#' @param adjacency a matrix that represents the graph of the underlying biological network.
#' 
#' @return matrix with 1st comlumn probeIDs 2nd column gene IDs
summarizeProbes <- function(exprs, mapping, method="median", groups=NULL, adjacency=NULL){

  if(!all(c('graphID','probesetID') %in% colnames(mapping))) stop('\'mapping\' does not contain the needed columns.')
  if(is.null(adjacency)) stop("You have to provide the adjacency of the underlying network!\n")
  if(method == 'foldChange' & is.null(groups)) stop('Parameter \'groups\' is needed if \'method\' is \'foldChange\'.')

  genes <- graphID2probesetID <- gene.mean.exprs <- NULL
  
  probesets.with.pathway.knowledge <- unique(mapping[mapping[,"graphID"] %in% rownames(adjacency),"probesetID"])
  exprs <- exprs[rownames(exprs) %in% probesets.with.pathway.knowledge, ]
  graphID2probesetID <- sapply(rownames(adjacency), function(gene) mapping[mapping[,"graphID"] %in% gene,"probesetID"])

       if(method == "median")     gene.mean.exprs <- sapply(graphID2probesetID, function(probenames, exprs) median(rowMedians(exprs[as.character(probenames),,drop=F])),exprs=exprs)
  else if(method == "mean")       gene.mean.exprs <- sapply(graphID2probesetID, function(probenames, exprs) mean(rowMeans(exprs[as.character(probenames),,drop=F])),exprs=exprs)
  else if(method == "foldChange") gene.mean.exprs <- sapply(graphID2probesetID, function(probenames, exprs, groups) mean(rowMeans(exprs[as.character(probenames),groups == levels(groups)[1],drop=F]) - rowMeans(exprs[as.character(probenames),groups == levels(groups)[2],drop=F])),exprs=exprs, groups=groups)
  else if(method == "none")       gene.mean.exprs <- exprs
  else stop(paste("Summarization method ", method, " not implemented yet\n",sep=""))
  
  gene.mean.exprs
}

#' Desummarize GeneRanks back to the corresponding probesets
#'
#' Desummarize the GeneRanks which were previously calculated for each node in the
#' underlying biological network back to the corresponding probesets for the Reweighted
#' Recursive Feature Elimination (RRFE).
#'
#' @param ranks the previously calculated GeneRanks or PageRanks.
#' @param mapping a matrix or data.frame with 2 columns. The colnames of mapping have to contain at least 'graphID' and 'probesetID'.
#' 
#' @return matrix with 1st comlumn probeIDs 2nd column gene IDs
desummarize.ranks <- function(ranks, mapping){

  if(!all(c('graphID','probesetID') %in% colnames(mapping))) stop('\'mapping\' does not contain the needed columns.')
  
  ind            <- mapping[,"graphID"] %in% names(ranks)
  mapping        <- mapping[ind,]
  mapping        <- cbind(mapping, ranks=ranks[mapping[,"graphID"]])
  ranks.complete <- tapply(as.numeric(mapping[,"ranks"]), mapping[,"probesetID"], mean)
  ranks.complete
}

#' Parse the HPRD flat file
#'
#' This function parses the tab delimited flat file of protein-protein interactions
#' coming from the HPRD (\url{http://www.hprd.org/download}).
#'
#' @param fname path to the HPRD flat file.
#' @param chipProteins limit the resulting adjacency matrix to certain proteins.
#' 
#' @return An adjacency matrix
#' @export
#' @callGraphPrimitives
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' hprd <- read.hprd('BINARY_PROTEIN_PROTEIN_INTERACTIONS.txt')
#' }
read.hprd <- function(fname, chipProteins = NULL){
  if(!file.exists(fname))
    stop(paste("File: ",fname," does not exist.\n",sep=""))
  
  hprd=read.delim(fname, header=FALSE, as.is=TRUE)

  ## use RefSeq IDs
  hprd=as.matrix(hprd[,c(3,6)])

  ## convert: NP_000680.2 => NP_000680
  hprd=sub("\\..*$","",hprd)
  hprd=unique(hprd)

  proteins = NULL
  if(is.null(chipProteins))
    proteins=unique(as.character(hprd))
  else
    proteins = unique(chipProteins)
  
  hprd=matrix(match(hprd, proteins), ncol=2)
  hprd=hprd[!is.na(rowSums(hprd)),]


  hprdM=matrix(0,length(proteins),length(proteins))
  rownames(hprdM) = proteins
  colnames(hprdM) = proteins
  
  hprdM[hprd]=1
  hprdM[hprd[,c(2,1)]]=1
  diag(hprdM)=0

  hprdM
}


#' Matches the expression data to the adjacency matrix using the provided mapping.
#'
#' Usually the dimension of the graph and the expression data do not fit to each other.
#' Additionally most often the graph comprises another type of knowledge, i.e. the
#' expression matrix measures 10.000 genes represented as 15.000 probe sets and the graph
#' provides information on 7.000 proteins. Thus, a node (protein) of the graph might match to
#' two probe sets in the expression matrix (since both target the gene encoding the protein).
#' Therefore, this method uses the relationship between probe sets and i.e. proteins which is
#' encoded in the \code{mapping} to create a graph of probe sets rather then a graph of proteins.
#'
#' @param x the p x n expression matrix with p patients and n genes.
#' @param mapping a mapping which encodes the relationship between the colnames of x and the
#'        row/colnames of the adjacency matrix.
#' @param adjacency the adjacencymatrix of the underlying graph structure.
#' @return the matched input \item{x}{the expression matrix containing only the features which are also present in the adjacency matrix}
#'                           \item{mapping}{the mapping containing only necessary information}
#'                           \item{adjacency}{the adjacency matrix with the same number of nodes as features in x}
#' @export
#' @callGraphPrimitives
#' @author Marc Johannes \email{JohannesMarc@@gmail.com}
#' @examples
#' \dontrun{
#' library(Biobase)
#' data(sample.ExpressionSet)
#' x <- t(exprs(sample.ExpressionSet))
#' y <- factor(pData(sample.ExpressionSet)$sex)
#' # create the mapping
#' library('hgu95av2.db')
#' mapped.probes <- mappedkeys(hgu95av2REFSEQ)
#' refseq <- as.list(hgu95av2REFSEQ[mapped.probes])
#' times <- sapply(refseq, length)
#' mapping <- data.frame(probesetID=rep(names(refseq), times=times), graphID=unlist(refseq), row.names=NULL, stringsAsFactors=FALSE)
#' mapping <- unique(mapping)
#' library(pathClass)
#' data(adjacency.matrix)
#' matched <- matchMatrices(x=x, adjacency=adjacency.matrix, mapping=mapping)
#' }
matchMatrices <- function(x, mapping, adjacency){

  mapping <- unique(mapping)
  
  ## 1. fit mapping to chip. So dass das mapping nur noch probesets, gene oder was auch immer
  ##    auf dem chip ist enthaelt. Sollte ja eigentlich schon der Fall sein, nur ein sicherheits Schritt.
  int <- intersect(colnames(x), mapping[,'probesetID'])
  x <- x[, int, drop=FALSE]
  mapping <- mapping[mapping[,'probesetID'] %in% int,]

  ## 2. nur noch den intersect aus mapping (dass ja jetzt zu 100% den chip widerspiegelt) und der adjazenz matrix behalten
  int <- intersect(rownames(adjacency), mapping[,"graphID"])
  adjacency <- adjacency[int, int]
  mapping <- mapping[mapping[,'graphID'] %in% int,]
  
  ## 3. den chip nun dem mapping anpassen
  x <- x[,mapping[,'probesetID']]

  ## 4. a) dim von x > dim von adjacency -> dimensionen des chips zusammenfassen
  if(ncol(x) > ncol(adjacency)){
    graphID2probesetID <- sapply(rownames(adjacency), function(gene) mapping[mapping[,"graphID"] %in% gene,"probesetID"])
    gene.mean.exprs <- sapply(graphID2probesetID, function(probenames, exprs) rowMeans(exprs[,as.character(probenames),drop=F]),exprs=x)
  }
  ## 4. b) dim von adjacency > dim von x -> dimensionen des chips vergroessern
  else if(ncol(adjacency) > ncol(x)){
    probesetID2graphID <- sapply(colnames(x), function(gene) mapping[mapping[,"probesetID"] %in% gene,"graphID"])
    stop('dim von adjacency > dim von x -> dimensionen des chips vergroessern')
  }
  list(x=gene.mean.exprs, mapping=mapping, adjacency=adjacency)
}

## 

## ## get a mapping from probeset to protein to entrezID
## ## chip = Type of chip, e.g. HGU 133a
## get.mapping.for.chip <- function(chip = "hgu133a"){
## 
##   if(chip == "hgu133a"){
##     ## Get the probe identifiers that are mapped to an ENTREZ Gene ID
##     mapped.probes <- mappedkeys(hgu133a.db::hgu133aREFSEQ)
##     refseq <- as.list(hgu133a.db::hgu133aREFSEQ[mapped.probes])
## 
##     mapped.probes <- mappedkeys(hgu133a.db::hgu133aENTREZID)
##     entrez        <- as.list(hgu133a.db::hgu133aENTREZID[mapped.probes])
## 
##     map <- NULL
##     
##     for(i in 1:length(refseq)){
## 
##       probesetID      <- names(refseq[i])
##       refseqProteinID <- grep(".*?P.*",refseq[[i]],value=T)
## 
##       if(length(refseqProteinID) > 0){
##         entrezID <- entrez[probesetID]
##         tmp      <- cbind(probesetID, refseqProteinID, entrezID)
##         map      <- rbind(map, tmp)
##       }
##     }
##   }
##   else
##     stop(paste("Chip type: ", chip, " not implemented yet.\n", sep=""))
## 
##   map = apply(map, 2, as.character)
##   map
## }
