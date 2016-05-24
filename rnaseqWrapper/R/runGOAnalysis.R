runGOAnalysis <- function(
  sigGenes, # character string of sig genes
  expGenes, # character string of all genes
  goAnno, # GO annotation, as read by readMappings
  pValThresh = 1,
  plotGO = FALSE, ## Should a plot be generated
  ontology = "BP", ## Can be BP, MF, CC, or all or vector with multiple of BP/MF/CC
  algorithm = "weight", ## which algorithm should be used
  statistic = "fisher", ## Stat to use for GO test
  description = NULL # Description for the GO data object
  ){
  
  if(!require(topGO)){
    stop("The Bioconductor package 'topGO' is required to run this function. ",
          "Please install it now using:\n\n",
          'source("http://bioconductor.org/biocLite.R")\n',
          'biocLite("topGO")\n')
  }
  
  ## Ask whether each expGenes is in the sigGenes set
  ## This returns a TRUE/FALSE, that we modify into 1's and 0's
  genesComp <- factor(as.integer(expGenes %in% sigGenes))
  names(genesComp) <- expGenes ## Assigns names
  
  goData <- new("topGOdata",
                description = "aDescription", ## Tell it what we are analyzing
                ontology = ontology, ## Which ontology (BP, MF, CC)?
                allGenes = genesComp, ## Tells it which genes to use
                nodeSize = 5, # prunes GO terms with < 5 genes for greater stability
                annot = annFUN.gene2GO, # The function to do the mapping
                gene2GO = goAnno) #Sets the actual annotation
  
  
  
  ## Run the GO analyis  
  goAnalysis <- runTest(goData,
                        algorithm = algorithm, ## which algorithm Should be used
                        statistic = statistic) ## Which test stat should be kept?
  
  ## Make and save a nice table
  
  # Identify significant GO terms
  # Note this currently uses un-corrected p-values
  sigGO <- names(score(goAnalysis))[score(goAnalysis) <= pValThresh] 
  length(sigGO)
  
  ## Make a nicer Table
  goInfo <- termStat(goData,sigGO)
  goInfo$pval <- score(goAnalysis)[sigGO] 
  goInfo$goDesc <- Term(GOTERM[row.names(goInfo)])
  
  
  if(plotGO){
    ## Show a plot
    showSigOfNodes(goData,score(goAnalysis),firstSigNodes =length(sigGO),useInfo='none')
  }
 
  return(goInfo)
}
