#' Determine significance of genes with multiple \emph{de novos}
#'
#' Are there more genes containing >1 \emph{de novos} than expected?
#'
#' See vignette (denovostats_intro) for more information.
#'
#' @inheritParams denovolyze
#'
#' @param nperms Number of permutations
#' @param nVars Select whether expected number of multihits is determined
#'   by "expected" total number of variants , or "actual" total.  Actual
#'   (default) is more conservative.
#'
#' @return Returns a data.frame
#'
#' @keywords keywords
#'
#' @export
#'
#' @examples
#' denovolyzeMultiHits(genes=autismDeNovos$gene,
#'                     classes=autismDeNovos$class,
#'                     nsamples=1078)
#'

##### ----------------------------------------

denovolyzeMultiHits <- function(genes,classes,nsamples,
                                nperms=100,
                                includeGenes="all",
                                includeClasses=c("syn","mis","lof","prot","all"),
                                nVars="actual",
                                geneId="geneName",
                                probTable=NULL,
                                misD=NULL,
                                signifP=3,
                                roundExpected=1) {

  # 2 options: the simulation draws N DNMs from the gene list.
  # N could be the actual number of variants seen in the population (case or control), or the expected number (based on DNM model).
  # The former is more conservative.  Samocha et al used the latter.
  # Set nVars="actual" or "expected

  if(is.null(probTable)){probTable <- pDNM}
  if(!is.null(misD)){names(probTable)[names(probTable)==misD] <- "misD"}

  # Record all variant classes in probability table.  Only calculate stats for annotated variant classes.
  allVariantClasses <- probTable %>% select(class) %>% unique %>% unlist

  # Use specified gene ID
  names(probTable)[names(probTable)==geneId] <- "gene"
  probTable$gene <- toupper(as.character(probTable$gene))
  includeGenes <- toupper(as.character(includeGenes))

  # If a list of genes for inclusion is specified, restrict analysis to these genes
  if(includeGenes[1]!="ALL"){
    probTable <- probTable[probTable$gene %in% includeGenes,]
    excludedgenes <- sum(!genes %in% includeGenes)
    if(excludedgenes > 0) {
      warning("De novo list includes ",excludedgenes," genes not specified for inclusion. These will not be analysed.")
    }
    classes <- classes[genes %in% includeGenes]
    genes <- genes[genes %in% includeGenes]
  }

  ## define internal function to perform permutation --------------------

  doPermute <- function(class,classgroup=class){
    nextvars <- genes[classes %in% classgroup]
    if(nVars == "actual") {
      x=length(nextvars)
    } else if(nVars == "expected") {
      x=2*sum(probTable$value[probTable$class==class])*nsamples
    }
    y=length(unique(nextvars[duplicated(nextvars)]))
    output <- PermuteMultiHits(x,y,nperms=nperms,class=class,geneId=geneId,includeGenes=includeGenes,probTable=probTable)
    #rownames(output) <- class
    return(output)
  }

  # --------------------
  ### Calculate probabilities for non-overlapping classes represented in data
  output <- list()
  myclasses=c("syn","misD","startloss",
              "stoploss","non","splice","frameshift")
  myclasses <- myclasses[myclasses %in% classes]
  myclasses <- myclasses[myclasses %in% allVariantClasses]
  for (class in myclasses){
    output[[class]] <- doPermute(class)
  }

  # --------------------
  ### Calculate probabilities for "aggregate classes"
  output[["mis"]] <- doPermute(class="mis",classgroup=c("mis","misD"))

  output[["lof"]] <- doPermute(class="lof",classgroup=c("non","splice","frameshift","startloss","stoploss","lof"))

  output[["prot"]] <- doPermute(class="prot",classgroup=c("non","splice","frameshift","startloss","stoploss","lof",
                                                          "mis","misD"))
  if("protD" %in% allVariantClasses){
    output[["protD"]] <- doPermute(class="protD",classgroup=c("non","splice","frameshift","startloss","stoploss","lof",
                                                              "misD"))
    }
  output[["all"]] <- doPermute(class="all",classgroup=c("non","splice","frameshift","startloss","stoploss","lof",
                                                        "mis","misD","syn"))

  # --------------------
  output <- output[names(output) %in% includeClasses] %>%
    data.frame %>%
    t %>%
    data.frame
  output$expMean <- round(output$expMean,roundExpected)
  output$pValue <- signif(output$pValue,signifP)
  return(output)
}
