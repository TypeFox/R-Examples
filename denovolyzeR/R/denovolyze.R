#' Evaluates burden of \emph{de novo} variation against expectation
#'
#' Determines whether the test population carry more \emph{de novo} variants
#' than expected. Variants may be grouped by variant class (e.g. are there more
#' LOF variants than expected, across the whole dataset?), or by gene (are there
#' more variants of a given class in SCN2A?).
#'
#' Analyses can be restricted to a subset of genes, and/or a subset of variant
#' classes
#'
#' See vignette("denovolyzeR_intro") for more information.
#'
#' @param genes A vector of genes containing de novo variants.
#' @param classes A vector of classes of de novo variants.  Standard supported
#'   classes are "syn" (synonymous), "mis" (missense), "non" (nonsense),
#'   "splice" (splice), "frameshift" (frameshift) and "lof" (loss of function =
#'   non + splice + frameshift). Additional classes that are supported by the
#'   code, but are not included in the built-in probability tables, are
#'   "stoploss","startloss", "misD" (damaging missense).  These labels may be
#'   used for user-supplied probability tables. If "misD" is present, then "mis"
#'   (in the input) implies non-damaging missense.
#' @param nsamples Number of individuals considered in de novo analysis.
#' @param groupBy Results can be tabulated by "gene", or by variant "class"
#' @param includeClasses Determines which variant classes are tabulated in
#'   output.  In addition to the input classes, summaries can be produced for
#'   "prot" (protein-altering = mis + lof), "all", and "protD" (protein damaging
#'   = misD + lof, only available if misD included in user-specified probability
#'   table).  If "misD" is present, then "mis" will return statistics for all
#'   missense.  Non-damaging missense are not analysed separately.
#' @param includeGenes Genes to include in analysis. "all" or a vector of gene
#'   names.
#' @param geneId Gene identifier used. One of "hgncID", "hgncSymbol",
#'   "enstID", "ensgID" or "geneName" (default, equals ensembl "external_gene_name")
#' @param signifP Number of significant figures used to round p-values in
#'   output.
#' @param roundExpected Number of decimal places used to round expected burdens
#'   in output.
#' @param probTable Probability table. A user-defined table of probabilities can
#'   be provided here, to replace the probability table included in the package.
#' @param misD If the user-specified probability table contains probabilities
#'   for a sub-category of missense variants (e.g. predicted to be damaging by
#'   an in silico algorithm), this column should be called misD, or the
#'   alternative name should be specified here.
#'
#'
#' @return Returns a data frame
#'
#' @export denovolyze denovolyzeByClass denovolyzeByGene
#'
#' @import dplyr
#'
#' @examples
#'
#' ### denovolyze
#'
#' denovolyze(genes=autismDeNovos$gene,
#'            classes=autismDeNovos$class,
#'            nsamples=1078)
#'
#' ### denovolyzeByClass
#'
#' denovolyzeByClass(genes=autismDeNovos$gene,
#'                   classes=autismDeNovos$class,
#'                   nsamples=1078)
#'
#' # this convenience function is identical to:
#'
#' denovolyze(genes=autismDeNovos$gene,
#'            classes=autismDeNovos$class,
#'            nsamples=1078,
#'            groupBy="class",
#'            includeClasses=c("syn","mis","lof","prot","all"),
#'            includeGenes="all"
#'            )
#'
#' ### denovolyzeByGene
#'
#' denovolyzeByGene(genes=autismDeNovos$gene,
#'                  classes=autismDeNovos$class,
#'                  nsamples=1078)
#'
#' # this is identical to:
#'
#' denovolyze(genes=autismDeNovos$gene,
#'            classes=autismDeNovos$class,
#'            nsamples=1078,
#'            groupBy="gene",
#'            includeClasses=c("lof","prot"),
#'            includeGenes="all"
#'            )
#'

##### ----------------------------------------------------

denovolyze <- function(genes,classes,nsamples,
                       groupBy="class",
                       includeGenes="all",
                       includeClasses=c("syn","mis","misD",
                                        "non","stoploss","startgain",
                                        "splice","frameshift","lof","prot",
                                        "protD", "all"),
                       geneId="geneName",
                       signifP=3,
                       roundExpected=1,
                       probTable=NULL,
                       misD=NULL) {

  # this line defines variables in order to pass R CMD check
  # these are column names used in dplyr::select(x) statement, but R CMD CHECK interprets them as global variables without visible binding
  gene <- value <- enrichment <- Row.names <- ends_with <- NULL

  # check inputs --------------------------
  parseInput(genes,
             classes,
             nsamples,
             groupBy,
             includeGenes,
             includeClasses,
             geneId,
             signifP,
             roundExpected,
             probTable)

  # With an external probTable table, label the damaging missense column as misD
  # --------------------------
  if(!is.null(misD)){names(probTable)[names(probTable)==misD] <- "misD"}

  #  # Include all genes if indicated --------------------------
  #  # If byClass - include all genes.  if byGene - include only genes containing at least one variant
  if(toupper(includeGenes[1])=="ALL" & length(includeGenes==1)){
    if(groupBy=="class"){
      includeGenes <- toupper(probTable$gene)
    } else if (groupBy=="gene"){
      includeGenes <- unique(genes)
    }

  }

  # generate meta-classes: "lof", "prot", "protD" and "all".  if "misD" is present, "mis" in input = not-damaging mis.  In output mis will refer to all missense.
  # --------------------------
  input <- data.frame(gene=genes,class=classes,stringsAsFactors=F)
  input$class.1[input$class %in% c("splice","frameshift","non","stoploss","startloss")] <- "lof"
  #temporarily relabel non-damaging missense as mis_notFilter.
  input$class[input$class == "mis"] <- "mis_notFilter"
  input$class.2[input$class %in% c("mis_notFilter","misD")] <- "mis"
  input$class.3[input$class %in% c("splice","frameshift","non","stoploss","startloss",
                                   "lof",
                                   "mis_notFilter","misD")] <- "prot"
  input$class.4[input$class %in% c("splice","frameshift","non","stoploss","startloss",
                                   "lof","misD")] <- "protD"
  input$class.5 <- "all"

  input <- reshape::melt.data.frame(input,id.vars="gene") %>%
    dplyr::select(gene, class = value)  %>%
    filter(!is.na(class)) %>%
    filter(class!="mis_notFilter")

  # tabulate observed & expected numbers, either by gene or by class
  # --------------------------
  if(groupBy=="class"){
    observed <-
      input %>%
      filter(gene %in% includeGenes, class %in% includeClasses) %>%
      group_by(class) %>%
      summarise(
        observed = n()
      )
    observed$class <- factor(observed$class, levels=c(c("syn","misD","mis",
                                                        "non","stoploss","startloss",
                                                        "splice","frameshift","lof","prot","protD","all")))
    observed <- observed[order(observed$class),]


    expected <- probTable %>%
      filter(gene %in% includeGenes, class %in% includeClasses) %>%
      group_by(class) %>%
      summarise(
        expected = 2*sum(value, na.rm=T)*nsamples
      )
    expected$class <- factor(expected$class, levels=c(c("syn","misD","mis",
                                                        "non","stoploss","startloss",
                                                        "splice","frameshift","lof","prot","protD","all")))


    output <- left_join(observed,expected,by=c("class"))

  } else if(groupBy=="gene"){

    obs <- input %>%
      filter(gene %in% includeGenes, class %in% includeClasses) %>%
      group_by(gene,class) %>%
      summarise(
        obs = n()
      )

    exp <- probTable %>%
      filter(gene %in% includeGenes, class %in% includeClasses) %>%
      group_by(gene,class) %>%
      summarise(
        exp = 2*sum(value, na.rm=T)*nsamples
      )

    output <- merge(obs,exp,by=c("gene","class"),all=T)
  }

  # calculate poisson stats and enrichments ____________________
  output[is.na(output)] <- 0
  output$enrichment <- signif(output$obs/output$exp,signifP)
  output$pValue <- signif(ppois(output$obs-1,lambda=output$exp,lower.tail=F),signifP)
  if("exp" %in% names(output)){
    output$exp <- round(output$exp,roundExpected)
  } else if ("expected" %in% names(output)){
    output$expected <- round(output$expected,roundExpected)
  }

  # when analysing by gene, apply additional formatting to arrange results for different variant classes side by side
  #  --------------------------
  if(groupBy=="gene"){

    #For single gene, format 1 row per class.  Else 1 row per gene.
    if(length(includeGenes)==1){
      output <-  dplyr::select(output,-enrichment)
    } else {
      output <- output %>%
        dplyr::select(-enrichment) %>%
        reshape::recast(id.var=c("gene","class"), formula = gene ~ variable ~ class)

      classNotRepresented <- includeClasses[!includeClasses %in% dimnames(output)$class]
      if (length(classNotRepresented)!=0){
        warning(paste(classNotRepresented,"is not found in data"))
        includeClasses[includeClasses %in% unique(output$class)]
      }

      output2 <- list()
      for (i in seq(along=includeClasses)){
        output2[[i]] <- as.data.frame(output[,,i])
        names(output2[[i]]) <- paste(includeClasses[i],names(output2[[i]]),sep=".")
      }

      output3 <- output2[[1]]
      for (i in seq(along=includeClasses)[-1]){
        output3 <- merge(output3,output2[[i]],by="row.names",all=T)
        rownames(output3) <- output3$Row.names
        output3 <- dplyr::select(output3,-Row.names)
      }

      myIndex <- output3 %>% dplyr::select(ends_with("pValue")) %>% apply(MARGIN=1,min, na.rm=T) %>% order()

      output <- output3[myIndex,]
    }
  }
  # --------------------------
  class(output) <- "data.frame"
  return(output)
}

##### ----------------------------------------------------

#' @describeIn denovolyze denovolyzeByClass

denovolyzeByClass <- function(genes,classes,nsamples,
                              groupBy="class",
                              includeGenes="all",
                              includeClasses=c("syn","mis","lof","prot","all"),
                              geneId="geneName",
                              signifP=3,roundExpected=1,
                              probTable=NULL){
  denovolyze(genes,classes,nsamples,groupBy,includeGenes,includeClasses,geneId,signifP,roundExpected,probTable)
}

##### ----------------------------------------------------

#' @describeIn denovolyze denovolyzeByGene

denovolyzeByGene <- function(genes,classes,nsamples,
                             groupBy="gene",
                             includeGenes="all",
                             includeClasses=c("lof","prot"),
                             geneId="geneName",
                             signifP=3,roundExpected=1,
                             probTable=NULL){
  denovolyze(genes,classes,nsamples,groupBy,includeGenes,includeClasses,geneId,signifP,roundExpected,probTable)
}


