#' Checks input for errors
#'
#' An internal function to check inputs
#'
#' @inheritParams denovolyze
#'
#' @return warning or error if any invalid input, else assigns variables back to
#'   parent function
#'

parseInput <- function(genes=genes,
                       classes=classes,
                       nsamples=nsamples,
                       groupBy=groupBy,
                       includeGenes=includeGenes,
                       includeClasses=includeClasses,
                       geneId=geneId,
                       signifP=signifP,
                       roundExpected=roundExpected,
                       probTable=NULL){

  ## Define a function to verify integer inputs --------------------
  is.wholenumber <- function(x, tol = .Machine$double.eps^0.5) {abs(x - round(x)) < tol}

  ## Use specified probTable --------------------
  if(is.null(probTable)){probTable <- pDNM}

  ## checks on geneId --------------------
  if(!geneId %in% names(probTable)){
    stop(paste(geneId,"is not a valid geneId"))
  }
  # use geneID if valid
  names(probTable)[names(probTable)==geneId] <- "gene"
  probTable$gene <- toupper(probTable$gene)
  #pass probTable to parent function
  assign("probTable",probTable,pos=sys.frame(sys.parent()))

  ## check inputs have same length --------------------
  if(length(genes)!=length(classes)){
    stop('The number of genes (genes) and number of variant consequences (classes) do not match')
  }

  ## checks on genenames --------------------

  # character & capitalisation
  genes <- toupper(genes)

  # recognised
  noMatchGenes <- (!genes %in% probTable$gene) %>% sum
  if(noMatchGenes>0){
    warning(
      paste(noMatchGenes,"gene identifiers in input list do not match the probability table, and are excluded from analysis."
      ))
  }

  # no NA
  if(sum(is.na(genes))>0){
    stop('gene name can not be NA')
  }

  #pass capitalised gene names to parent function
  assign("genes",genes,pos=sys.frame(sys.parent()))


  ## checks on variant classes --------------------
  # capitalisation - not implemented
  # as.character
  classes <- as.character(classes)
  # recognised
  validClasses <- unique(probTable$class)
  noMatchClasses <- classes[!classes %in% validClasses] %>% unique
  if(length(noMatchClasses)>0){
    stop(paste("The following are not recognised variant classes:",noMatchClasses))
  }
  # no NA
  if(sum(is.na(classes))>0){
    stop('variant class can not be NA')
  }
  # match to SO? - not implemented

  assign("classes",classes,pos=sys.frame(sys.parent()))

  ## checks on nsamples --------------------
  # is integer
  if(!is.wholenumber(nsamples)){
    stop('nsamples must be an integer')
  }

  ## checks on groupBy --------------------
  groupBy <- tolower(groupBy)
  if(!groupBy %in% c("gene","class")){
    stop(paste("\"",groupBy,"\" is not a valid groupBy option",sep=""))
  }
  # passes lower case groupBy to parent function
  assign("groupBy",groupBy,pos=sys.frame(sys.parent()))


  ## checks on includeGenes --------------------
  includeGenes <- toupper(includeGenes)
  #apply only if includeGenes is not "all"
  if(!includeGenes[1]=="ALL" & length(includeGenes==1)){
    noMatchGenes <- (!includeGenes %in% probTable$gene) %>% sum
    if(noMatchGenes>0){
      warning(
        paste(noMatchGenes,"gene identifiers in \"includeGene\" are not in the probability table, and are excluded from analysis."
        ))
    }
  }
  assign("includeGenes",includeGenes,pos=sys.frame(sys.parent()))

  ## checks on includeClasses --------------------
  validClasses  <- c("syn","mis","misD",
                                    "non","stoploss","startgain",
                                    "splice","frameshift","lof","prot",
                                    "protD", "all")
  noMatchClasses <- includeClasses[!includeClasses %in% validClasses] %>% unique
  if(length(noMatchClasses)>0){
    stop(paste("The following are not recognised variant classes:",noMatchClasses))
  }

  ## checks on roundExpected --------------------
  # is integer
  if(!is.wholenumber(roundExpected)){
    stop('roundExpected must be an integer')
  }

  ## checks on signifP --------------------
  # is integer
  if(!is.wholenumber(signifP)){
    stop('signifP must be an integer')
  }


} #end of function

##### --------------------

# TESTS
# genes=autismDeNovos$gene;classes=autismDeNovos$class;nsamples=1078
# denovolyze(genes,classes,nsamples)
# denovolyze(genes,classes[-1],nsamples)
# denovolyze(genes=c("KCNQ1","RYR2",NA), classes=c("lof","lof","lof"), nsamples=10)
# denovolyze(genes=c("KCNQ1","RYR2","RYR2"), classes=c("lof","lof","invalidClass"), nsamples=10)
# denovolyze(genes,classes,nsamples=1078.3)
# denovolyze(genes,classes,nsamples,groupBy="magic")
# denovolyze(genes,classes,nsamples,geneId="invalidId")
# denovolyze(genes,classes,nsamples,roundExpected=0.6)
# denovolyze(genes,classes,nsamples,signifP=0.6)
# denovolyze(genes,classes,nsamples,includeGenes=c("MARK1","misspelledGene"))
# denovolyze(genes,classes,nsamples,includeClasses="madeUpClass")
