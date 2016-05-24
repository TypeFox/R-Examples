#' Performs a split-apply-combine on an ffdf
#'
#' Performs a split-apply-combine on an ffdf. 
#' Splits the x ffdf according to split and applies FUN to the data, stores the result of the FUN in an ffdf.\cr
#' Remark that this function does not actually split the data. In order to reduce the number of times data is put into RAM for situations with a lot
#' of split levels, the function extracts groups of split elements which can be put into RAM according to BATCHBYTES. Please make sure your FUN covers the
#' fact that several split elements can be in one chunk of data on which FUN is applied.\cr
#' Mark also that NA's in the split are not considered as a split on which the FUN will be applied.
#'
#' @example ../examples/ffdfplyr.R
#' @param x an ffdf
#' @param split an ff vector which is part of the ffdf x
#' @param FUN the function to apply to each split. This function needs to return a data.frame
#' @param BATCHBYTES integer scalar limiting the number of bytes to be processed in one chunk
#' @param RECORDBYTES optional integer scalar representing the bytes needed to process one row of x
#' @param trace logical indicating to show on which split the function is computing
#' @param ... other parameters passed on to FUN
#' @return 
#' an ffdf
#' @export
#' @importFrom fastmatch fmatch
#' @seealso \code{\link{grouprunningcumsum}, \link{table}}
ffdfdply <- function (x, split, FUN, BATCHBYTES = getOption("ffbatchbytes"), 
                      RECORDBYTES = sum(.rambytes[vmode(x)]), trace = TRUE, ...) {
    
  MAXSIZE = BATCHBYTES/RECORDBYTES  
  force(split)
  
  if(nrow(x) != length(split)){
    stop("split needs to be the same length as the number of rows in x")
  }
  ##
  ## Detect optimal split size -> 2 runs over the split data
  ##
  if(trace) message(sprintf("%s, calculating split sizes", Sys.time()))
  if(!is.factor.ff(split)) {
    warning("split needs to be an ff factor, converting using as.character.ff to an ff factor")
    splitby <- as.character.ff(split)
  }else{
    splitby <- split
  }
  tmp <- splitby
  levels(tmp) <- NULL
  missings <- is.na(tmp)
  if(sum(missings) > 0){
    splitby <- splitby[!missings]
  }
  splitgroups <- list()
  splitgroups$tab <- binned_tabulate.ff(x = ff(factor("data", levels="data"), length = length(splitby)), bin = splitby, nbins = max(tmp), nlevels = 1)
  splitgroups$tab <- as.table(splitgroups$tab[,'data'])
  splitgroups$tab <- splitgroups$tab[order(splitgroups$tab, decreasing = TRUE)]
  if (max(splitgroups$tab) > MAXSIZE) {
    warning("single split does not fit into BATCHBYTES")
  }
  splitgroups$tab.groups <- as.integer(grouprunningcumsum(x = as.integer(splitgroups$tab), max = MAXSIZE))
  
  ##
  ## Identify split locations -> 1 run over the split data
  ##
  if(trace) message(sprintf("%s, building up split locations", Sys.time()))
  splitpositions <- list()
  splitpositions$rowidxgroups <- ffdf(pos = ffseq_len(as.integer(length(split))), split = split)
  splitpositions$rowidxgroups$group <- ffdfwith(splitpositions$rowidxgroups, splitgroups$tab.groups[fmatch(as.character(split), names(splitgroups$tab))])
  splitpositions$rowidxgroups <- splitpositions$rowidxgroups[c("pos","group")]  
  splitpositions$nrsplits <- max(splitgroups$tab.groups)
  
  splitpositions$grouprowidx <- list()  
  splitpositions$runninggrouppos <- list()
  for(idx in 1:splitpositions$nrsplits){
    splitpositions$grouprowidx[[as.character(idx)]] <- ff(as.integer(NA), length = sum(splitgroups$tab[splitgroups$tab.groups == idx]), vmode = "integer")
    splitpositions$runninggrouppos[[as.character(idx)]] <- 0
    close(splitpositions$grouprowidx[[as.character(idx)]])
  }  
  for(idx in chunk(splitpositions$rowidxgroups)){
    tmp <- splitpositions$rowidxgroups[idx, ]
    tmp <- split(tmp$pos, tmp$group)
    for(group in names(tmp)){      
      loc <- ri(from=splitpositions$runninggrouppos[[group]]+1, to = splitpositions$runninggrouppos[[group]] + length(tmp[[group]]))      
      open(splitpositions$grouprowidx[[group]])
      splitpositions$grouprowidx[[group]][loc] <- tmp[[group]]
      splitpositions$runninggrouppos[[group]] <- splitpositions$runninggrouppos[[group]] + length(tmp[[group]])
      close(splitpositions$grouprowidx[[group]])
    }
  }
  
  ##
  ## Apply the function on each group of split elements
  ##
  allresults <- NULL
  for (idx in 1:splitpositions$nrsplits) {
    fltr <- splitpositions$grouprowidx[[idx]]
    open(fltr)
    if(trace) {
      message(sprintf("%s, working on split %s/%s, extracting data in RAM of %s split elements, totalling, %s GB, while max specified data specified using BATCHBYTES is %s GB", 
                      Sys.time(), idx, splitpositions$nrsplits, length(splitgroups$tab[splitgroups$tab.groups == idx]), round(as.numeric(RECORDBYTES) * as.numeric(length(fltr))/2^30, 5), round(BATCHBYTES/2^30, 5)))
    }
    inram <- ffdfget_columnwise(x, fltr)
    close(fltr)
    if(trace) message(sprintf("%s, ... applying FUN to selected data", Sys.time()))
    result <- FUN(inram, ...)
    if (!inherits(result, "data.frame")) {
      stop("FUN needs to return a data frame")
    }
    if(trace) message(sprintf("%s, ... appending result to the output ffdf", Sys.time()))
    rownames(result) <- NULL
    if (!is.null(allresults) & nrow(result) > 0) {
      rownames(result) <- (nrow(allresults) + 1):(nrow(allresults) + nrow(result))
    }
    if (nrow(result) > 0) {
      allresults <- ffdfappend(x = allresults, dat = result, recode = FALSE)
    }
  }
  allresults
}