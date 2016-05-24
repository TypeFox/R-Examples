readGenome <- function(fn.genome, ex.sh.aa=0, rm.first.aa=0)
{
  ## read sequences
  seq.data <- read.seq(fn.genome)
  
  ## remove leading aa
  if(rm.first.aa > 0)
  {
    seq.names <- names(seq.data)
    seq.data <- lapply(1:length(seq.data), function(i)
    {
      seq <- seqinr::getSequence.SeqFastadna(seq.data[[i]])[-(4:((rm.first.aa+1)*3))] # +1 is necessary since I leave the start codon and remove rm.first.aa AA after the start codon
      return(seqinr::as.SeqFastadna(seq, name=seqinr::getName.SeqFastadna(seq.data[[i]])))
    })
    names(seq.data) <- seq.names # names are lost thus set them here again
  }
  
  ## remove sequences shorter than threshold
  if(ex.sh.aa)
  {
    ind <- unlist(lapply(1:length(seq.data), function(i)
    {
      return(length(seq.data[[i]]) > (ex.sh.aa*3))
    }))
    ## eleminate short sequences (length < ex.sh.aa)
    seq.data <- seq.data[ind]  
  }
  
  ## convert sequences into proper format
  seq.string <- convert.seq.data.to.string(seq.data)
  seq.string <- seq.string[order(names(seq.string))]
  
  return(seq.string)
}

## set mean of dataset to 1
normalizeDataSet <- function(data)
{
  data <- data/mean(data)
}