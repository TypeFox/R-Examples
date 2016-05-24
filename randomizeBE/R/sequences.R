###############################################################################
# Function to get the vector of sequences for the crossover designs
# 
# Author: dlabes
###############################################################################

sequences <- function(design, tmts=NULL)
{ 
  seqs <- character()
  ntmt <- as.numeric(substr(design,1,1))
  if (!is.null(tmts) & length(tmts)!=ntmt){
    stop ("Treatment codes must have ", ntmt, " entries!")
  }
  
  # classical 2x2 crossover
  if(design=="2x2" | design=="2x2x2")  seqs <- c("AB","BA")
  # 2-group parallel
  if(design=="parallel")               seqs <- c("A","B")
  # 3x3 crossover
  if(design=="3x3" | design=="3x3x3"){
    seqs <- c("ABC","BCA","CAB")
    # randomize it
    seqs <- randLS(seqs)
  }
  # 3x3 crossover: Williams design with 6 sequencess
  if(design=="3x6x3"){
    seqs <- c("ABC",
              "BCA",
              "CAB",
              "ACB",
              "BAC",
              "CBA")
  }
  # 4x4 crossover
  if(design=="4x4" | design=="4x4x4"){
    # this is one of the four standard latin squares, it is a williams design
    seqs <- c("ABCD",
              "BDAC",
              "CADB",
              "DCBA")
    # eventually we should randomize it? yes
    seqs <- randLS(seqs)
    # although the standard Latin square is a williams design
    # the randLS gives not always a Williams design back
  }
  # partial replicate (reference replicate)
  if(design=="2x3x3") seqs <- c("ABB","BAB","BBA")
  # 2-sequence-3-period full replicate
  if(design=="2x2x3") seqs <- c("ABA","BAB")
  # 2-sequence-4-period full replicate: FDA design
  if(design=="2x2x4") seqs <- c("ABAB","BABA")
  # 4-sequence-4-period  full replicate
  # Chen, Chow, Li: "SAMPLE SIZE higher order crossover"
  if(design=="2x4x4") seqs <- c("ABBA","BAAB","AABB","BBAA")
  # Baalams design
  if(design=="2x4x2") seqs <- c("AB","BA","AA","BB")
  
  # if tmts are given
  if (! is.null(tmts)){
    #replace the a codes by tmts
    for (i in seq_along(tmts)){
      st <- LETTERS[i]
      seqs <- gsub(st,tmts[i],seqs)
    }
  }
  if (length(seqs)==0) stop("Design not implemented!")
  return(seqs)
}

# ------------------------------------------------------------------------
# internal function to randomize latin square
randLS <- function(seqs)
{
  seqm <- strsplit(seqs, split="")
  n    <- length(seqm)
  seqm2 <- matrix(nrow=n,ncol=n)
  for (i in seq_along(seqm)){
    seqm2[i,] <- seqm[[i]]
  }
  rows <- sample(1:n)
  cols <- sample(1:n)
  # order by rows
  seqm2 <- seqm2[rows,]
  seqm2 <- seqm2[,cols]
  # sort by first col?
  seqm2 <- seqm2[order(seqm2[,1]),] 
  for (i in seq_along(seqm2[,1])){
    seqs[i] <- paste(seqm2[i,],sep="",collapse="")
  }
  seqs
}