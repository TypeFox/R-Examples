###############################################################################
# function RL4 with randomness control via runs test
# 
# Author: dlabes
###############################################################################
# does not require(tseries) any longer
# function runs.pvalue implemented by our self

# function for (block) randomization
RL4 <- function(nsubj, seqs=c("TR","RT"), blocksize, seed=runif(1,max=1E7), 
                randctrl=TRUE, pmethod=c("normal","exact","cc"), alpha=0.025)
{
  pmethod <- match.arg(pmethod)
  
  if (length(nsubj)>1) {
    subj <- nsubj
    nsubj <- length(subj)
  } else {
    subj <- 1:nsubj
  }
  # if we use as.numeric here the re-use of the printed seed 
  # will give different random lists! (if seed comes from Sys.time())
  seed <- as.integer(seed)
  if(is.na(seed)) seed <- 0
  set.seed(seed)

  nseq <- length(seqs)
  
  if(missing(blocksize)) blocksize <- 2*nseq
    
  if(blocksize[1]>nsubj) {
    blocksize <- 0
    warning("Blocksize > # of subjects!", 
        " Blocksize adapted to ", nsubj,".", call. = FALSE )
  }
  
  if (blocksize==0) {
    blocksize <- nsubj
  } else {
    # silent adaption or warning?
    blocksize[blocksize<nseq] <- nseq 
    # is blocksize a multiple of sequences?
    # if not, a highly unbalanced design may be the result
    # therefore we adapt the blocksize to a multiple if # of seqs
    if (any(nseq*(blocksize%/%nseq) != blocksize)) {
      blocksize <- nseq*(blocksize%/%nseq)
      # may be that due to truncation some elements are <nseq
      blocksize[blocksize<nseq] <- nseq # vector form
      warning("Blocksize is not a multiple of sequences!", 
              " Blocksize adapted to ", blocksize,".", call. = FALSE )
      blocksize <- unique(blocksize)
    }
  }  
  # call the helper function
  rlv <- rlv(nsubj, nseq, blocksize)
  rl  <- rlv$rl
  # runs test of randomnes is only possible if 2 sequences
  # if more than 2 sequences we use the dichotomization by the median
  # see wikipedia entry or package lawstat
  # -----------------
  # the simple normal approximation seems the 'best' test to reject randomness
  runs.p <- runs.pvalue(rl, pmethod=pmethod)
  maxrlen <- max(rle(rl)$lengths)
  # randomness ctrl doesn't work good in case of # of seqs>2
  # thus we test recurrent patterns
  # blocks of nseq length here implemented as matrix columns
  # gives a warning if not balanced, last column padded
  # omit the unbalanced part?
  xm <- matrix(rl[1:(nseq*(nsubj%/%nseq))], nrow=nseq)
  allblocksequal <- ncol(unique(xm,MARGIN=2))
  # randomness control
  if (randctrl){
    iter <- 0
    while(runs.p < alpha | allblocksequal==1){
      msg <- paste("runs.p= ", format(runs.p, digits=4),
                   ". Recreating randomlist.", sep="")
      if (allblocksequal==1) msg <- 
            "All nseq blocks equal. Recreating randomlist."        
      message(msg)
      # the seed must be adapted???
      seed <- seed + as.integer(runif(1,max=100))
      set.seed(seed)
      iter <- iter + 1
      rlv  <- rlv(nsubj, nseq, blocksize)
      rl   <- rlv$rl
      runs.p  <- runs.pvalue(rl, pmethod=pmethod)
      maxrlen <- max(rle(rl)$lengths)
      if (iter>=2) break
    }
  } else {
    # create a warning
    if (allblocksequal==1){
      warning("Recurrent pattern of sequences detected.")
    }
  }
  bsv <- rlv$bsv
  # create the randomlist with character representation of seqs
  rlc <- seqs[rl]
  # check if design is balanced
  ns <- table(rlc)
  nsequal <- (ns - ns[1]) == 0
  if (!all(nsequal)){
    msg  <- paste(" ", names(ns))
    msg2 <- paste(" ", ns)
    warning("Unbalanced design!", " # of subj. in sequences", 
             msg, ":", msg2, call. = FALSE)
  }
  # the random list itself
  rl <- data.frame(subject=subj, seqno=rl, sequence=rlc, stringsAsFactors=FALSE)
  # number of subjects in groups
  ns    <- t(as.matrix(ns))
  nsv   <- as.vector(ns)
  names(nsv) <- colnames(ns)
  # blocksize(s)
  if (length(unique(bsv))==1) bsv <- bsv[1]
  rlret <- list(rl=rl, seed=seed, blocksize=bsv, 
                ninseqs=nsv, runs.pvalue=runs.p, date=Sys.time())
  class(rlret) <- "rl4"
  return(rlret)
}

# internal function for (block) randomisation - working horse
# returns a list with the components
#   rl = random list (sequences numeric 1:nseq)
#   bsv = blocksize vector (actual)
rlv <- function(nsubj, nseq, blocksize)
{
  # we are working with the numeric coding of sequences
  seqsn <- 1:nseq
  n <- 0
  rl  <- vector(mode="numeric") # numeric random sequence 
  bsv <- vector(mode="numeric") # blocksize vector
  # loop over blocks
  while (n<nsubj){
    # choose a blocksize by random
    # if blocksize has only one element the sample function
    # chooses integers up to blocksize
    bs <- ifelse(length(blocksize)>1, sample(blocksize,1), blocksize[1])
    # last block may be smaller
    bs <- ifelse((nsubj-n) < bs, nsubj-n, bs)
    # random list of a block
    rpp <- ifelse(nseq*bs%/%nseq != bs, bs%/%nseq+1, bs%/%nseq)
    rlb <- sample(rep(seqsn, rpp), bs)
    # debug prints
    #cat("bs:");print(bs)
    #cat("rlb");print(rlb)
    rl <- c(rl, rlb)
    bsv <- c(bsv, bs)
    n <- n + bs
  }
  rl <- rl[1:nsubj]
  return(list(rl=rl, bsv=bsv))
}

