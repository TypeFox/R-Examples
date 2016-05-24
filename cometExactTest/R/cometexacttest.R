library(dplyr)

# tbl ordering as LSB, e.g. 2^3 contingency table, the ordering is 
# x000, x001, x010, x011, x100, x101, x110, x111. 
# For example, x101 represents the number of samples with mutations 
# in the first gene and the third gene
comet_exact_test <- function(tbl, pvalthresh = 1.1, mutmatplot = T) {
	N <- as.integer(sum(tbl))
	k <- as.integer(log2(length(tbl)))
	stopifnot(is.numeric(tbl), 2**k == length(tbl), is.integer(k), is.integer(N), is.numeric(pvalthresh))

	out_pval <- call_comet_exact_test(k, N, tbl, pvalthresh)
	#.Call('cometexacttest', k, N, as.integer(tbl), pvalthresh)
	if (mutmatplot)
	  mutmat(tbl, out_pval, k, N)
  out_pval
}

call_comet_exact_test <- function(k, N, tbl, pvalthresh){
  .Call('cometexacttest', k, N, as.integer(tbl), pvalthresh)
}

# turn input contingency table as mutation matrix
tbl2mutmat <- function(tbl, k, N){
  
  alts <- matrix(0, ncol = N, nrow = k)
  sample_index <- 1
  for (i in 0:(2^k-1) ) {
    index <- i + 1
    for ( j in 0:(tbl[index]-1)){
      col_index <- sample_index + j
      alts[,col_index] <-  tail(rev(as.integer(intToBits(i))) , k)
    }
    sample_index <- sample_index + tbl[index]
  }
  
  coverage <- sum(tbl[1:length(tbl)])
  M <- data.frame(alts)
  rownames(M) <- sapply(1:k, function(x){paste0("gene", toString(x))})
  alts <- memoSort(M)
  return(list(M=M, alts=alts, coverage=coverage))
}

# print mutation matrix
mutmat <- function(tbl, cp, k, N) {
  freq <- NULL
  key <- NULL
  gene <- NULL
  . <- NULL
  
  nM <- tbl2mutmat(tbl, k, N)
  # retrive all numbers and matrices
  M <- nM$M
  alts <- nM$alts
  ngenes <- k
  nsamples <- N
  coverage <- nM$coverage
  
  ### OncoPrint
  numOfOncos <- ngenes*nsamples;
  oncoCords <- matrix( rep(0, numOfOncos * 6), nrow=numOfOncos );
  colnames(oncoCords) <- c("xleft", "ybottom", "xright", 
                           "ytop", "altered", "cooccurrence");
  
  # create rows for indicating co-occurrence or not
  tmp_alts <- as.data.frame(alts)
  sumColumnsCooc <- tmp_alts %>% summarise_each(funs(sum))
  coverage <- sum(sumColumnsCooc > 0)
  
  xpadding <- .01;
  ypadding <- .01;
  cnt <- 1;
  for(i in 1:ngenes) {
    for(j in 1:nsamples) {
      xleft <- j-1 + xpadding;
      ybottom <- ((ngenes-i+1) -1) + ypadding;
      xright <- j - xpadding;
      ytop <- (ngenes-i+1) -ypadding;
      cooc <- sumColumnsCooc[[j]]
      
      altered <- alts[i, j]
      
      oncoCords[cnt, ] <- c(xleft, ybottom, xright, ytop, altered, cooc);
      cnt <- cnt+1;
    }
  }
  
  agene <- rownames(M)
  
  tmp_M <- as.data.frame(M)
  alts_freq <- bind_cols(
    (tmp_M %>%
       mutate(freq=rowSums(.)) %>% select(freq)),
    (tmp_M %>% mutate(gene=agene) %>% select(gene) )
  ) %>% mutate(key=paste0(gene,"[",freq,"]")) %>% select(key) %>% unlist()
  
  colors <- rep("lightgray", cnt);
  colors[ which(oncoCords[, "altered"] == 1 & 
                  oncoCords[, "cooccurrence"] == 1 ) ] <- "steelblue1";
  colors[ which(oncoCords[, "altered"] == 1 & 
                  oncoCords[, "cooccurrence"] > 1 ) ] <- "darkorange2";
  
  pdf('cometExactTest_plot.pdf', width = 10, height = ngenes+1)
  
  plot(c(0, nsamples), c(0, ngenes), type="n", 
       main= sprintf("Gene set altered in %.2f%%: %d of %d cases. 
                      CoMEt pvalue: %s.", coverage/nsamples*100, coverage, 
                      nsamples, format(cp, scientific = T, digits = 3)), 
       xlab="Samples", ylab="", yaxt="n");
  rect(oncoCords[, "xleft"], 
       oncoCords[, "ybottom"],
       oncoCords[, "xright"], 
       oncoCords[, "ytop"], 
       col=colors, border="white");
  axis(2, at=(ngenes:1)-.5, labels=alts_freq, las=3, cex.axis = 0.5);
  dev.off()
}


# This function sort the matrix to better visualize Mutual Exclusivity
# source: https://github.com/dakl/oncoprint/blob/master/R/memosort.R
memoSort <- function(M) {
  
  scoreCol <- function(x) {
    score <- 0;
    for(i in 1:length(x)) {
      if(x[i]) {
        score <- score + 2^(length(x)-i);
      }
    }
    return(score);
  }
  
  if (nrow(M) == 1){
    scores <- M[1,]
    sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
    return(M[1, sampleOrder])
  }
  else{
    geneOrder <- sort(rowSums(M), decreasing=TRUE, index.return=TRUE)$ix;
    scores <- apply(M[geneOrder, ], 2, scoreCol);
    sampleOrder <- sort(scores, decreasing=TRUE, index.return=TRUE)$ix;
    return(M[geneOrder, sampleOrder]);
  }
  
}
