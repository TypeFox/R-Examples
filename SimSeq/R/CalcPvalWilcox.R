CalcPvalWilcox <-
function(counts, treatment, replic = NULL, sort.method, sorted = FALSE, norm.factors, exact = FALSE){
    #calculate p-value of differential depression for each gene
    #if sort.method == paired, use signed rank test
    #if sort.method == unpaired, use ranksum test
    if(sorted == FALSE) {
      sort.list <- SortData(counts, treatment, replic, sort.method, norm.factors)
      counts <- sort.list[[1]]
      treatment <- sort.list[[3]]
      norm.factors <- sort.list[[4]]
    }
    n.row <- dim(counts)[1]
    n.col <- dim(counts)[2]
    probs <- numeric(dim(counts)[1])
    if (sort.method == "paired"){
      odds <- seq(from = 1, to = n.col, by = 2)
      diff1 <- t(counts[, odds, drop = FALSE] + 1)/norm.factors[odds] 
      diff2 <- t(counts[, (odds + 1), drop = FALSE] + 1)/norm.factors[(odds + 1)]
      diff <- log(t(diff1/diff2))
      probs <- apply(diff, 1 , function(x) wilcox.test(x, exact = exact)$p.value)        
    }
    if (sort.method == "unpaired"){
      seq.trt1 <- 1:table(treatment)[1]
      trt1 <- log(t(t(counts[, seq.trt1, drop = FALSE] + 1)/norm.factors[seq.trt1]))
      trt2 <- log(t(t(counts[, -seq.trt1, drop = FALSE] + 1)/norm.factors[-seq.trt1]))
      probs <- numeric(n.row)
      for (jj in 1:n.row) probs[jj] <-
        wilcox.test(trt1[jj, ], trt2[jj, ], exact = exact)$p.value
    }
    return(probs)
  }
