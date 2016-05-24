SimData <-
function(counts, treatment, replic = NULL, sort.method, 
                     k.ind, n.genes = NULL, n.diff = NULL, 
                     norm.factors = NULL, samp.independent = FALSE,
                     genes.select = NULL, genes.diff = NULL, switch.trt = FALSE, 
                     probs = NULL, weights = NULL, exact = FALSE, power = 1){
  #given matrix of gene expression counts with two treatment groups, 
  #simulates matrix of gene expression data with known list of DE
  #and EE genes. See help file for more info.
  
  #check inputs
  counts <- as.matrix(counts)
  if (any(!is.numeric(counts))) 
    stop("Error: counts matrix contains non-numeric values.")
  if (any(round(counts) != counts)) 
    stop("Error: counts matrix contains non-integer values.")

  n.row <- dim(counts)[1]
  n.col <- dim(counts)[2]
  
  genes1.diff <- NULL
  if (is.null(n.genes) && is.null(genes.select)) 
    stop("Error: Both n.genes and genes.select are set to NULL.")
  
  if (!is.null(genes.select)){
    if (!is.numeric(genes.select) && !is.logical(genes.select)) 
      stop("Error: genes.select must be a NULL, numeric or logical vector.")
    if (is.logical(genes.select)){
      if (length(genes.select) != n.row) 
        stop("Error: When genes.select is a logical vector, 
             its length must equal the number of rows of the counts matrix.") 
      else {
        genes.select <- which(genes.select)
      }
    }
    if (!is.null(n.genes)){
      if (length(genes.select) != n.genes) 
        stop("Error: n.genes must equal number of genes selected in 
                 simulation through genes.select vector.")
    }
    if (is.null(n.genes)) n.genes <- length(genes.select) 
  }
  
  if (is.null(genes.select) && !is.null(genes.diff)) 
    stop("Error: genes.diff vector cannot be specified without 
         first specifying genes.select vector.")
  
  if (is.null(n.diff) && is.null(genes.diff)) 
    stop("Error: Both n.diff and genes.diff are set to NULL.")
  
  if (!is.null(genes.select) && !is.null(genes.diff)){
    if (!is.numeric(genes.diff) && !is.logical(genes.diff)) 
      stop("Error: genes.diff must be a NULL, numeric or logical vector.")
    if (is.numeric(genes.diff)){
      if (any(!genes.diff %in% genes.select)) 
        stop("Error: genes.diff vector must be a subset of genes.select.")
    } 
    if (is.logical(genes.diff)){
      if (length(genes.diff) != length(genes.select)) 
        stop("Error: When genes.diff is a logical vector, length of genes.diff must equal 
             number of genes selected in genes.select vector.")
      else {
        genes.diff <- genes.select[genes.diff]
      }
    }
    if (!is.null(n.diff)){
      if (length(genes.diff) != n.diff) 
        stop("Error: n.diff must equal number of genes to be differentially expressed
              as selected through genes.diff vector.")
    }
    if (is.null(n.diff)) n.diff <- length(genes.diff)
  }
  
  if (!is.numeric(n.genes)) 
    stop("Error: Number of genes selected, n.genes, must be a positive integer 
          less than the total number of rows in the counts matrix.")
  else if (round(n.genes) != n.genes || n.genes > n.row || n.genes <= 0) 
    stop("Error: Number of genes selected, n.genes, must be a positive integer 
         less than the total number of rows in the counts matrix.")

  if (!is.numeric(n.diff)) 
    stop("Error: Number of DE genes, n.diff, must be a positive integer 
          less than n.genes.")
  else if (round(n.diff) != n.diff || n.diff > n.genes || n.diff < 0) 
    stop("Error: Number of DE genes, n.diff, must be a positive integer 
          less than n.genes.")
  
  if (!is.numeric(k.ind)) 
    stop("Error: Number of replicates simulated per treatment group, k.ind, 
          must be a positive integer less than total the number of 
          columns in the counts matrix divided by two.")
  else if (round(k.ind) != k.ind || k.ind > n.col/2 || k.ind <= 0) 
    stop("Error: Number of replicates in each treatment group, k.ind, 
          must be a positive integerless than the total number of 
          columns in counts matrix divided by two.")

  if (!is.null(probs)){
    if (!is.numeric(probs))
      stop("Error: probs must be a numeric vector of length equal to the 
           number of rows in the counts matrix with each entry between 0 and 1.") 
    if (any(probs > 1) || any(probs < 0) || length(probs) != n.row) 
      stop("Error: probs must be a numeric vector of length equal to the 
           number of rows in the counts matrix with each entry between 0 and 1.") 
  }
  
  if (!is.null(weights)){
    if (!is.numeric(weights))
      stop("Error: weights must be a numeric vector of length equal to the 
           number of rows of the counts matrix.")
    if (length(weights) != n.row) 
      stop("Error: weights must be a numeric vector of length equal to the 
           number of rows of the counts matrix.")
  }
  
  if (!is.null(replic)){
    if(length(replic) != n.col) 
      stop("Error: Length of replic vector must equal number of 
           columns in counts matrix")
    if (any(tabulate(replic) > 2)) 
      stop("Error: Number of observations per replicate in counts 
           matrix must not be greater than 2.")
  }
  
  if(length(treatment) != n.col)
    stop("Error: Length of treatment vector must equal number of 
         columns in counts matrix")
  if (length(unique(treatment)) != 2) 
    stop("Error: Number of treatment groups in counts matrix 
         must be equal to two.")
  if (!sort.method %in% c("paired", "unpaired")) 
    stop("Error: sort.method must be set to either 'paired' or 'unpaired'.")
  if (sort.method == "paired" && is.null(replic)) 
    stop("Error: Must specify replic vector when sort.method equals 'paired'.")
  
  if(!is.null(norm.factors)){
    if (!is.numeric(norm.factors)) 
      stop("Error: norm.factors must be a positive numeric vector with 
           length equal to the number of columns in the counts matrix.")
    if (is.numeric(norm.factors)){
      if (any(norm.factors <= 0) || length(norm.factors) != n.col) 
        stop("Error: norm.factors must be a positive numeric vector with 
             length equal to the number of columns in the counts matrix.")
      }
    }
  
  if (is.null(norm.factors)) {
    norm.factors <- apply(counts, 2, quantile, 0.75)
    norm.factors <- norm.factors/exp(mean(log(norm.factors))) #normalize so that geometric mean is 1
  }
  #sort necessary inputs
  sort.list <- SortData(counts = counts, treatment = treatment,
                        replic = replic, sort.method = sort.method, norm.factors = norm.factors)
  counts <- sort.list[[1]]
  treatment <- sort.list[[3]]
  norm.factors <- sort.list[[4]]
  sorting <- sort.list[[5]]
  
  #reset n.row and n.col variables since SortData function may have trimmed
  #counts matrix
  n.row <- dim(counts)[1]
  n.col <- dim(counts)[2]
  
  #calculate p-value of differential depression for each gene
  if (is.null(probs) && is.null(weights) && is.null(genes.diff)){
    probs <- CalcPvalWilcox(counts = counts, treatment = treatment,
                            sort.method = sort.method, sorted = TRUE, norm.factors = norm.factors)
  }  
  
  #calulate local fdr (empirical bayes probability of differential expression) for each gene
  if (is.null(weights)){
    if (is.null(genes.select) | is.null(genes.diff)) weights <- 
      1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr
  } 
  
  SampGenes <- function(genes.select, genes.diff, n.genes, n.diff, 
                         weights, power, exact){
    #sample genes to be used in simulation and sample which genes are
    #to be differentially expressed
    
    genes <- 1:n.row
    if (is.null(genes.diff)){
      if(n.diff > 0){
        if(is.null(genes.select)){
          genes.diff <- 
            sample(genes, n.diff, prob = (weights)^power, replace = FALSE)
        } else {
          genes.diff <- 
            sample(genes.select, n.diff, prob = (weights[genes.select])^power, replace = FALSE)
        }
      }
      else genes.diff <- NULL
    }     
    if (is.null(genes.select)) {
      if(is.null(genes.diff)) genes.subset <- sample(genes, n.genes, replace = FALSE)
        else {
          genes.subset <-
        c(sample(genes[-genes.diff], n.genes - n.diff, replace = FALSE), genes.diff)
      }
    }
    else genes.subset <- genes.select

    genes.subset <- sort(genes.subset)
    if (!is.null(genes.diff)) genes.diff <- sort(genes.diff)
    DE.genes <- genes.subset %in% genes.diff
    
    return(list(genes.subset = genes.subset, genes.diff = genes.diff,  
                DE.genes = DE.genes))
  }
  
  #sample EE and DE genes
  samp.genes.list <- SampGenes(genes.select, genes.diff, n.genes, n.diff,
                                weights, power, exact)

  genes.subset <- samp.genes.list$genes.subset
  genes.diff <- samp.genes.list$genes.diff
  DE.genes <- samp.genes.list$DE.genes
  
  
  SampCol <- function(n.col, k.ind, sort.method, treatment){
    #sample columns (replicates) to be used
    if (sort.method == "paired" && switch.trt == FALSE){
      odds <- seq(1, n.col, by = 2)
      samp <- sample(odds, 2*k.ind, replace = FALSE)
      samp <- c(samp, samp[-(1:k.ind)] + 1)
    }
    else if (sort.method == "paired" && switch.trt == TRUE){
      odds <- seq(1, n.col, by = 2)
      samp <- sample(odds, 2*k.ind, replace = FALSE)
      samp <- c(samp[1:k.ind], samp + 1)
    }
    else if (sort.method == "unpaired" && switch.trt == FALSE){  
      trt1 <- 1:table(treatment)[1]
      trt2 <- (table(treatment)[1] + 1):length(treatment)
      samp <- c(sample(trt1, 2*k.ind), sample(trt2, k.ind))
    }
    else if (sort.method == "unpaired" && switch.trt == TRUE){  
      trt1 <- 1:table(treatment)[1]
      trt2 <- (table(treatment)[1] + 1):length(treatment)
      samp <- c(sample(trt1, k.ind), sample(trt2, 2*k.ind))
    }
    samp <- as.numeric(samp)
    return(samp)
  }
  
  #swap in appropriate data to create matrix with known DE and EE genes
  #multiply by appropriate norm.factors and round

  
  samp.col <- NULL
  
  if (samp.independent == FALSE){
    samp <- SampCol(n.col, k.ind, sort.method, treatment)
    samp.col <- sorting[samp]
    norm.factors.col <- norm.factors[samp]
    
    #perform swaps
    data <- counts[genes.subset, samp, drop = FALSE]
    
    if(switch.trt == FALSE){
      data.table <- data[, 1:(2*k.ind), drop = FALSE]
      data.table[DE.genes, (k.ind+1):(2*k.ind)] <- 
        round(t(t(data[DE.genes, (2*k.ind + 1):(3*k.ind)])/norm.factors.col[(2*k.ind + 1):(3*k.ind)]*norm.factors.col[(k.ind+1):(2*k.ind)]))
    } 
    else if (switch.trt == TRUE){
    data.table <- data[, (k.ind + 1):(3*k.ind), drop = FALSE]
    data.table[DE.genes, 1:k.ind] <- 
      round(t(t(data[DE.genes, 1:k.ind])/norm.factors.col[1:k.ind]*norm.factors.col[(k.ind + 1):(2*k.ind)]))
    }
  }
  
  #swap in appropriate data to create matrix with known DE and EE genes
  #multiply by appropriate norm.factors and round
  #under independent method, sample columns separately and independently
  #for each gene

  if (samp.independent == TRUE){
    samp <- t(replicate(n.genes, SampCol(n.col, k.ind, sort.method, treatment)))
    norm.factors.col <- apply(samp, 2, function(x) norm.factors[x])
    
    #perform swaps
      data.temp <- counts[genes.subset, , drop = FALSE]
      data <- matrix(mapply(function(x, y) data.temp[x, y, drop = FALSE], x = 1:n.genes, 
                             y = samp, SIMPLIFY = TRUE), ncol = 3*k.ind)
      normalized <- data/norm.factors.col
      if(switch.trt == FALSE){
        data.table <- normalized[, 1:(2*k.ind), drop = FALSE]
        if (n.diff > 0){
          data.table[DE.genes, (k.ind + 1):(2*k.ind)] <- normalized[DE.genes, (2*k.ind + 1):(3*k.ind), drop = FALSE]
        }
      }
      if(switch.trt == TRUE){
        data.table <- normalized[, (k.ind + 1):(3*k.ind), drop = FALSE]
        if (n.diff > 0){
          data.table[DE.genes, 1:k.ind] <- normalized[DE.genes, 1:k.ind, drop = FALSE]
        }
      }
    
    #unnormalize dataset
    data.table <- round(t(t(data.table)*norm.factors[sample(1:n.col, 2*k.ind, replace = FALSE)]))
  }
  trt <- c(rep(1, k.ind), rep(2, k.ind))
  
  ### remove column attributes
  colnames(data.table) <- NULL
  
  return(list(counts = data.table, treatment = c(rep(0, k.ind), rep(1, k.ind)), genes.subset = genes.subset, 
              DE.genes = genes.diff, DE.ind = genes.subset %in% genes.diff, col = samp.col))
}