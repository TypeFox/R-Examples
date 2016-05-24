data(kidney)
counts <- kidney$counts # Matrix of read counts from KIRC dataset
replic <- kidney$replic # Replic vector indicating paired columns
treatment <- kidney$treatment # Treatment vector indicating Non-Tumor or Tumor columns

nf <- apply(counts, 2, quantile, 0.75)

require(fdrtool)

  
### Example 1: Simulate Matrix with 1000 DE genes and 4000 EE genes
data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 5000, n.diff = 1000,
                    norm.factors = nf)

### Example 2: Calculate weights vector beforehand to save run time in
### repeated simulations
sort.list <- SortData(counts = counts, treatment = treatment, replic = replic,
                      sort.method = "paired", norm.factors = nf)
counts <- sort.list$counts
replic <- sort.list$replic
treatment <- sort.list$treatment
nf <- sort.list$norm.factors

probs <- CalcPvalWilcox(counts, treatment, sort.method = "paired", 
                        sorted = TRUE, norm.factors = nf, exact = FALSE)
weights <- 1 - fdrtool(probs, statistic = "pvalue", plot = FALSE, verbose = FALSE)$lfdr 

data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 5000, n.diff = 1000,
                    weights = weights, norm.factors = nf)

### Example 3: Specify which genes you want to use in the simulation

# Randomly sample genes or feed in the exact genes you wish to use
genes.diff <- sample(1:nrow(counts), size = 1000, prob = weights)
genes <- c(sample(1:nrow(counts)[-genes.diff], 4000), genes.diff)

data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, genes.select = genes,
                    genes.diff = genes.diff, weights = weights, norm.factors = nf)

### Example 4: Simulate matrix with DE genes having log base 2 fold change greater than 1

# add one to counts matrix to avoid infinities when taking logs
tumor.mean <- rowMeans(log2((counts[, treatment == "Tumor"] + 1) %*% diag(1/nf[treatment == "Tumor"])))
nontumor.mean <- rowMeans(log2((counts[, treatment == "Non-Tumor"] + 1) %*% diag(1/nf[treatment == "Non-Tumor"])))
lfc <- tumor.mean - nontumor.mean
weights.zero <- abs(lfc) < 1
weights[weights.zero] <- 0

data.sim <- SimData(counts = counts, replic = replic, treatment = treatment, 
                    sort.method = "paired", k.ind = 5, n.genes = 5000, n.diff = 1000,
                    weights = weights, norm.factors = nf)

### Example 5: Simulate three treatment groups:
### 3 Different types of Differential Expression Allowed
### First Group Diff, Second and Third group Equal
### Second Group Diff, First and Third group Equal
### Third Group Diff, First and Second group Equal

k <- 5 # Sample Size in Each treatment group

### Sample DE genes beforehand
N <- nrow(counts)
genes.de <- sample(1:N, size = 1000, prob = weights) # Sample all DE genes
DE1 <- genes.de[1:333] # Sample DE genes with first trt diff
DE2 <- genes.de[334:666] # Sample DE genes with sec trt diff
DE3 <- genes.de[667:1000] # Sample DE genes with third trt diff
EE <- sample( (1:N)[-genes.de], size = 4000) #Sample EE genes

genes.tot <- c(EE, genes.de)
genes.de1 <- union(DE2, EE) #Assign DE genes for first sim
genes.de2 <- union(DE2, DE3) #Assign DE genes for second sim

data.sim1 <- SimData(counts = counts, replic = replic, treatment = treatment, 
                     sort.method = "paired", k.ind = k, genes.select = genes.tot,
                     genes.diff = genes.de1, weights = weights, norm.factors = nf)

#remove pairs of columns used in first simulation
cols.rm <- c(data.sim1$col[1:(2*k)], data.sim1$col[1:(2*k)] + 1)
counts.new <- counts[, -cols.rm]
nf.new <- nf[-cols.rm]
replic.new <- replic[-cols.rm]
treatment.new <- treatment[-cols.rm]

### Set switch.trt = TRUE for second sim
data.sim2 <- SimData(counts = counts.new, replic = replic.new, treatment = treatment.new, 
                     sort.method = "paired", k.ind = k, genes.select = genes.tot,
                     genes.diff = genes.de2, weights = weights, norm.factors = nf.new,
                     switch.trt = TRUE)

### Remove first k.ind entries from first sim and combine two count matrices
counts.sim <- cbind(data.sim1$counts[, -(1:k)],  data.sim2$counts)

### treatment group levels for simulated matrix
trt.grp <- rep(NA, 5000)
trt.grp[is.element(data.sim1$genes.subset, DE1)] <- "DE_First_Trt"
trt.grp[is.element(data.sim1$genes.subset, DE2)] <- "DE_Second_Trt"
trt.grp[is.element(data.sim1$genes.subset, DE3)] <- "DE_Third_Trt"
trt.grp[is.element(data.sim1$genes.subset, EE)] <- "EE"