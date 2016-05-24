# Relational class analysis
RCA <- function(matrix, max=NULL, min=NULL, num=1000, alpha=0.05) {
  UseMethod("RCA")
}

RCA.default <- function(matrix, max=NULL, min=NULL, num=1000, alpha=0.05) {
  matrix <- as.matrix(matrix)
  K <- dim(matrix)[2]
  R <- relationality(matrix, max, min)
  R <- removeInsignificantEdges(R, num, alpha)
  graph <- createIGraph(R, abs=TRUE)
  eigcom <- leading.eigenvector.community(graph)
  modules <- partitionData(matrix, eigcom$membership) 
  val <- list(membership=eigcom$membership, modules=modules, R=R)
  class(val) <- "RCA"
  return(val)
}

# Check if object is of class "RCA"
is.RCA <- function(x) inherits(x, "RCA")

# Print method
print.RCA <- function(x, ...) {
  cat("RCA found", length(unique(x$membership)), "relational classes. Sizes:", table(x$membership), "\n")
}

# Plot method
plot.RCA <- function(x, module=NULL, colorblind=FALSE, heatmap=TRUE, heat_labels=FALSE, drop_neg_ties=TRUE, 
                     layout=layout.kamada.kawai, edge_color="gray", vertex_color="white", vertex_frame_color="black",
                     vertex_size=20, vertex_label_color="black", vertex_label_cex=0.8, margin=0, ...) {
  if (is.null(module)) {
    warn <- warning("Please specify a module to plot.")
  } else {
    palette <- colorRampPalette(c("red", "white", "mediumseagreen"))(n=299)
    if (colorblind==TRUE) palette <- colorRampPalette(c("red", "white", "royalblue"))(n=299)
    col_breaks = c(seq(-1, 1, length=300)) 
    mod <- x$modules[module][[1]]
    mod <- data.frame(mod$matrix)
    mod <- round(cor(mod), digits=2)
    if (heatmap) {
      if (heat_labels) {
        heatmap.2(mod, cellnote=mod, notecol="black", density.info="none", trace="none", col=palette, 
                  breaks=col_breaks, symm=TRUE, revC=TRUE)
      } else {
        heatmap.2(mod, notecol="black", density.info="none", trace="none", col=palette, 
                  breaks=col_breaks, symm=TRUE, revC=TRUE)
      } 
    } else {
      ig <- createIGraph(cor=mod, abs=TRUE)
      if (drop_neg_ties) ig <- delete_edges(ig, E(ig)[E(ig)$weight < 0])	
      width <- (E(ig)$weight) * 10
      plot(ig, layout=layout, edge.lty=1, edge.width=width, edge.color=edge_color, vertex.color=vertex_color, 
           vertex.frame.color=vertex_frame_color, vertex.size=vertex_size, vertex.label.color=vertex_label_color, 
           vertex.label.cex=vertex_label_cex, margin=margin-0.3)
    } 
  }
}

# Relationality function
# Takes 3 inputs, only the first of which is necessary for the code to run. Inputs are: 
# (1) a matrix of data where the rows correspond to observations and the columns correspond to variables
# (2) a max value for the normalization - if not specified, then the max value in the matrix is used
# (3) a min value for the normalization - if not specified, then the min value in the matrix is used
relationality <- function(matrix, max=NULL, min=NULL) {
  
  # Define a set of constants to use in the code
  nrows = dim(matrix)[1]
  K = dim(matrix)[2]
  constant = 2/(K*(K-1))
  
  # Check if a min & max value have been specified
  # If not, assign them the max & min values of the input matrix
  # If only a single min & max is specified, use this for all columns
  if (is.null(max)) max <- apply(matrix, 2, max)
  if (is.null(min)) min <- apply(matrix, 2, min)
  if (length(max)==1) max <- rep(max, times=K)
  if (length(min)==1) min <- rep(min, times=K)

  # Standardize each variable in the dataset using the maximum and minimum values 
  maxMinNorm <- function(X) {
    if (all(max > min)) {
      for (i in 1:length(max)) {
        if (max[i] != min[i]) X[,i] <- (X[,i] - min[i]) / (max[i] - min[i])
      }
      return(X)
    } else {
      return(X - X)
    }
  }
  matrix <- maxMinNorm(matrix) 
  
  # Turn the matrix into a list, where each index, [[i]], corresponds to the ith row of the matrix 
  # (This is required to use the lapply function later)
  # X^k_i, k = column/var, i = row/obs
  X <- lapply(seq_len(nrow(matrix)), function(i) matrix[i,])
  
  # Distance between the values of variables k and l for observation i 
  # \DeltaX^(kl)_i, i = list index [[]], kl pair = column 
  deltaXFunction <- function(i) {
    x <- X[[i]]
    combn <- combn(x, 2) # all combinations of the elements of x taken 2 at a time
    return(combn[1,] - combn[2,])
  }
  deltaX <- lapply(1:length(X), function(i) deltaXFunction(i))
  
  # Change deltaX from a list of matrix rows to a list of matrix columns, take it's absolute value
  # \DeltaX^(kl)_i, kl pair = list index [[]], i = column 
  deltaX <- matrix(unlist(deltaX), nrow=length(deltaX), byrow=TRUE)
  deltaX <- lapply(seq_len(ncol(deltaX)), function(i) deltaX[,i])
  absDeltaX <- lapply(deltaX, abs)
  
  # Schematic similarity for the variable pair kl between observations i & j
  # \delta^(kl)_ij, kl pair = list #, ij pair = column
  deltaFunction <- function(i) {
    x <- absDeltaX[[i]]
    combn <- combn(x, 2) # all combinations of the elements of x taken 2 at a time
    absDiff <- abs(combn[1,] - combn[2,])
    return(1 - absDiff)
  }
  delta <- lapply(1:length(absDeltaX), function(i) deltaFunction(i))
  
  # Binary coefficient that changes the sign of the schematic similarity if both distances are in the opposite direction
  # \lambda^(kl)_ij, kl pair = list #, ij pair = column
  lambdaFunction <- function(i) {
    x <- deltaX[[i]]
    combn <- combn(x, 2) # all combinations of the elements of x taken 2 at a time
    diff <- combn[1,] * combn[2,]
    return(ifelse(diff >= 0, yes=1, no=-1))
  }
  lambda <- lapply(1:length(deltaX), function(i) lambdaFunction(i))
  
  # Lambda * delta (binary coefficient * schematic similarity)
  # \lambda^(kl)_ij * \delta^(kl)_ij, ij pair = list #, kl pair = column
  lambdaDelta <- lapply(1:length(delta), function(i) delta[[i]] * lambda[[i]])
  lambdaDelta <- matrix(unlist(lambdaDelta), nrow=length(lambdaDelta), byrow=TRUE)
  lambdaDelta <- lapply(seq_len(ncol(lambdaDelta)), function(i) lambdaDelta[,i])
  
  # Relationality between observations i & j
  # R_(ij) = 2/(K(K-1)) * sum_(k=1)^(K-1)(sum_(l=k+1)^K(lambda_(ij)^(kl)*delta_(ij)^(kl)))
  r <- constant * unlist(lapply(lambdaDelta, sum))
  R1 <- matrix(0, ncol=nrows, nrow=nrows)
  R2 <- matrix(0, ncol=nrows, nrow=nrows)
  R1[lower.tri(R1, diag=FALSE)] <- r
  R2[lower.tri(R2, diag=FALSE)] <- r
  R2 <- t(R2)
  R <- matrix(ncol=nrows, nrow=nrows)
  R <- ifelse(R1==0, yes=R2, no=R1)
  return(R)
}

# Bootstrap function for relationality matrix
bootstrapR <- function(R) {
  N <- seq(1:dim(R)[1])
  sample <- sample(x=N, size=length(N), replace=TRUE)
  bootstrap <- R[sample,sample]
  return(bootstrap)
}

# Use bootstrapping to eliminate non-significant graph edges
# 'num' bootstrap samples
# (1) For each bootstrap sample, randomly sample w/ replacement N people and recreate the NxN relationality 
# matrix with the corresponding rows (e.g., c(1,1,2,3,4) x c(1,1,2,3,4))
# (2) Find the average relationality between all pairs & the average standard deviation 
# (3) Compute the z-score for all observed relationalities
# (4) Keep the elements of R that are significant at alpha = 0.05, 1.96 above/below the mean 
removeInsignificantEdges <- function(R, num, alpha) {
  means <- numeric(num)
  sds <- numeric(num)
  for (i in 1:num) {
    bootR <- bootstrapR(R) 
    means[i] <- mean(bootR)
    sds[i] <- sd(bootR)
  }
  bootMean <- mean(means)
  bootSD <- mean(sds)
  R <- R - bootMean
  z <- qnorm(1 - alpha / 2)
  zscores <- R / bootSD
  R[zscores > -z & zscores < z] <- 0
  return(R)
}

# Create IGraph
createIGraph <- function(cor, abs=TRUE) { 
  if (abs) cor <- abs(cor)
  diag(cor) <- 0
  graph <- graph.adjacency(cor, mode="undirected", weighted=TRUE, diag=FALSE)
  return(graph)
}

# Partition data according to IGraph membership vector
partitionData <- function(matrix, membership) {
  ids <- sort(unique(membership))
  modules <- list()
  
  for (i in 1:length(ids)) {
    mod <- list()
    class(mod) <- "RCA.module"
    mod$matrix <- matrix[membership==ids[i],]
    modules[[i]] <- mod
  }
  return(modules)
}







