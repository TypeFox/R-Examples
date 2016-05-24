# A. Creating a matrix
aa <- 1:1000
bb <- matrix(data = aa, nrow = 1); class(bb); mode(bb)
cc <- data.frame(bb)
library(erer); lss()  # size

# diag(): four usuages
ma <- matrix(data = 1:20, nrow = 4, ncol = 5, byrow = FALSE); ma
diag(x = ma)                     # 1. extract the diagonal values
diag(diag(x = ma))               #    extract the diagnonal matrix
diag(nrow = 4)                   # 2. create an identity matrix
diag(x = 4)                      # 3. create an identify matrix
diag(x = c(3, 9, 10))            # 4. a matrix with the given diagonal
diag(x = c(3, 9, 10), nrow = 4)  #    a matrix with more rows
mb <- ma 
diag(x = mb) <- c(NA, 0, 1, NA); mb

# coercion
as.matrix(1:10) 
as.vector(ma); c(ma); identical(as.vector(ma), c(ma))
v <- c(40, 3, 2); names(v) <- c("H", "I", "K"); v
as.vector(v); c(v); identical(as.vector(v), c(v))
rbind(ma, mb); cbind(ma, mb)  # combine
            
# B. Matrix attributes 
mc <- ma
dim(mc); nrow(mc); ncol(mc)
dimnames(mc) <- list(myrow = letters[1:nrow(mc)], mycol = LETTERS[1:5])
dimnames(mc); rownames(mc); colnames(mc)
attributes(mc)

# C1. Matrix indexing: two indexes
mc[1:3, c("A", "E")]
mc[c(TRUE, TRUE, TRUE, FALSE), c(TRUE, FALSE, FALSE, FALSE, TRUE)]
mc[-c(1, 2), ]
ia <- mc[1, ]; ia; class(ia)
ib <- mc[1, , drop = FALSE]; ib; class(ib)

# C2. Matrix indexing: one single index
mc[1]; mc[19:20]
identical(as.vector(mc[3:4, 5]), ma[19:20])
mc[mc > 18]; mc[as.vector(mc > 18)]; mc[c(mc > 18)]
upp <- upper.tri(x = mc, diag = TRUE); upp
low <- lower.tri(x = mc, diag = FALSE); low
ta <- mc[upp]; ta
tb <- mc; tb[low] <- 0; tb
mc[row(mc) == (col(mc) - 1)]