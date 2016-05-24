`ix.mat2df` <-
function(d) {
	.Deprecated("as.mlds.df", "MLDS")
# d, index matrix for glm approach	
	stim1 <- ifelse(rowSums(d[, -1]) == 0, 0, 1)
	ix.mat <- cbind(stim.1 = stim1, d[, -1])
	dd <- t(apply(ix.mat, 1, function(x) which(x != 0)))
	dd <- as.data.frame(cbind(d[, 1], dd))
	names(dd) <- c("resp", "S1", "S2", "S3", "S4")
	dd	
	}

