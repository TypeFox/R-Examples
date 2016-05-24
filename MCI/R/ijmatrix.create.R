ijmatrix.create <-
function (rawdataset, submarkets, suppliers) {
mcirawdata <- rawdataset   
submarkets_single <- levels(as.factor(mcirawdata[[submarkets]]))
suppliers_single <- levels(as.factor(mcirawdata[[suppliers]]))
matrix_ij <- merge (submarkets_single, suppliers_single)
submarkets_colname <- names(rawdataset[submarkets])
suppliers_colname <- names(rawdataset[suppliers])
names(matrix_ij) <- c(submarkets_colname, suppliers_colname)
matrix_ij$interaction <- paste(matrix_ij[[submarkets_colname]], "-", matrix_ij[[suppliers_colname]], sep="")   
mcirawdata$interaction <- paste(mcirawdata[[submarkets]], "-", mcirawdata[[suppliers]], sep="")
interactions <- mcirawdata$interaction
interactions_count <- as.data.frame(table(interactions))
names(interactions_count) <- c("interaction", "freq_ij_abs")

mciworkfile <- merge (matrix_ij, interactions_count, by="interaction", all=TRUE)
mciworkfile$freq_ij_abs[is.na(mciworkfile$freq_ij_abs)] <- 0

submarkets_count <- nlevels(as.factor(mcirawdata[[submarkets]]))   
suppliers_count <- nlevels(as.factor(mcirawdata[[suppliers]]))   
p_ij_obs <- vector()
freq_ij_rel <- vector()
freq_i_abs <- vector()
freq_ij_rel_j <- vector()
submarket_i_total <- 0

for(i in 1:submarkets_count){   
submarket_i <- subset (mciworkfile, mciworkfile[[submarkets_colname]] == submarkets_single[i])  
submarket_i_total[i] <- sum (submarket_i$freq_ij_abs)
for(j in 1:suppliers_count) {   
freq_ij_rel_j[j] <- submarket_i$freq_ij_abs[j]/submarket_i_total[i]
freq_i_abs <- rbind(freq_i_abs, as.numeric(list(submarket_i_total[i])))
}   
}   
mciworkfile$freq_i_total <- freq_i_abs 
mciworkfile$p_ij_obs <- mciworkfile$freq_ij_abs/mciworkfile$freq_i_total
return(mciworkfile)
}
