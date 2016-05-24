ijmatrix.shares <-
function (rawmatrix, submarkets, suppliers, observations) {
mcirawmatrix <- rawmatrix   
submarkets_single <- levels(as.factor(mcirawmatrix[[submarkets]]))
suppliers_single <- levels(as.factor(mcirawmatrix[[suppliers]]))
submarkets_colname <- names(rawmatrix[submarkets])
suppliers_colname <- names(rawmatrix[suppliers])
submarkets_count <- nlevels(as.factor(mcirawmatrix[[submarkets]]))   
suppliers_count <- nlevels(as.factor(mcirawmatrix[[suppliers]]))   

p_ij_obs <- vector()
freq_ij_rel <- vector()
freq_i_abs <- vector()
freq_ij_rel_j <- vector()
submarket_i_total <- 0

for(i in 1:submarkets_count){   
submarket_i <- subset (mcirawmatrix, mcirawmatrix[[submarkets_colname]] == submarkets_single[i])  
submarket_i_total[i] <- sum (submarket_i[[observations]])
for(j in 1:suppliers_count) {   
freq_ij_rel_j[j] <- submarket_i[[observations]][j]/submarket_i_total[i]
freq_i_abs <- rbind(freq_i_abs, as.numeric(list(submarket_i_total[i])))
}   
}   

mcirawmatrix$freq_i_total <- freq_i_abs 
mcirawmatrix$p_ij_obs <- mcirawmatrix[[observations]]/mcirawmatrix$freq_i_total
return(mcirawmatrix)
}
