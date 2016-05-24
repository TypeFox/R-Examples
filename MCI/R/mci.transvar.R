mci.transvar <-
function (mcidataset, submarkets, suppliers, mcivariable, output_ij = FALSE, output_var = "numeric") {   
sort_i_j <- order(mcidataset[[submarkets]], mcidataset[[suppliers]])
mciworkfile <- mcidataset[sort_i_j,]   

if (var.check(mciworkfile[[mcivariable]]) == "valid_d") {
return(mciworkfile[mcivariable]) 
}

if (var.check(mciworkfile[[mcivariable]]) == "valid_n") {
logvarnewname <- paste(names(mciworkfile[mcivariable]), "_t", sep="")   
submarkets_single <- levels(as.factor(mciworkfile[[submarkets]]))   
suppliers_single <- levels(as.factor(mciworkfile[[suppliers]]))
submarkets_count <- nlevels(as.factor(mciworkfile[[submarkets]]))
suppliers_count <- nlevels(as.factor(mciworkfile[[suppliers]]))   

submarket_i_geom <- 0   
submarket_i_rel <- 0   
submarket_i_rel_log <- 0   
mcivariablelog <- vector()  

for(i in 1:submarkets_count){   
submarket_i <- subset (mciworkfile, mciworkfile[[submarkets]] == submarkets_single[i])  
submarket_i_geom[i] <- geom (submarket_i[[mcivariable]])  
for(j in 1:suppliers_count) {   
submarket_i_rel[j] <- submarket_i[[mcivariable]][j]/submarket_i_geom[i]  
submarket_i_rel_log[j] <- log10(submarket_i_rel[j])   
mcivariablelog <- rbind(mcivariablelog, list(submarket_i_rel_log[j]))  
}   
}   

mcilinvar <- as.data.frame(mcivariablelog)   
names(mcilinvar) <- logvarnewname   

if (output_ij == TRUE) {
mcilinoutput <- cbind(mciworkfile[submarkets], mciworkfile[suppliers], mcilinvar)
if (output_var == "numeric") { mcilinoutput[3] <- as.numeric(unlist(mcilinoutput[3])) } 
return(mcilinoutput)
}
else {
if (output_var == "numeric") { return(as.numeric(unlist(mcilinvar))) }
if (output_var == "list") { return(mcilinvar) }  
}
}   
else {
return(0)
}   
}
