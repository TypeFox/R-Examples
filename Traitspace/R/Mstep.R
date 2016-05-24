Mstep <-
function(categ, data){
# Using Mclust to compute P(data/categ)
par <- list()
summary.pdf <- list()

for (i in 1:length(unique(categ))) { 
temp = log(data[which(categ == unique(categ)[i]),])
#mclust.options ONLY because there is a bug in mclust code currently (June2015)
#mclust.options(emModelNames = c("EII", "VII", "EEI", "VEI", "EVI", "VVI", "EEE")) #, "EEV", "VEV", "VVV"))
pdf <- Mclust(temp, modelNames=c("EII", "VII", "EEI", "EVI", "VEI", "VVI", "EEE", "EEV", "VEV", "VVV"))
par[i] <- pdf[12]
summary.pdf[i] <- summary(pdf) $ bic 
}
result <- list(par = par, summary.pdf = summary.pdf)
return(result)
}
