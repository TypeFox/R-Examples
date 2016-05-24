test.rms <- function(observed,expected){
#Computes root-mean-square test statistic between observed and model 
#genotypic counts
rms <- 0
obs <- as.matrix(observed,rownames.force=T)
exp <- as.matrix(expected,rownames.force=T)
r <- nrow(obs)
count <- 1
while(count <= r){
for(s in count:r){
rms <- rms + (observed[s,count] - expected[s,count])^2
}
count <- count + 1
}
return(rms)
}
