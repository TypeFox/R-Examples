test.gsq <- function(observed,expected){
#Computes log likelihood-ratio test statistic between observed and model 
#genotypic counts
gsq <- 0
obs <- as.matrix(observed,rownames.force=T)
exp <- as.matrix(expected,rownames.force=T)
r <- nrow(obs)
count <- 1
while(count <= r){
for(s in count:r){
if((exp[s,count] == 0) || (obs[s,count] == 0)){
gsq <- gsq
}
if((exp[s,count] != 0) && (obs[s,count] != 0)){
gsq <- gsq + 2*obs[s,count]*log(obs[s,count]/exp[s,count])
}
}
count <- count + 1
}
return(gsq)
}
