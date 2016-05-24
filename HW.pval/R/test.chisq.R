test.chisq <- function(observed,expected){
#Computes Pearsons chi-square test statistic between observed and model 
#genotypic counts
chisq <- 0
obs <- as.matrix(observed,rownames.force=T)
exp <- as.matrix(expected,rownames.force=T)
r <- nrow(obs)
count <- 1
while(count <= r){
for(s in count:r){
if(exp[s,count] == 0){
chisq <- chisq
}
if(exp[s,count] != 0){
chisq <- chisq + (obs[s,count] - exp[s,count])^2/exp[s,count]
}
}
count <- count + 1
}
return(chisq)
}
