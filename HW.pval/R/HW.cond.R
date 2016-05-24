HW.cond <- function(obs_dist,model_dist,rms,chisq,gsq,T,n){
#Computes the fully-conditional p-value for the observed genotype counts
r <- ncol(obs_dist)
rarities_rms <- 0
rarities_chisq <- 0
rarities_gsq <- 0
pval_rms <- 0
pval_chisq <- 0
pval_gsq <- 0

#Find individual allele counts
allele_count <- array(0,c(1,r))
k <- 1
for(j in k:r){
count <- 1
while(count <= j){
allele_count[1,j] <- allele_count[1,j] + obs_dist[j,count] 
count <- count + 1
}
count <- j
while(count <= r){
allele_count[1,j] <- allele_count[1,j] + obs_dist[count,j]
count <- count + 1
}
}
#Create vec representing no. of each allele present
allele_vec <- array(0,c(1,2*n))
count <- 1
g <- 1
while(count <= r){
for(h in 1:allele_count[1,count]){
if(allele_vec[1,g] == 0){
allele_vec[1,g] <- count
}
g <- g+1
}
count <- count + 1
}

#Run T monte carlo simulations
for(h in 1:T){
#Find model distribution
samp_model_dist <- mat.or.vec(r,r)
samp_allele_vec <- sample(allele_vec) #Permute allele_vec
j <- 1
while(j <= (2*n-1)){
k <- j+1
h <- samp_allele_vec[j]
i <- samp_allele_vec[k]
if(h < i){
temp_h <- h
h <- i
i <- temp_h
}
samp_model_dist[h,i] <- samp_model_dist[h,i] + 1
j <- j + 2
}

#Test stats with samp_model_dist and model_dist
sim_rms <- test.rms(samp_model_dist,model_dist)
sim_chisq <- test.chisq(samp_model_dist,model_dist)
sim_gsq <- test.gsq(samp_model_dist,model_dist)

if(sim_rms > rms){
rarities_rms <- rarities_rms + 1
}
if(sim_chisq > chisq){
rarities_chisq <- rarities_chisq + 1
}
if(sim_gsq > gsq){
rarities_gsq <- rarities_gsq + 1
}
}
pval_rms <- rarities_rms/T
pval_chisq <- rarities_chisq/T
pval_gsq <- rarities_gsq/T

return(c(pval_rms,pval_chisq,pval_gsq))
}
