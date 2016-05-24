HW.plain <- function(model_dist,rms,chisq,gsq,T,n){
#Computes the plain p-value for the observed genotype counts
r <- ncol(model_dist)
rarities_rms <- 0
rarities_chisq <- 0
rarities_gsq <- 0
pval_rms <- 0
pval_chisq <- 0
pval_gsq <- 0

#Run T monte carlo simulations
for(h in 1:T){
#Draw n genotypes i.i.d. from model_dist by first converting it to linear array
count <- 1
count_lmd <- 1
lin_model_dist <- array(0,c(1,r*(r+1)/2))
while(count <= r){
for(s in count:r){
lin_model_dist[1,count_lmd] <- model_dist[s,count]
count_lmd <- count_lmd + 1
}
count <- count + 1
}
lin_model_dist <- lin_model_dist/n
g <- ncol(lin_model_dist)
cdvec <- cumsum(lin_model_dist)
lin_samp_model_dist <- array(0, c(1,g))
for (i in 1:n){
t <- runif(1)
ind <- which(cdvec > t)
first <- ind[1]
lin_samp_model_dist[1,first] <- lin_samp_model_dist[1,first] + 1
}
samp_model_dist <- mat.or.vec(r,r)
count <- 1
count_smd <- 1
while(count <= r){
for(s in count:r){
samp_model_dist[s,count] <- lin_samp_model_dist[1,count_smd]
count_smd <- count_smd + 1
}
count <- count + 1
}

#Create model distribution based on samp_model_dist
tilde_model_dist <- create.model(samp_model_dist,n)

#Test stats with samp_model_dist and tilde_model_dist
sim_rms <- test.rms(samp_model_dist,tilde_model_dist)
sim_chisq <- test.chisq(samp_model_dist,tilde_model_dist)
sim_gsq <- test.gsq(samp_model_dist,tilde_model_dist)

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
