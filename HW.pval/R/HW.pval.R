HW.pval <- function(genotype_count,num_simulations=10000,type="both"){
#Main function of the package.  Returns plain and/or fully conditional p-value for 
#observed genotype counts
#ensure genotype_count is in matrix form
obs_dist <- genotype_count
T <- num_simulations
if(!is.matrix(obs_dist)){
obs_dist <- as.matrix(obs_dist, rownames.force=T)
}

r <- nrow(obs_dist) #number of alleles
#find number of observed genotypes
n <- 0
count <- 1
while(count<=r){
for(s in count:r){
n <- n + obs_dist[s,count]
}
count<-count+1
}

#create model distribution based on obs_dist
model_dist <- create.model(obs_dist,n) 

rms <- test.rms(obs_dist,model_dist)
chisq <- test.chisq(obs_dist,model_dist)
gsq <- test.gsq(obs_dist,model_dist)

#Plain P-Values
if(type=="plain"){
plain <- HW.plain(model_dist,rms,chisq,gsq,T,n)
cat("Plain P-value: \n \t root-mean-square=",plain[1], "\n \t chi-square=",plain[2], "\n \t likelihood-ratio=",plain[3])
}

#Fully Conditional P-Values
if(type=="cond"){
cond <- HW.cond(obs_dist,model_dist,rms,chisq,gsq,T,n)
cat("Fully conditional P-value: \n \t root-mean-square=",cond[1], "\n \t chi-square=",cond[2], "\n \t likelihood-ratio=",cond[3])
}

#Both, Plain and Fully Conditional P-Values
if(type=="both"){
plain <- HW.plain(model_dist,rms,chisq,gsq,T,n)
cond <- HW.cond(obs_dist,model_dist,rms,chisq,gsq,T,n)
cat("Plain P-value: \n \t root-mean-square=",plain[1], "\n \t chi-square=",plain[2], "\n \t likelihood-ratio=",plain[3], "\n 
Fully conditional P-value: \n \t root-mean-square=",cond[1], "\n \t chi-square=",cond[2], "\n \t likelihood-ratio=",cond[3])
}
}
