pval <-
function(result, obs = NULL, byrow = TRUE, permutations = 999){
if(is.null(obs)& byrow) obs <- result$true.p$p
if(is.null(obs)& !byrow) obs <- result$true.p$p_dist
if(byrow) predicted <- result$predicted.p$P_level_1_level_3
if(!byrow) predicted <- result$predicted.p$P_level_1_level_3_dist

# reference distances
dist <- distTraitspace(result, byrow =  byrow)
mean.dist <- lapply(dist, mean)

# Permutation
set.seed(12345)
#library(permute)
if (byrow){
nr <- nrow(obs)
perm <- vector("list", permutations)
perm_site  <- vector("list", permutations)
for (i in 1:permutations) {
take <- shuffle(nr)
result$predicted.p$P_level_1_level_3 <- predicted[take,]
perm_site[i] <- list(i = distTraitspace(result, byrow =  byrow))
perm[i] <- list(i = (lapply(distTraitspace(result,  byrow =  byrow),mean)))
}

}else {
nc <- ncol(obs)
perm <- vector("list", permutations) # ?numPerms(N)
perm_site  <- vector("list", permutations)
for (i in 1:permutations) {
take <- shuffle(nc)
result$predicted.p$P_level_1_level_3_dist <- predicted[,take]
perm_site[i] <- list(i = distTraitspace(result, byrow =  byrow))
perm[i] <- list(i = (lapply(distTraitspace(result, byrow =  byrow),mean)))
}
}
warnings()

# P-value
Pval <- vector("list", length(mean.dist))
#x11()
par(mfrow=c(3,2))
for (j in 1:length(mean.dist)){
Pval[[j]] <- (sum(sapply(perm, "[[", j) <= mean.dist[[j]])+1)/(permutations + 1)
min.p <- min(sapply(perm, "[[", j), mean.dist[[j]])
max.p <- min(max(sapply(perm, "[[", j), mean.dist[[j]]),10)
hist(sapply(perm, "[[", j), main = names(perm[[1]])[j], xlim = c(min.p, max.p))
abline(v = mean.dist[[j]], col = "red")
}
names(Pval) <- names(mean.dist)
# P-value for each site
# sapply(perm_site, "[[", 1)[1,]
Pval_site <- vector("list", length(mean.dist))

if(byrow)site.num <- nrow(obs)
if(!byrow)site.num <- ncol(obs)

for (j in 1:length(mean.dist)){
for (i in 1:site.num){
Pval_site[[j]] <- c(Pval_site[[j]],(sum(sapply(perm_site, "[[", j)[i,] <= dist[[j]][i])+1)/(permutations + 1))
}}
names(Pval_site) <- names(mean.dist)
result= list(Pval = Pval, mean.dist = mean.dist, Pval_site = Pval_site, dist_site = dist) # perm_site = perm_site, perm = perm,
return(result)
}
