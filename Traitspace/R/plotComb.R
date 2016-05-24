plotComb <-
function(result, obs = NULL, byrow = TRUE, compare = TRUE){
if(is.null(obs)) obs <- result$true.p$p_dist
if(!all(dim(result$predicted.p$P_level_1_level_3)==dim(obs))) stop("Different dimensions!")

if(length(as.matrix(unique(result$check$level_3))) == nlevels(result$check$site.name)){
env <- unique(result$check$level_3)
}else{
env <- as.matrix(as.integer(unique(result$check$site)$site))
}
species <- as.integer(unique(result$check$level_1))
if(byrow){
for (i in 1:nrow(result$predicted.p$P_level_1_level_3)){
if(i %% 6 == 1 ){
#x11()
par(mfrow=c(3,2))
}
plot(species, result$predicted.p$P_level_1_level_3[i,], pch = 1, axes = FALSE, frame.plot = TRUE, col = "red", ylim=c(0,1) ,
ylab = "Relative abundances", main = rownames(result$predicted.p$P_level_1_level_3)[i])
axis(1,at=1:ncol(result$predicted.p$P_level_1_level_3), labels = colnames(result$predicted.p$P_level_1_level_3))
axis(2,at= seq(0,1,0.1), labels = seq(0,1,0.1))
if(compare){points(species, obs[i,], pch = 3)}
}

}else{
for(j in 1:ncol(env)){
for (i in 1:ncol(result$predicted.p$P_level_1_level_3_dist)){
if(i %% 6 == 1 ){
#x11()
par(mfrow=c(3,2))
}
plot(env[,j], result$predicted.p$P_level_1_level_3_dist[,i], pch = 1, frame.plot = TRUE, col = "red", 
ylim = c(0,1),xlim = range(env[,j]), ylab = "Species distributions",xlab = names(env)[j], main = colnames(result$predicted.p$P_level_1_level_3_dist)[i])
if(compare){points(env[,j], obs[,i], pch = 3)}
}
}
}
par(mfrow=c(1,1))
}
