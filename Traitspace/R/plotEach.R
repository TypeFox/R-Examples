plotEach <-
function(result, obs = NULL, byrow = TRUE, compare = TRUE){

if(is.null(obs)& byrow) obs <- result$true.p$p
if(is.null(obs)& !byrow) obs <- result$true.p$p_dist

op <- par(ask=TRUE)
if(byrow){
for (i in 1:nrow(result$predicted.p$P_level_1_level_3)){
plot(result$predicted.p$P_level_1_level_3[i,], pch = 1, axes = FALSE, frame.plot = TRUE, ylim = c(0,1), ylab = "Relative abundances",xlab = "Species", main = rownames(result$step_2D$P_level_1_level_3)[i])
title(main = rownames(result$predicted.p$P_level_1_level_3)[i])
axis(1,at=1:ncol(result$predicted.p$P_level_1_level_3), labels = colnames(result$predicted.p$P_level_1_level_3))
axis(2,at= seq(0,1,0.1), labels = seq(0,1,0.1))
if(compare == TRUE){
points(obs[i,], pch = 3)
legend(1, 1, c("prediction", "observation"), pch = c(1,3))
}
}
}else{
for (i in 1:ncol(result$predicted.p$P_level_1_level_3_dist)){
plot(result$predicted.p$P_level_1_level_3_dist[,i], pch = 1, axes = FALSE, frame.plot = TRUE, ylim = c(0,1), ylab = "Species distributions",xlab = "Sites", main = colnames(result$step_2D$P_level_1_level_3)[i])
title(main = colnames(result$predicted.p$P_level_1_level_3)[i])
axis(1,at=1:nrow(result$predicted.p$P_level_1_level_3_dist), labels = rownames(result$predicted.p$P_level_1_level_3_dist))
axis(2,at= seq(0,1,0.1), labels = seq(0,1,0.1))
if(compare == TRUE){
points(obs[,i], pch = 3)
legend(1, 1, c("prediction", "observation"), pch = c(1,3))
}
}
}
par(op)
}
