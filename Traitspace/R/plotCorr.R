plotCorr <-
function(result, obs = NULL, byrow = TRUE){

if(is.null(obs)& byrow) obs <- result$true.p$p
if(is.null(obs)& !byrow) obs <- result$true.p$p_dist
if(!all(dim(result$predicted.p$P_level_1_level_3)==dim(obs))) stop("Different dimensions!")

op <- par(ask=TRUE)
if(byrow){
for (i in 1:nrow(result$predicted.p$P_level_1_level_3)){
plot(result$predicted.p$P_level_1_level_3[i,], obs[i,], pch = 1, axes = FALSE, frame.plot = TRUE, ylim = c(0,1),xlim = c(0,1), ylab = "Observed relative abundances",xlab = "Predicted relative abundances", main = rownames(result$step_2D$P_level_1_level_3)[i])
title(main = rownames(result$predicted.p$P_level_1_level_3)[i])
text(result$predicted.p$P_level_1_level_3[i,], obs[i,],labels = colnames(result$predicted.p$P_level_1_level_3),pos = 3)
axis(1,at= seq(0,1,0.1), labels = seq(0,1,0.1))
axis(2,at= seq(0,1,0.1), labels = seq(0,1,0.1))
text(0.2,0.9,labels = paste("Correlation:" ,round(cor(result$predicted.p$P_level_1_level_3[i,], obs[i,]),3)))
}
}else{
for (i in 1:ncol(result$predicted.p$P_level_1_level_3_dist)){
plot(result$predicted.p$P_level_1_level_3_dist[,i], obs[,i], pch = 1, axes = FALSE, frame.plot = TRUE, ylim = c(0,1),xlim = c(0,1), ylab = "Observed species distributions",xlab = "Predicted species distributions", main = colnames(result$step_2D$P_level_1_level_3)[i])
title(main = colnames(result$predicted.p$P_level_1_level_3)[i])
text(result$predicted.p$P_level_1_level_3_dist[,i], obs[,i], labels = rownames(result$predicted.p$P_level_1_level_3),pos = 3)
axis(1,at= seq(0,1,0.1), labels = seq(0,1,0.1))
axis(2,at= seq(0,1,0.1), labels = seq(0,1,0.1))
text(0.2,0.9,labels = paste("Correlation:" ,round(cor(result$predicted.p$P_level_1_level_3_dist[,i], obs[,i]),3)))
}
}
par(op)
}
