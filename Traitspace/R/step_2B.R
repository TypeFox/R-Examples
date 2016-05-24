step_2B <-
function(par_2_1, level_2_sample,level_1){
P_level_2_level_1 <- matrix(0,nrow(level_2_sample),nlevels(level_1))
level_2_sample <- matrix(data = level_2_sample, nrow = nrow(level_2_sample), ncol = ncol(level_2_sample))
for (i in 1:nlevels(level_1)) {                
  P_level_2_level_1[,i] <- dens(par_2_1[[i]]$variance$modelName, data = level_2_sample, parameters=par_2_1[[i]])         
}
result <- list(P_level_2_level_1 = P_level_2_level_1)
return(result)
}
