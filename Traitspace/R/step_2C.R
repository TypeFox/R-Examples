step_2C <-
function(P_level_2_level_1,level_1){                                     
P_level_1_level_2_level_3 <- matrix(0,nrow(P_level_2_level_1), nlevels(level_1))             
for (i in 1:nrow(P_level_2_level_1)) {
  P_level_1_level_2_level_3[i,] = exp(log(P_level_2_level_1[i,]) - log(apply( P_level_2_level_1, 1, sum)[i]))
}
result <- list(P_level_1_level_2_level_3 = P_level_1_level_2_level_3)
return(result)
}
