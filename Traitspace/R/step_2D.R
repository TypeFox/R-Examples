step_2D <-
function(P_level_2_level_3, P_level_1_level_2_level_3, N,level_1,level_3,pred.site){

P_level_1_level_3_all <- matrix(0,length(P_level_2_level_3),nlevels(level_1))                    
P_level_1_level_3_unnorm <- matrix(0, nrow(unique(level_3)), nlevels(level_1))  ##nlevels(site.name)

#Computing the integr and (with log)
for (i in 1:length(P_level_2_level_3)) {
  P_level_1_level_3_all[i,] <- exp(log(P_level_2_level_3[i])+log(P_level_1_level_2_level_3[i,]))
}

#MC integration and normalisation

for(i in 1:nlevels(level_1)){
  P_level_1_level_3_unnorm[,i] <- as.matrix(tapply(P_level_1_level_3_all[,i], list(pred.site),mean))
}
  P_level_1_level_3 <- t(apply(P_level_1_level_3_unnorm,1,norm<-function(x){return (x/sum(x))}))


colnames(P_level_1_level_3) <- levels(level_1)
rownames(P_level_1_level_3) <- as.character(t(unique(level_3))) #as.character(unique(site.name))

# calculating results in another form
P_level_1_level_3_dist <- apply(P_level_1_level_3_unnorm,2,norm<-function(x){return (x/sum(x))})

colnames(P_level_1_level_3_dist) <- levels(level_1)
rownames(P_level_1_level_3_dist) <- as.character(t(unique(level_3))) #as.character(unique(site.name))

result <- list(P_level_1_level_3 = round(P_level_1_level_3,5), P_level_1_level_3_dist = round(P_level_1_level_3_dist,5))
return(result)
}
