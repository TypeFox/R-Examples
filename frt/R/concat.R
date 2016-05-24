concat <-
function(x1, x2){
lst <- NA
for (i in 1:dim(x1)[1]){
for(j in 1:dim(x2)[1]){
lst <- rbind(lst,c(x1[i,],x2[j,]))
}
}
return(lst[-1,])
}

