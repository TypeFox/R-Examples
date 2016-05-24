ctoc <-
function(y,x,data,min.obs){
for (i in seq(1,length(x),1)) {
if (class(data[,x[i]]) %in% c("factor", "character")){
cnt <- NULL
temp1 <- aggregate(data[,y], list(data[,x[i]]), mean)
names(temp1)[2] <- "avg_mean"
temp2 <- aggregate(data[,y], list(data[,x[i]]), length)
names(temp2)[2] <- "cnt"
temp <- merge(temp1, temp2, by = "Group.1")
temp1 <- subset(temp, cnt >= min.obs)
temp2 <- subset(temp, cnt < min.obs)
if(nrow(temp2) > 0) 
{
temp2$avg_mean <- sum(temp2$avg_mean*temp2$cnt)/sum(temp2$cnt)
temp2$cnt <- sum(temp2$cnt)
}
temp <- rbind(temp1, temp2)
data <- merge(data, temp, by.x = x[i], by.y = "Group.1", all.x = T)
data[,paste(x[i],"cont", sep="_")] <- ((data$avg_mean * data$cnt) - (data[,y]))/(data$cnt - 1)
data <- data[, !(colnames(data) %in% c("avg_mean","cnt"))]
}
else{
warning(paste(x[i], " is not a factor/character variable", sep = ""))
}
}
return(data)
}
