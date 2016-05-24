percentilsArchetypoid <- function(column,indiv,data,digits){
 aux1 <- data[,colnames(data)[column]]  
 
 aux2 <- as.matrix(aux1)  
 
 Fn <- ecdf(aux2)
 
 aux3 <- Fn(data[indiv,colnames(data)[column]])
 
 round(aux3 * 100, digits = digits)
}
