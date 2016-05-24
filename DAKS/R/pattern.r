#####################
# pattern frequency #
#####################

############################################
#                                          #
# This function computes the absolute      # 
# frequencies of the response patterns,    # 
# and optionally, the absolute frequencies #
# of a collection of specified knowledge   #
# states in a dataset.                     #
#                                          #
############################################

pattern<-function(dataset, n = 5, P = NULL){

pattern<-sort(table(apply(dataset,1, function(x) paste(x, collapse = ""))), decreasing = TRUE)

if(n < 1) stop("Number of patterns must be greater than zero.\n")
if(n > length(pattern)) n = length(pattern)

if(is.null(P)){
out<-list(response.patterns = pattern[1:n], states = P, n = n)
}else{
states<-cbind(P, 0)
states[,ncol(states)] <- sapply(apply(P, 1, function(x) pattern[names(pattern) == paste(x, collapse = "")]), function(y) max(0, y))
if(is.null(names(dataset)) == FALSE){
colnames(states)<-c(names(dataset), "size")
rownames(states)<-c(paste("State", 1:(nrow(states))))
}else{
colnames(states)<-c(paste("Item", 1:(ncol(dataset))), "size")
rownames(states)<-c(paste("State", 1:(nrow(states))))
}
out<-list(response.patterns = pattern[1:n],states = states, n = n)
}
class(out)<-"pat"
out
}
