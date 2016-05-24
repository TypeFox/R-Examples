############
# Plotting #
############

################################################################
#                                                              #
# This function plots the Hasse diagram of a surmise relation. #
#                                                              #
################################################################

hasse<-function(imp, items){
struct<-relation(domain = list(1:items,1:items),graph = imp)

#computation of parallel items
parallel<-list()
for(i in 1:items){
for(j in i:items){
if(relation_incidence(struct)[i,j] ==1 && relation_incidence(struct)[j,i] ==1){
if(length(parallel) == 0){parallel[[i]]<-set(i,j)} 
if(length(parallel) > 0 && length(parallel) < i && sum(sapply(parallel, function(x) is.element(i,x))) == 0){parallel[[i]]<-set(i,j)}
if(is.null(parallel[i][[1]])){}else{
if(sum(sapply(parallel[[i]], function(x) is.element(i,x))) > 0){parallel[[i]]<-set_union(parallel[[i]], set(i,j))}
}
}
}
}
parallel<-parallel[!sapply(parallel, is.null)]

#collapsing of parallel items
if(length(parallel) > 0){
pardrop<-0
for(i in 1:length(parallel)){
if(pardrop[1] == 0){pardrop<-sapply(parallel[[i]],invisible)[2:length(parallel[[i]])]}
else{
pardrop<-c(pardrop,sapply(parallel[[i]],invisible)[2:length(parallel[[i]])])
}
}

nparitems<-1:items
nparitems<-nparitems[-pardrop]
struct<-relation(domain = list(nparitems, nparitems), incidence = relation_incidence(struct)[-pardrop, -pardrop])
}

#plotting
plot(struct)

#returning a list of parallel items
return(parallel)
}
