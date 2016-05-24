##################
# tranformations #
##################

####################################################
#                                                  #
# This function transforms a set of knowledge      #
# states to the corresponding set of implications. #
#                                                  #
####################################################

state2imp<-function(P){

#Base
base_list<-list()
for(i in 1:ncol(P)){
base_list[[i]]<-set(i)
tmp<-P[which(P[,i] == 1),]
for(j in (1:ncol(P))[-i]){
if(length(which(P[,i] == 1)) == 1){if(sum(tmp[j]) == 1){base_list[[i]]<-set_union(base_list[[i]], set(j))}}
if(length(which(P[,i] == 1)) > 1){if(sum(tmp[,j]) == nrow(tmp)){base_list[[i]]<-set_union(base_list[[i]], set(j))}}
}
}

imp<-set()
for(i in 1:ncol(P)){
for(j in 1:ncol(P)){
if(i != j && set_is_subset(base_list[[i]], base_list[[j]])){imp<-set_union(imp,set(tuple(i,j)))}
}
}

return(imp)
}