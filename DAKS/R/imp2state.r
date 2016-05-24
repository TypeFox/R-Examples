##################
# tranformations #
##################

##################################################
#                                                #
# This function transforms a set of implications #
# to the corresponding set of knowledge states.  #
#                                                #
##################################################

imp2state<-function(imp, items){
R_2<-matrix(1, ncol = items, nrow = items)
for(i in 1:items){
for(j in 1:items){
if(!is.element(set(tuple(i,j)), imp) && i != j){R_2[j,i]<-0}
}
}

#Base
base<-list()

for(i in 1:items){
tmp<-vector()
for(j in 1:items){
if(R_2[i,j] == 1){tmp<-c(tmp, j)}
}
base[[i]]<-sort(tmp)
}

base_list<-list()
for(i in 1:items){
base_list[[i]]<-set()
for(j in 1:length(base[[i]]))
base_list[[i]]<-set_union(base_list[[i]], set(base[[i]][j]))
}

#span of base
G<-list()
G[[1]]<-set(set())
G[[2]]<-set()
for(i in 1:length(base[[1]])){G[[2]]<-set_union(G[[2]], base[[1]][i])}
G[[2]]<-set(set(), G[[2]])

for(i in 2:items){
H<-set(set())
for(j in G[[i]]){
check<-0
if(set_is_subset(base_list[[i]], j) == FALSE){
for(d in 1:i){
if(set_is_subset(base_list[[d]], set(j, base_list[[i]])) == TRUE){
if(set_is_subset(base_list[[d]], j)){
H<-set_union(H,set(set_union(j,base_list[[i]])))
}
}
if(set_is_subset(base_list[[d]], set(j, base_list[[i]])) == FALSE){
H<-set_union(H,set(set_union(j,base_list[[i]]))) 
}
}
}
}
G[[i+1]]<-set_union(G[[i]], H)
}

#Patterns

P<-matrix(0, ncol = items, nrow = length(G[[items+1]]))
i<-1

for(k in (G[[items+1]])){
for(j in 1:items){
if(is.element(j, k)){P[i,j]<-1}
}
i<-i+1
}

return(P)
}