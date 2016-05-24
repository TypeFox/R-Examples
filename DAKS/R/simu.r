###################
# data simulation #
###################

#########################################################
#                                                       #
# This function simulates a dataset using a basic local #
# independence model. The number of items, the sample   # 
# size, and two parameters for the careless error and   #
# lucky guess probabilities can be set explicitly. The  #
# underlying combinatorial structure can either be      #
# specified manually or is generated randomly.          #
#                                                       #
#########################################################

simu<-function(items, size, ce, lg, imp = NULL, delta){

R<-set()

if(is.null(imp)){
#computation of transitive relations
for(i in 1:items){
for(j in 1:items){
if(i != j && delta > runif(1,0,1)){R<-set_union(R, set(tuple(i,j)))}
if(i == j) {R<-set_union(R, set(tuple(i,j)))}
}
}

R_2<-relation_incidence(transitive_closure(relation(domain = list(1:items,1:items), graph = R)))

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

#implications

imp<-set()
for(i in 1:items){
for(j in 1:items){
if(i != j && set_is_subset(base_list[[i]], base_list[[j]])){imp<-set_union(imp,set(tuple(i,j)))}
}
}

#for specified imp
}else{
#Patterns
P<-imp2state(imp, items)
}

#simulating the dataset
sim<-matrix(ncol = items, nrow = size)

for(i in 1:size){
sim[i,]<-P[sample(1:nrow(P), 1),]
for(j in 1:items){
if(sim[i,j] == 1 && runif(1,0,1) < ce) {sim[i,j]<-0}
if(sim[i,j] == 0 && runif(1,0,1) < lg) {sim[i,j]<-1}
}
}

list(dataset = sim, implications = imp, states = P)
}
