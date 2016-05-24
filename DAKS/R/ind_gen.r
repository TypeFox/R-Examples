########################
# inductive generation #
########################

#######################################
#                                     #
# This function generates inductively #
# a set of competing quasi orders.    #
#                                     #
#######################################

ind_gen<-function(b){
m<-nrow(b)
#set of all pairs with a maximum of k-1 counterexamples
S<-list()

#constructed relation for a maximum of k-1 counterexamples
A<-list()

#set of non-transitive triples
M<-list()
M[[1]]<-set()
S[[1]]<-set()
for(i in 1:m){
for(j in 1:m){
if(b[i,j] == min(b) & i != j) {S[[1]]<-set_union(S[[1]],set(tuple(i,j)))}
}
}

A[[1]]<-S[[1]]

#inductive gneration process
elements<-sort(b)[!duplicated(sort(b))]
if(is.element(0,elements)){elements<-elements[2:(length(elements))]}

k<-2

for(elem in elements){
S[[k]]<-set()
A[[k]]<-set()
M[[k]]<-set()
#building of S
for(i in 1:m){
for(j in 1:m){
if(b[i,j] <= elem && i !=j && is.element(set(tuple(i,j)), A[[k-1]]) == FALSE){S[[k]]<-set_union(S[[k]], set(tuple(i,j)))}
}
}
#transitivity test
if(set_is_empty(S[[k]]) == FALSE){
M[[k]]<-S[[k]]
brake_test<-1
while(brake_test != 0){
brake<-M[[k]]
for(i in M[[k]]){
for(h in 1:m){
if(h != as.integer(i)[1] && h != as.integer(i)[2] && is.element(set(tuple(as.integer(i)[2],h)), set_union(A[[k-1]], M[[k]])) == TRUE && is.element(set(tuple(as.integer(i)[1],h)), set_union(A[[k-1]], M[[k]])) == FALSE){M[[k]]<-M[[k]] - set(i)}
if(h != as.integer(i)[1] && h != as.integer(i)[2] && is.element(set(tuple(h,as.integer(i)[1])), set_union(A[[k-1]], M[[k]])) == TRUE && is.element(set(tuple(h,as.integer(i)[2])), set_union(A[[k-1]], M[[k]])) == FALSE){M[[k]]<-M[[k]] - set(i)}
}
}
if(brake == M[[k]]){brake_test<-0}
}
A[[k]]<-set_union(A[[k-1]], (M[[k]]))
}
k<-k+1
}

#deletion of empty and duplicated quasi orders
A<-A[(!duplicated(A))]
A<-A[!set_is_empty(A)]
return(A)
}