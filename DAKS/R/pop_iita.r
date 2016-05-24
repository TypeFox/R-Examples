##################### 
# Population values #
#####################

###################################################### 
#                                                    #
# This function can be used to perform one of        #
# the three inductive item tree analysis algorithms  #
# in population quantities selectively.              #
#                                                    # 
######################################################

pop_iita<-function(imp, ce, lg, items, dataset = NULL, A = NULL, v){

if(v != 1 && v != 2 && v !=3){
stop("IITA version must be specified")
}

if(is.null(dataset) == FALSE && is.null(A) == FALSE){
stop("Either the dataset or the selection set can be specified")
}

p_pop<-vector(length = items)
b_pop<-matrix(0,ncol = items, nrow = items)
bs_pop_alt<-matrix(0,ncol = items, nrow = items)
bs_pop_neu<-matrix(0,ncol = items, nrow = items)
bs_pop_num<-matrix(0,ncol = items, nrow = items)
error_pop_theo<-matrix(ncol = items, nrow = items)

#Population matrix
pop_matrix<-matrix(0,nrow = 2^items, ncol = items+1)
for(i in 1:items){
pop_matrix[,i]<-c(rep(0, 2^(items-i)), rep(1, 2^(items-i)))
}

#Probabilities of all patterns

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

#computation of population matrix
r_pop<-matrix(1, ncol = length(G[[items+1]]), nrow = 2^items)

for(i in 1:2^items){
for(j in 1:(length(G[[items+1]]))){
for(k in 1:items){
if(pop_matrix[i,k] == 0 && P[j,k] == 1){r_pop[i,j]<-r_pop[i,j] * ce}
if(pop_matrix[i,k] == 1 && P[j,k] == 1){r_pop[i,j]<-r_pop[i,j] * (1-ce)}
if(pop_matrix[i,k] == 1 && P[j,k] == 0){r_pop[i,j]<-r_pop[i,j] * lg}
if(pop_matrix[i,k] == 0 && P[j,k] == 0){r_pop[i,j]<-r_pop[i,j] * (1-lg)}
}
pop_matrix[i,items+1]<-pop_matrix[i,items+1] + r_pop[i,j]/length(G[[items+1]])
}
}

#item probabilities and b_ij
for(i in 1:items){
p_pop[i]<-sum(pop_matrix[pop_matrix[,i] == 1,items+1 ])
for(j in 1:items){
tmp1<-which(pop_matrix[,i] == 0)
tmp2<-which(pop_matrix[,j] == 1)
if(i != j){b_pop[i,j]<-sum(pop_matrix[c(tmp1,tmp2)[duplicated(c(tmp1, tmp2))] ,items+1])}
}
}

if(is.null(dataset) && is.null(A)){
#inductive generation process on population values
S<-list()
A<-list()
M<-list()

S[[1]]<-set()
A[[1]]<-set()
M[[1]]<-set()

elements<-sort(b_pop)[!duplicated(sort(b_pop))]
elements<-elements[2:(length(elements))]

k<-2

for(elem in elements){
S[[k]]<-set()
A[[k]]<-set()
M[[k]]<-set()
for(i in 1:items){
for(j in 1:items){
if(b_pop[i,j] <= elem && i !=j && is.element(set(tuple(i,j)), A[[k-1]]) == FALSE){S[[k]]<-set_union(S[[k]], set(tuple(i,j)))}
}
}
#transitivity test
if(set_is_empty(S[[k]])){A[[k]]<-A[[k-1]]}
if(set_is_empty(S[[k]]) == FALSE){
M[[k]]<-S[[k]]
brake_test<-1
while(brake_test != 0){
brake<-M[[k]]
for(i in M[[k]]){
for(h in 1:items){
if(h != as.integer(i)[1] && h != as.integer(i)[2] && is.element(set(tuple(as.integer(i)[2],h)), set_union(A[[k-1]], M[[k]])) == TRUE && is.element(set(tuple(as.integer(i)[1],h)), set_union(A[[k-1]], M[[k]])) == FALSE){M[[k]]<-set_intersection(M[[k]], set_symdiff(M[[k]],set(i)))}
if(h != as.integer(i)[1] && h != as.integer(i)[2] && is.element(set(tuple(h,as.integer(i)[1])), set_union(A[[k-1]], M[[k]])) == TRUE && is.element(set(tuple(h,as.integer(i)[2])), set_union(A[[k-1]], M[[k]])) == FALSE){M[[k]]<-set_intersection(M[[k]], set_symdiff(M[[k]],set(i)))}
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
}else{
if(is.null(dataset) == FALSE){
A<-ind_gen(ob_counter(dataset))
}
}

error_pop<-vector(length = (length(A)))
error_pop_num<-vector(length = (length(A)))

#Gamma_ij
for(i in 1:items){
error_pop_theo[,i]<-b_pop[,i] / p_pop[i] 
}

if(v == 3 | v == 2){

all_imp<-set()

for(i in 1:(items-1)){
for(j in (i+1):items){
all_imp<-set_union(all_imp, set(tuple(i,j), tuple(j,i)))
}
}

#Gamma_L
for(k in 1:length(A)){
for(i in A[[k]]){
error_pop[k]<-error_pop[k] + ((b_pop[as.integer(i[1]), as.integer(i[2])]) / sum(p_pop[as.integer(i[2])])) 
}
if(set_is_empty(A[[k]])){error_pop[k]<-NA}
if(set_is_empty(A[[k]]) == FALSE) {error_pop[k]<-error_pop[k] / length(A[[k]])}
}
}

if(v == 1){
#Gamma_min
for(k in 1:length(A)){
x<-rep(0,4)
for(i in 1:items){
for(j in 1:items){
if(is.element(set(tuple(i,j)), A[[k]]) == TRUE && i != j){
x[2]<-x[2]-2*b_pop[i,j] * p_pop[j]
x[4]<-x[4]+2 * p_pop[j]^2
}
if(is.element(set(tuple(i,j)), A[[k]]) == FALSE && is.element(set(tuple(j,i)), A[[k]]) == TRUE && i != j){
x[1]<-x[1]-2*b_pop[i,j]*p_pop[i] + 2 * p_pop[i] * p_pop[j] - 2 * p_pop[i]^2  
x[3]<-x[3]+2*p_pop[i]^2 
}
}
}
error_pop_num[k]<- -(x[1] + x[2]) / (x[3] + x[4])
}
}

#bs_ij and diff

if(v == 3){
#original
diff_value_pop_alt<-vector(length = length(A))
for(k in 1:length(A)){
for(i in 1:items){
for(j in 1:items){
if(is.element(set(tuple(i,j)), A[[k]]) == TRUE && i != j){bs_pop_alt[i,j]<-error_pop[k] * p_pop[j]}
if(is.element(set(tuple(i,j)), A[[k]]) == FALSE && i != j){bs_pop_alt[i,j]<-(1-p_pop[i]) * p_pop[j] * (1-error_pop[k])}
}
}
if(set_is_empty(A[[k]])){diff_value_pop_alt[k]<-NA}
if(set_is_empty(A[[k]]) == FALSE) {diff_value_pop_alt[k]<-sum((b_pop - bs_pop_alt)^2) / (items^2 - items)}
}
obj<-list(pop.diff = diff_value_pop_alt, pop.matrix = pop_matrix, error.pop = error_pop, selection.set = A, v = v)
class(obj)<-"popiita"	
return(obj)
}

if(v == 2){
#corrected
diff_value_pop_neu<-vector(length = length(A))
for(k in 1:length(A)){
for(i in 1:items){
for(j in 1:items){
if(is.element(set(tuple(i,j)), A[[k]]) == TRUE && i != j){bs_pop_neu[i,j]<-error_pop[k] * p_pop[j]}
if(is.element(set(tuple(i,j)), A[[k]]) == FALSE && is.element(set(tuple(j,i)), A[[k]]) == FALSE && i != j){bs_pop_neu[i,j]<-(1-p_pop[i]) * p_pop[j]}
if(is.element(set(tuple(i,j)), A[[k]]) == FALSE && is.element(set(tuple(j,i)), A[[k]]) == TRUE && i != j){bs_pop_neu[i,j]<-p_pop[j] - p_pop[i] + p_pop[i] * error_pop[k]}
}
}
if(set_is_empty(A[[k]])){diff_value_pop_neu[k]<-NA}
if(set_is_empty(A[[k]]) == FALSE) {diff_value_pop_neu[k]<-sum((b_pop - bs_pop_neu)^2) / (items^2 - items)}
}
obj<-list(pop.diff = diff_value_pop_neu, pop.matrix = pop_matrix, error.pop = error_pop, selection.set = A, v = v)
class(obj)<-"popiita"	
return(obj)
}

if(v == 1){
#minimized corrected
diff_value_pop_num<-vector(length = length(A))
for(k in 1:length(A)){
for(i in 1:items){
for(j in 1:items){
if(is.element(set(tuple(i,j)), A[[k]]) == TRUE && i != j){bs_pop_num[i,j]<-error_pop_num[k] * p_pop[j]}
if(is.element(set(tuple(i,j)), A[[k]]) == FALSE && is.element(set(tuple(j,i)), A[[k]]) == FALSE && i != j){bs_pop_num[i,j]<-(1-p_pop[i]) * p_pop[j]}
if(is.element(set(tuple(i,j)), A[[k]]) == FALSE && is.element(set(tuple(j,i)), A[[k]]) == TRUE && i != j){bs_pop_num[i,j]<-p_pop[j] - p_pop[i] + p_pop[i] * error_pop_num[k]}
}
}
if(set_is_empty(A[[k]])){diff_value_pop_num[k]<-NA}
if(set_is_empty(A[[k]]) == FALSE) {diff_value_pop_num[k]<-sum((b_pop - bs_pop_num)^2) / (items^2 - items)}
}
obj<-list(pop.diff = diff_value_pop_num, pop.matrix = pop_matrix, error.pop = error_pop_num, selection.set = A, v = v)
class(obj)<-"popiita"	
return(obj)
}
}
