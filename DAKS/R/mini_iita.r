############################ 
# minimized corrected IITA #
############################

##################################################
#                                                # 
# This function performs the minimized corrected # 
# inductive item tree analysis procedure and     #
# returns the corresponding diff values.         #
#                                                # 
##################################################

mini_iita<-function(dataset, A){
b<-ob_counter(dataset)
m<-ncol(dataset)
n<-nrow(dataset)

bs_num<-list()
for(i in 1:length(A)){
bs_num[[i]]<-matrix(0,ncol = m, nrow = m)
}

p<-rep(0,m)
for(i in 1:m){p[i]<-sum(dataset[,i])}

error_num<-rep(0,length(A))
diff_value_num<-rep(0,length(A))

#computation of error rate
for(k in 1:length(A)){
x<-rep(0,4)
for(i in 1:m){
for(j in 1:m){
if(is.element(set(tuple(i,j)), A[[k]]) == TRUE && i != j){
x[2]<-x[2]-2*b[i,j] * p[j]
x[4]<-x[4]+2 * p[j]^2
}
if(is.element(set(tuple(i,j)), A[[k]]) == FALSE && is.element(set(tuple(j,i)), A[[k]]) == TRUE && i != j){
x[1]<-x[1]-2*b[i,j]*p[i] + 2 * p[i] * p[j] - 2 * p[i]^2  
x[3]<-x[3]+2*p[i]^2 
}
}
}
error_num[k]<- -(x[1] + x[2]) / (x[3] + x[4])
}

#computation of diff values
all_imp<-set()

for(i in 1:(m-1)){
for(j in (i+1):m){
all_imp<-set_union(all_imp, set(tuple(i,j), tuple(j,i)))
}
}

for(k in 1:length(A)){
if(set_is_empty(A[[k]])){diff_value_num[k]<-NA}
else{
for(i in all_imp){
if(is.element(set(i), A[[k]])) {bs_num[[k]][as.integer(i[1]),as.integer(i[2])]<-error_num[k] * sum(dataset[,as.integer(i[2])])}
if(is.element(set(i), A[[k]]) == FALSE && is.element(set(tuple(as.integer(i[2]),as.integer(i[1]))), A[[k]]) == FALSE){bs_num[[k]][as.integer(i[1]),as.integer(i[2])]<-(1- sum(dataset[,as.integer(i[1])]) / n) * sum(dataset[,as.integer(i[2])])}
if(is.element(set(i), A[[k]]) == FALSE && is.element(set(tuple(as.integer(i[2]),as.integer(i[1]))), A[[k]]) == TRUE){bs_num[[k]][as.integer(i[1]),as.integer(i[2])]<-sum(dataset[,as.integer(i[2])]) - sum(dataset[,as.integer(i[1])]) + sum(dataset[,as.integer(i[1])]) * error_num[k]}
}
diff_value_num[k]<-sum((b - bs_num[[k]])^2) / (m^2 - m)
}
}

return(list(diff.value = diff_value_num, error.rate = error_num))
}
