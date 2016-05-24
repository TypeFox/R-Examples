#######################
# population variance #
#######################

########################################################
#                                                      #
# This function  computes the population asymptotic    #
# variances of the maximum likelihood estimators diff, #
# assuming a multinomial probability distribution on   #
# the set of all response patterns.                    #
#                                                      # 
########################################################

pop_variance<-function(pop_matrix, imp, error_pop, v){
if(length(imp) == 0){
stop("Number of implications must be greater than zero.\n")
}

if(v != 1 && v != 2){
stop("IITA version must be specified")
}

items<-ncol(pop_matrix)-1

#expected fisher information
exp_fish<-matrix(0, ncol = 2^items -1, nrow = 2^items -1)

for(i in 2:2^items){
for(j in 2:2^items){
if(i == j){exp_fish[i-1,j-1]<-pop_matrix[i,items+1] * (1-pop_matrix[i,items+1])}
if(i != j){exp_fish[i-1,j-1]<- (-1) * pop_matrix[i,items+1] * pop_matrix[j,items+1]}
}
}

#Sum of the rho, to avoid computaion every single time
rho_sum<-vector(length = items)
for(i in 1:items){
rho_sum[i]<-sum(pop_matrix[which(pop_matrix[,i] == 1),items+1])
}

rho_sum_counter<-matrix(0, ncol = items, nrow = items)
for(i in 1:items){
for(j in 1:items){
if(i !=j){rho_sum_counter[i,j]<-sum(pop_matrix[which(pop_matrix[,i] == 0 & pop_matrix[,j] == 1), items+1])}
}
}

#gamma derivative

#corrected and original
 
if(v == 2){
gamma_deriv<-rep(0,2^items -1)
for(i in 2:2^items){
for(j in imp){
if(pop_matrix[i,as.integer(j[1])] == 0 && pop_matrix[i,as.integer(j[2])] == 1){
gamma_deriv[i-1]<-gamma_deriv[i-1] + (rho_sum[as.integer(j[2])] - rho_sum_counter[as.integer(j[1]), as.integer(j[2])])/(rho_sum[as.integer(j[2])]^2)
}else{
if(pop_matrix[i,as.integer(j[2])] ==1)
gamma_deriv[i-1]<-gamma_deriv[i-1] + rho_sum_counter[as.integer(j[1]), as.integer(j[2])]/(rho_sum[as.integer(j[2])]^2)
}
}
}

gamma_deriv<-gamma_deriv/length(imp)
}

#minimized corrected
if(v == 1){
gamma_deriv<-rep(0, 2^items -1)
x<-rep(0,4)

for(k in 1:items){
for(h in 1:items){
if(is.element(set(tuple(k,h)), imp) == FALSE && is.element(set(tuple(h,k)), imp) && k !=h){
x[1]<- x[1] + (-2) * rho_sum_counter[k,h] * rho_sum[k] + 2 * rho_sum[k] * rho_sum[h] - 2 * (rho_sum[k])^2
x[3]<- x[3] + 2 * (rho_sum[k])^2
}
if(is.element(set(tuple(k,h)), imp) && k !=h){
x[2]<- x[2] + (-2) * rho_sum_counter[k,h] * rho_sum[h]
x[4]<- x[4] + 2 * (rho_sum[h])^2
}
}
}

tmp1<-0
tmp2<-0
for(i in 2:2^items){
for(k in 1:items){
for(h in 1:items){
if(is.element(set(tuple(k,h)), imp) == FALSE && is.element(set(tuple(h,k)), imp) && k !=h){
if(pop_matrix[i,k] == 0 && pop_matrix[i,h] == 1){
tmp1<-tmp1 + (-2) * rho_sum[k]
}
if(pop_matrix[i,k] == 1){
tmp1<-tmp1 + (-2) * rho_sum_counter[k,h] + 2 * rho_sum[h] + (-4) * rho_sum[k]
tmp2<-tmp2 + 4 * rho_sum[k]
}
if(pop_matrix[i,h] == 1){
tmp1<-tmp1 + 2 * rho_sum[k]
}
}
if(is.element(set(tuple(h,k)), imp) && k !=h){
if(pop_matrix[i,k] == 0 && pop_matrix[i,h] == 1){
tmp1<-tmp1 + (-2) * rho_sum[h]
}
if(pop_matrix[i,h] == 1){
tmp1<-tmp1 + (-2) * rho_sum_counter[k,h]
tmp2<-tmp2 + 4 * (rho_sum[h])
}
}
}
}
gamma_deriv[i-1]<- (-1) * (tmp1 * (x[3] + x[4]) - tmp2 * (x[1] + x[2])) / ((x[3] + x[4])^2 ) 
}
}


#gradient of diff for corrected and minimized corrected
grad<-rep(0, 2^items -1)
for(i in 2:2^items){
for(k in 1:items){
for(h in 1:items){
if(is.element(set(tuple(k,h)), imp) && k !=h){
if(pop_matrix[i,k] == 0 && pop_matrix[i,h] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] * error_pop) * (1 - error_pop - (rho_sum[h])*gamma_deriv[i-1])
}else{
if(pop_matrix[i, h] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] * error_pop) * (-error_pop - rho_sum[h] * gamma_deriv[i-1] )
}else{
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] * error_pop) * (-rho_sum[h] * gamma_deriv[i-1])
}
}
}
if(is.element(set(tuple(k,h)), imp) == FALSE && is.element(set(tuple(h,k)), imp) && k !=h){
if(pop_matrix[i,k] == 0 && pop_matrix[i,h] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error_pop) * (-rho_sum[k] * gamma_deriv[i-1])
}else{
if(pop_matrix[i,h] == 1){
if(pop_matrix[i,k] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error_pop) * (-1 + 1 - error_pop - rho_sum[k] * gamma_deriv[i-1])
}else{
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error_pop) * (-1 - rho_sum[k] * gamma_deriv[i-1])
}
}else{
if(pop_matrix[i,k] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error_pop) * (1 - error_pop - rho_sum[k] * gamma_deriv[i-1])
}else{
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error_pop) * (-rho_sum[k] * gamma_deriv[i-1])
}
}
}
}
if(is.element(set(tuple(k,h)), imp) == FALSE && is.element(set(tuple(h,k)), imp) == FALSE && k !=h){
if(pop_matrix[i,k] == 0 && pop_matrix[i,h] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - (1-rho_sum[k]) * rho_sum[h]) * rho_sum[k] 
}else{
if(pop_matrix[i,h] == 1){
if(pop_matrix[i,k] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - (1-rho_sum[k]) * rho_sum[h]) * (rho_sum[h] - 1 + rho_sum[k]) 
}else{
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - (1-rho_sum[k]) * rho_sum[h]) * rho_sum[h]
}
}
}
}
}
}
}
grad<-grad / (items * (items-1))

#final computation
variance<- grad%*%exp_fish%*%grad
variance<-as.vector(variance)
return(variance)
}
