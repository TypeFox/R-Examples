#################################
# estimated asymptotic variance #
#################################

###################################################  
#                                                 #
# This function computes estimated asymptotic     #
# variances of the maximum likelihood estimators  #
# diff from data, assuming a multinomial          #
# probability distribution on the set of all      #
# response patterns.                              #
#                                                 #
###################################################


variance<-function(dataset, imp, v){
if(length(imp) == 0){
stop("Number of implications must be greater than zero.\n")
}

if(v != 1 && v != 2){
stop("IITA version must be specified")
}

items<-ncol(dataset)

#Number of times a pattern occurs 
pat<-matrix(0,ncol = ncol(dataset), nrow = 2^ncol(dataset))
for(j in 1:ncol(pat)){
pat[,j]<-c(rep(0,2^(ncol(dataset)-j)), rep(1,2^(ncol(dataset)-j)))
}
patterns<-pattern(dataset, P = pat)$states

#relative frequencies
rho_sum<-vector(length = items)
for(i in 1:items){
rho_sum[i]<-sum(patterns[which(patterns[,i] == 1),items+1]) / nrow(dataset)
}

rho_sum_counter<-matrix(0, ncol = items, nrow = items)
for(i in 1:items){
for(j in 1:items){
if(i !=j){rho_sum_counter[i,j]<-sum(patterns[which(patterns[,i] == 0 & patterns[,j] == 1), items+1]) / nrow(dataset)}
}
}

#expected fisher information

exp_fish<-matrix(0, ncol = 2^items -1, nrow = 2^items -1)

for(i in 2:2^items){
for(j in 2:2^items){
if(i == j){exp_fish[i-1,j-1]<-(patterns[i,ncol(patterns)] / nrow(dataset)) * (1 - patterns[i,items+1]/nrow(dataset))}
if(i != j){exp_fish[i-1,j-1]<- (-1) * (patterns[i,ncol(patterns)] / nrow(dataset)) * (patterns[j,items+1] / nrow(dataset))}
}
}

#error

#original and corrected
error<-0
if(v == 2){
for(i in imp){
error<-error + ((rho_sum_counter[as.integer(i[1]), as.integer(i[2])]) * ncol(dataset) / sum(dataset[,as.integer(i[2])])) 
}
error<-error / length(imp)
}

#minimized corrected
if(v == 1){
x<-rep(0,4)
for(i in 1:items){
for(j in 1:items){
if(is.element(set(tuple(i,j)), imp) == TRUE && i != j){
x[2]<-x[2]-2*rho_sum_counter[i,j] * rho_sum[j]
x[4]<-x[4]+2 * rho_sum[j]^2
}
if(is.element(set(tuple(i,j)), imp) == FALSE && is.element(set(tuple(j,i)), imp) == TRUE && i != j){
x[1]<-x[1]-2*rho_sum_counter[i,j]*rho_sum[i] + 2 * rho_sum[i] * rho_sum[j] - 2 * rho_sum[i]^2  
x[3]<-x[3]+2*rho_sum[i]^2 
}
}
}
error<- -(x[1] + x[2]) / (x[3] + x[4])
}

#gamma derivative

# original and corrected
if(v == 2){
gamma_deriv<-rep(0,2^items -1)
for(i in 2:2^items){
for(j in imp){
if(patterns[i,as.integer(j[1])] == 0 && patterns[i,as.integer(j[2])] == 1){
gamma_deriv[i-1]<-gamma_deriv[i-1] + (rho_sum[as.integer(j[2])] - rho_sum_counter[as.integer(j[1]), as.integer(j[2])])/(rho_sum[as.integer(j[2])]^2)
}else{
if(patterns[i,as.integer(j[2])] ==1)
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
if(patterns[i,k] == 0 && patterns[i,h] == 1){
tmp1<-tmp1 + (-2) * rho_sum[k]
}
if(patterns[i,k] == 1){
tmp1<-tmp1 + (-2) * rho_sum_counter[k,h] + 2 * rho_sum[h] + (-4) * rho_sum[k]
tmp2<-tmp2 + 4 * rho_sum[k]
}
if(patterns[i,h] == 1){
tmp1<-tmp1 + 2 * rho_sum[k]
}
}
if(is.element(set(tuple(h,k)), imp) && k !=h){
if(patterns[i,k] == 0 && patterns[i,h] == 1){
tmp1<-tmp1 + (-2) * rho_sum[h]
}
if(patterns[i,h] == 1){
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
if(patterns[i,k] == 0 && patterns[i,h] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] * error) * (1 - error - (rho_sum[h])*gamma_deriv[i-1])
}else{
if(patterns[i, h] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] * error) * (-error - rho_sum[h] * gamma_deriv[i-1] )
}else{
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] * error) * (-rho_sum[h] * gamma_deriv[i-1])
}
}
}
if(is.element(set(tuple(k,h)), imp) == FALSE && is.element(set(tuple(h,k)), imp) && k !=h){
if(patterns[i,k] == 1){
if(patterns[i,h] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error) * ( - error - rho_sum[k] * gamma_deriv[i-1])
}else{
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error) * (1 - error - rho_sum[k] * gamma_deriv[i-1])
}
}else{
if(patterns[i,k] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error) * (-rho_sum[k] * gamma_deriv[i-1])
}else{
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - rho_sum[h] + rho_sum[k] - rho_sum[k] * error) * (-rho_sum[k] * gamma_deriv[i-1])
}
}
}
if(is.element(set(tuple(k,h)), imp) == FALSE && is.element(set(tuple(h,k)), imp) == FALSE && k !=h){
if(patterns[i,k] == 0 && patterns[i,h] == 1){
grad[i-1]<-grad[i-1] + 2 * (rho_sum_counter[k,h] - (1-rho_sum[k]) * rho_sum[h]) * rho_sum[k] 
}else{
if(patterns[i,h] == 1){
if(patterns[i,k] == 1){
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
