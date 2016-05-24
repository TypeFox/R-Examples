DEEPR <-
function (group_1,group_2,perm_number = 999) {

n_group_1<-nrow(group_1)
group_1_model<-dirmult(group_1,trace=FALSE)
n_group_2<-nrow(group_2)
group_2_model<-dirmult(group_2,trace=FALSE)
group_0<-rbind(group_1,group_2)

obs_teststat<-group_1_model$loglik+group_2_model$loglik

perm_function<-function(group_0){
perm_group_0<-group_0[sample(nrow(group_0)),]
perm_group_1<-head(perm_group_0,n_group_1)
perm_group_1_model<-dirmult(perm_group_1,trace=FALSE)
perm_group_2<-tail(perm_group_0,n_group_2)
perm_group_2_model<-dirmult(perm_group_2,trace=FALSE)
perm_teststat<-perm_group_1_model$loglik+perm_group_2_model$loglik
}
perm_data<-replicate(perm_number,perm_function(group_0))
total_teststat<-c(perm_data,obs_teststat)
p_value<-mean(total_teststat >= obs_teststat)
list(p_value = p_value, group_1_pi = group_1_model$pi, group_2_pi = group_2_model$pi, group_1_alphas = group_1_model$gamma, group_2_alphas = group_2_model$gamma, group_1_theta = group_1_model$theta, group_2_theta = group_2_model$theta)
}
