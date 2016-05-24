##############################################
#                                            #
# This function performs a Z-test            #
# for testing hypotheses on the diff values  #
#                                            #
##############################################

z_test<-function(dataset, imp, imp_alt = NULL, alternative = c("two.sided", "less", "greater"), mu = 0, conf.level = 0.95, v){
if(is.null(imp_alt)){
if(v == 1){iita_imp<-mini_iita(dataset, list(imp))}
if(v == 2){iita_imp<-corr_iita(dataset, list(imp))}	
#if(v == 3){var_imp<-orig_iita(dataset, list(imp))}		

var_imp<-variance(dataset, imp, v)

Z<-sqrt(nrow(dataset)) * (iita_imp$diff/nrow(dataset)^2 - mu) / sqrt(var_imp)

if(alternative[1] == "two.sided"){
out<-list(Z.value = Z, p.value = 2 * (min(pnorm(Z), 1-pnorm(Z))), conf = c(iita_imp$diff/nrow(dataset)^2 - qnorm(1-(1-conf.level)/2) * sqrt(var_imp/nrow(dataset)),iita_imp$diff/nrow(dataset)^2 + qnorm(1-(1-conf.level)/2) * sqrt(var_imp/nrow(dataset))),diff_value = iita_imp$diff/nrow(dataset)^2, alternative = alternative[1], mu = mu, conf.level = conf.level)
}
if(alternative[1] == "greater"){
out<-list(Z.value = Z, p.value = pnorm(Z),conf = c(iita_imp$diff/nrow(dataset)^2 - qnorm(conf.level) * sqrt(var_imp/nrow(dataset)),Inf), diff_value = iita_imp$diff/nrow(dataset)^2, alternative = alternative[1], mu = mu, conf.level = conf.level)
}
if(alternative[1] == "less"){
out<-list(Z.value = Z, p.value = 1-pnorm(Z),conf = c(iita_imp$diff/nrow(dataset)^2 - qnorm(conf.level) * sqrt(var_imp/nrow(dataset)),Inf), diff_value = iita_imp$diff/nrow(dataset)^2, alternative = alternative[1], mu = mu, conf.level = conf.level)
}
}

else{
if(v == 1){iita_imp<-mini_iita(dataset, list(imp,imp_alt))}
if(v == 2){iita_imp<-corr_iita(dataset, list(imp, imp_alt))}	
#if(v == 3){var_imp<-orig_iita(dataset, list(imp, imp_alt))}	

var_imp<-sapply(list(imp, imp_alt), function(x){variance(dataset,x,v)})

Z<-sqrt(nrow(dataset)) * ((iita_imp$diff[1] - iita_imp$diff[2])/nrow(dataset)^2 - mu) / sqrt(var_imp[1] + var_imp[2])

if(alternative[1] == "two.sided"){
out<-list(Z.value = Z, p.value = 2 * (min(pnorm(Z), 1-pnorm(Z))), conf = c((iita_imp$diff[1] - iita_imp$diff[2])/nrow(dataset)^2 - qnorm(1-(1-conf.level)/2) * sqrt(sum(var_imp)/nrow(dataset)),(iita_imp$diff[1] - iita_imp$diff[2])/nrow(dataset)^2 + qnorm(1-(1-conf.level)/2) * sqrt(sum(var_imp)/nrow(dataset))),diff_value = iita_imp$diff/nrow(dataset)^2, alternative = alternative[1], mu = mu, conf.level = conf.level)
}
if(alternative[1] == "greater"){
out<-list(Z.value = Z, p.value = pnorm(Z),conf = c((iita_imp$diff[1] - iita_imp$diff[2])/nrow(dataset)^2 - qnorm(conf.level) * sqrt(sum(var_imp)/nrow(dataset)),Inf), diff_value = iita_imp$diff/nrow(dataset)^2, alternative = alternative[1], mu = mu, conf.level = conf.level)	
}
if(alternative[1] == "less"){
out<-list(Z.value = Z, p.value = 1-pnorm(Z),conf = c((iita_imp$diff[1] - iita_imp$diff[2])/nrow(dataset)^2 - qnorm(conf.level) * sqrt(sum(var_imp)/nrow(dataset)),Inf), diff_value = iita_imp$diff/nrow(dataset)^2, alternative = alternative[1], mu = mu, conf.level = conf.level)
}
}
class(out)<-"ztest"
out
}
