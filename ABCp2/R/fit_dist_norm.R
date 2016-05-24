library(MASS)
fit_dist_norm <-
function(dist){
l<-length(dist)
j<-1
data_norm<-c()
fit_norm<-fitdistr(dist, "normal")

while (j<=l){
norm<-round(rnorm(1, mean=fit_norm$estimate[1], sd=fit_norm$estimate[2]))
if(norm>0){
	data_norm[j]<-norm
	j<-j+1
}
}
chi_norm<-chisq.test(dist,data_norm)
list(data_norm = data_norm, fit_norm = fit_norm, chi_norm = chi_norm)
}
