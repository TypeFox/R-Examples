library(MASS)
fit_dist_gamma <-
function(dist){
l<-length(dist)
j<-1
data_gamma<-c()
fit_gamma<-fitdistr(dist, "gamma")

while (j<=l){
gamma<-round(rgamma(1, fit_gamma$estimate[1], rate=fit_gamma$estimate[2]))
if(gamma>0){
	data_gamma[j]<-gamma
	j<-j+1
}
}
chi_gamma<-chisq.test(dist, data_gamma)
list(data_gamma = data_gamma, fit_gamma = fit_gamma, chi_gamma = chi_gamma)
}
