library(MASS)
fit_dist_pois <-
function(dist){
l<-length(dist)
fit_pois<-fitdistr(dist, "poisson")
data_pois<-rpois(l, fit_pois$estimate)
chi_pois<-chisq.test(dist,data_pois)
list(data_pois = data_pois, fit_pois = fit_pois, chi_pois = chi_pois)
}
